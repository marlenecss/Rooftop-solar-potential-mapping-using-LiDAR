################################################################################
# LiDAR course project
# Rooftop Solar Potential Mapping from ALS using OpenData Bavaria
#
# Data: Laserdaten (ALS point clouds, classified) + Hausumringe (building footprints)
# Source: https://geodaten.bayern.de/opengeodata/OpenDataDetail.html?pn=laserdaten
# Download: 14.09.2026
# CRS: ETRS89 / UTM zone 32N (EPSG:25832) -- Bavaria ALS tiles have no CRS in the
#      LAS header, so it is assigned manually throughout this script.
#
# v2: adds a cast-shadow model (neighbouring buildings, chimneys, dormers) to the
#     irradiation loop, reorders the facet segmentation, and uses one colour
#     palette for all figures.
#
# author: Marlene Sehrbrock
# 09.10.2026
################################################################################


library(lidR)
library(sf)
library(terra)
library(mapview)
library(ggplot2)
library(solrad)
library(suncalc)
library(rayshader)

mapviewOptions(fgb = FALSE, basemaps = c("OpenStreetMap", "CartoDB.Positron"))

# working folder and tile -- change once here, used throughout
data_dir <- "C:/Users/marle/OneDrive - Universität Würzburg/Sem2/Lidar/"
tile_id  <- "651_5478"


################################################################################
# COLOUR PALETTES (used for every figure below)
# Facet classes: code = aspect class * 10 + slope class
#   aspect class: 0 = flat, 1 = N, 2 = E, 3 = S, 4 = W
#   slope class : 1 = <15 deg, 2 = >=15 deg
facet_lut <- data.frame(
  id    = c(1,      2,              11,       12,        21,       22,        31,       32,        41,       42),
  label = c("Flat", "Flat (steep)", "N <15°", "N ≥15°",  "E <15°", "E ≥15°",  "S <15°", "S ≥15°",  "W <15°", "W ≥15°"),
  col   = c("grey85", "grey60",     "#c6dbef", "#2171b5", "#c7e9c0", "#238b45", "#fdd0a2", "#d94801", "#e1c3f0", "#6a3d9a")
)

# Suitability classes: grey -> orange family taken from the palette above
# (#fd8d3c is the one in-between orange added to get four steps)
suit_cols <- c(Low = "grey85", Moderate = "#fdd0a2", Good = "#fd8d3c", Excellent = "#d94801")

# Continuous irradiation ramp, same colour family
irr_ramp <- colorRampPalette(c("grey85", "#fdd0a2", "#fd8d3c", "#d94801"))


################################################################################
# 1. CHOOSE AOI TILE FROM THE FULL CATALOG

ctg <- readLAScatalog(data_dir)
sf::st_crs(ctg) <- 25832

tiles <- st_as_sf(ctg)
tiles$filename <- basename(tiles$filename)
tiles <- tiles[, c("filename", "geometry")]

mapview(tiles, zcol = "filename", label = tiles$filename)


################################################################################
# 2. LOAD CHOSEN TILE
aoi <- readLAS(paste0(data_dir, tile_id, ".laz"))
projection(aoi) <- 25832   # Bavaria ALS = ETRS89 / UTM zone 32N, header has no CRS

las_check(aoi)
table(aoi@data$Classification) # 2 = ground, 6 = building, etc.
plot(aoi, color = "Classification") # 3D view in a separate window


################################################################################
# 3. BUILD DTM AND DSM
dtm <- rasterize_terrain(aoi, res = 0.5, algorithm = tin())

dsm <- rasterize_canopy(aoi, res = 0.5,
                        algorithm = pitfree(thresholds = c(0, 2, 5, 10, 15),
                                            max_edge = c(0, 1.5)))
crs(dsm) <- "EPSG:25832"


################################################################################
# 4. TILE CENTROID LAT/LON (for solar geometry)
aoi_bbox <- st_as_sfc(st_bbox(aoi))
st_crs(aoi_bbox) <- 25832

aoi_center_wgs84 <- st_transform(st_centroid(aoi_bbox), 4326)
coords <- st_coordinates(aoi_center_wgs84)
lon <- coords[1]
lat <- coords[2]
lon; lat # check


################################################################################
# 5. BUILDING FOOTPRINTS (Hausumringe)
buildings <- st_read(paste0(data_dir, "hausumringe.shp"))
buildings <- st_transform(buildings, 25832)

# clip to aoi
buildings_tile <- buildings[st_intersects(buildings, aoi_bbox, sparse = FALSE), ]
buildings_tile$build_id <- seq_len(nrow(buildings_tile))

plot(st_geometry(buildings_tile), col = facet_lut$col[3], border = "grey40")

st_write(buildings_tile, paste0(data_dir, "buildings_", tile_id, ".shp"), delete_layer = TRUE)


################################################################################
# 6. MASK DSM TO ROOF FOOTPRINTS
roof_dsm <- mask(dsm, vect(buildings_tile))
roof_dsm <- crop(roof_dsm, vect(buildings_tile))
plot(roof_dsm, col = irr_ramp(100), main = "Roof DSM (m)")


################################################################################
# 7. SLOPE AND ASPECT
slope  <- terrain(roof_dsm, v = "slope",  unit = "degrees")
aspect <- terrain(roof_dsm, v = "aspect", unit = "degrees")

# Light smoothing -- used for slope (flat-vs-pitched distinction, noise-tolerant)
# and for the irradiation calc. Aspect is kept UNsmoothed for facet segmentation,
# since smoothing blurs the sharp ridge-line boundary between opposite roof pitches.
slope_smooth  <- focal(slope,  w = 3, fun = "median", na.rm = TRUE)
aspect_smooth <- focal(aspect, w = 3, fun = "median", na.rm = TRUE)


################################################################################
# 8. SHADOW MODEL -- cast shadows from the FULL DSM (not just the roofs), so
#    neighbouring buildings, chimneys and dormers shade each other.
# Shadow grid resolution = DSM resolution * shade_fact.
# 2 -> 1 m grid, ~4x faster (good for testing). Use 1 for the final run.
shade_fact <- 2

dsm_s <- if (shade_fact > 1) aggregate(dsm, fact = shade_fact, fun = "max", na.rm = TRUE) else dsm
dsm_s <- subst(dsm_s, NA, min(values(dsm_s), na.rm = TRUE))   # ray_shade() cannot handle NA

dsm_s_mat <- raster_to_matrix(dsm_s, verbose = FALSE)

# rayshader matrices are rotated/flipped relative to terra rasters. Rather than
# assuming how, test all 8 possible orientations against the DSM itself and keep
# the one that reproduces the raster exactly.
orientations <- list(
  function(m) m,
  function(m) m[nrow(m):1, ],
  function(m) m[, ncol(m):1],
  function(m) m[nrow(m):1, ncol(m):1],
  function(m) t(m),
  function(m) t(m)[ncol(m):1, ],
  function(m) t(m)[, nrow(m):1],
  function(m) t(m)[ncol(m):1, nrow(m):1]
)

dsm_s_ref <- as.matrix(dsm_s, wide = TRUE)
to_grid <- NULL
for (f in orientations) {
  cand <- f(dsm_s_mat)
  if (all(dim(cand) == dim(dsm_s_ref)) &&
      isTRUE(all.equal(unname(cand), unname(dsm_s_ref)))) { to_grid <- f; break }
}
if (is.null(to_grid)) stop("Could not match the rayshader matrix orientation to the raster")

# Shadow raster for one sun position: 1 = lit, 0 = in shadow
# alt_deg = sun elevation above horizon, az_compass = sun azimuth (0 = N, 90 = E, 180 = S)
shade_raster <- function(alt_deg, az_compass) {
  sm <- ray_shade(dsm_s_mat,
                  sunaltitude = unname(alt_deg), sunangle = unname(az_compass),
                  zscale = res(dsm_s)[1],        # horizontal cell size / vertical unit (m)
                  maxsearch = 400,               # max shadow length searched, in cells
                  lambert = FALSE,               # pure cast shadows, no surface-angle shading
                  multicore = TRUE, progbar = FALSE)
  r <- rast(to_grid(sm))
  ext(r) <- ext(dsm_s)
  crs(r) <- crs(dsm_s)
  r
}

# Sun position for day-of-year `doy` at `hour_cet` (local standard time, UTC+1,
# the same clock convention used in the solrad calls below).
# suncalc azimuth is measured from SOUTH towards west, so convert to compass bearing.
sun_pos <- function(doy, hour_cet) {
  t_utc <- as.POSIXct(as.Date(doy - 1, origin = "2026-01-01"), tz = "UTC") + (hour_cet - 1) * 3600
  p <- getSunlightPosition(date = t_utc, lat = lat, lon = lon, keep = c("altitude", "azimuth"))
  c(alt = p$altitude * 180 / pi,
    az  = (p$azimuth * 180 / pi + 180) %% 360)
}

# --- SANITY CHECKS (look at these once before trusting the full run) ----------
sun_pos(172, 12)    # summer noon: altitude ~60-65, azimuth ~175-180 (almost due south)

# Sun in the EAST -> shadows must fall to the WEST of each building (right-to-left
# in the map). If they fall east instead, the shadow direction is flipped.
test_sh <- shade_raster(alt_deg = 15, az_compass = 90)
plot(test_sh, col = c(facet_lut$col[4], "grey95"),
     main = "Shadow test, sun in the EAST (shadows should fall WEST)")
plot(vect(buildings_tile), add = TRUE, border = "black", lwd = 0.3)

# Same ray_shade call, but returns a plain vector in terra cell order (no raster building)
shade_vec <- function(alt_deg, az_compass, maxsearch = 150) {
  sm <- ray_shade(dsm_s_mat,
                  sunaltitude = unname(alt_deg), sunangle = unname(az_compass),
                  zscale = res(dsm_s)[1], maxsearch = maxsearch,
                  lambert = FALSE, multicore = TRUE, progbar = FALSE)
  as.vector(t(to_grid(sm)))
}

################################################################################
# 9. SOLAR GEOMETRY ORIENTATION CHECK (solrad: aspect is relative to SOUTH)
# South-facing vs north-facing, both 35 deg slope, local solar noon in summer.
# South should clearly outperform north.
DirectRadiation(c(172.5, 172.5), Lat = lat, Lon = lon, SLon = 15, DS = 0,
                Elevation = 500, Slope = c(35, 35), Aspect = c(0, 180))


################################################################################
# 10. CONVERT RASTER ASPECT (terra: 0=N clockwise) TO SOLRAD CONVENTION (0=S)
slope_deg        <- values(slope_smooth)[, 1]
aspect_deg_terra <- values(aspect_smooth)[, 1]

# flat cells have no aspect; any value works when slope ~ 0, so don't drop them
aspect_deg_terra[is.na(aspect_deg_terra) & !is.na(slope_deg) & slope_deg < 1] <- 180

valid     <- which(!is.na(slope_deg) & !is.na(aspect_deg_terra))
slp       <- slope_deg[valid]
asp_terra <- aspect_deg_terra[valid]
asp_solrad <- asp_terra - 180
asp_solrad <- ifelse(asp_solrad > 180, asp_solrad - 360, asp_solrad)

# OPTIONAL: panels on flat roofs sit on tilted racks, so a horizontal surface
# understates them. Evaluate flat pixels as if tilted 30° towards south.
# (Assumption to state in your methods; it ignores row-to-row shading and the
# lower usable area on flat roofs.)
is_flat <- slp < 5
slp[is_flat]        <- 30
asp_solrad[is_flat] <- 0

################################################################################
# 11. ANNUAL IRRADIATION WITH SHADING
rep_days      <- as.integer(format(as.Date(paste0("2026-", 1:12, "-15")), "%j"))
days_in_month <- c(31,28,31,30,31,30,31,31,30,31,30,31)

hour_step     <- 2                          # 2 = every 2nd hour (each sample counts for 2 h)
hours         <- seq(6, 18, by = hour_step)
min_shade_alt <- 5                          # below this sun elevation: direct beam treated as blocked
elev_mean     <- mean(values(dtm), na.rm = TRUE)

# Once: which shadow-grid cell does each valid roof pixel sit in?
cells <- cellFromXY(dsm_s, xyFromCell(roof_dsm, valid))

total_irr         <- numeric(length(valid))
total_irr_noshade <- numeric(length(valid))

t0 <- Sys.time()
for (i in seq_along(rep_days)) {
  d <- rep_days[i]
  w <- days_in_month[i]
  month_total  <- numeric(length(valid))
  month_noshad <- numeric(length(valid))
  
  for (h in hours) {
    doy_frac <- d + h / 24
    alt <- Altitude(doy_frac, Lat = lat, Lon = lon, SLon = 15, DS = 0)
    sp  <- sun_pos(d, h)
    if (alt <= 0 || sp["alt"] <= 0) next
    
    lit <- if (sp["alt"] < min_shade_alt) {
      rep(0, length(valid))
    } else {
      v <- shade_vec(sp["alt"], sp["az"])[cells]
      v[is.na(v)] <- 1
      v
    }
    
    doy_vec <- rep(doy_frac, length(valid))
    Sdir <- DirectRadiation(doy_vec, Lat = lat, Lon = lon, SLon = 15, DS = 0,
                            Elevation = elev_mean, Slope = slp, Aspect = asp_solrad)
    Sdif <- DiffuseRadiation(doy_vec, Lat = lat, Lon = lon, SLon = 15, DS = 0,
                             Elevation = elev_mean, Slope = slp)
    
    month_total  <- month_total  + (Sdir * lit + Sdif) * hour_step
    month_noshad <- month_noshad + (Sdir       + Sdif) * hour_step
  }
  
  total_irr         <- total_irr         + month_total  * w
  total_irr_noshade <- total_irr_noshade + month_noshad * w
  cat("Month", i, "of 12 done --", round(difftime(Sys.time(), t0, units = "mins"), 1), "min elapsed\n")
}

################################################################################
# 12. PER-BUILDING AGGREGATION (whole-roof summary stats)
px_area <- res(irr_raster)[1] * res(irr_raster)[2] # m^2 per pixel

buildings_tile$mean_irr_Wm2 <- extract(irr_raster, vect(buildings_tile), fun = mean, na.rm = TRUE)[, 2]
buildings_tile$roof_area_m2 <- extract(irr_raster, vect(buildings_tile),
                                       fun = function(x) sum(!is.na(x)) * px_area)[, 2]

buildings_tile$irr_kWh_m2_yr <- buildings_tile$mean_irr_Wm2 / 1000

performance_ratio <- 0.80 # inverter/soiling/temperature losses
panel_efficiency  <- 0.20 # modern mono-Si panels

buildings_tile$pv_potential_kWh_yr <- buildings_tile$irr_kWh_m2_yr *
  buildings_tile$roof_area_m2 * panel_efficiency * performance_ratio

summary(buildings_tile$irr_kWh_m2_yr) # re-check the breaks below after shading!

suit_breaks <- c(-Inf, 600, 900, 1200, Inf)
suit_labels <- c("Low", "Moderate", "Good", "Excellent")

buildings_tile$suitability <- cut(buildings_tile$irr_kWh_m2_yr,
                                  breaks = suit_breaks, labels = suit_labels)

ggplot(buildings_tile) +
  geom_sf(aes(fill = suitability), color = "grey40", linewidth = 0.1) +
  scale_fill_manual(values = suit_cols, drop = FALSE) +
  theme_minimal() +
  labs(title = "Rooftop Solar Suitability (per building)", fill = "Suitability")

mapview(buildings_tile, zcol = "suitability", col.regions = suit_cols)


################################################################################
# 13. PER-ROOF-FACET SEGMENTATION
#
# NOTE on method: terra::patches()/raster::clump() only split regions by NA gaps,
# NOT by value, so they would merge a whole gap-free roof into one blob. Instead:
# dissolve by value with as.polygons(dissolve = TRUE), then explode multi-part
# polygons into single-part features with st_cast("POLYGON").
min_facet_area <- 8   # m^2, drops slivers too small for a panel anyway

# Aspect: 4 compass sectors, UNsmoothed (keeps the ridge-line boundary sharp)
aspect_roof <- mask(aspect, vect(buildings_tile))
aspect_bin <- classify(aspect_roof, matrix(c(
  0,45,1,   45,135,2,   135,225,3,   225,315,4,   315,360,1   # N / E / S / W
), ncol = 3, byrow = TRUE))

# Slope: flat vs. pitched (smoothed is fine for a coarse distinction)
slope_roof <- mask(slope_smooth, vect(buildings_tile))
slope_bin <- classify(slope_roof, matrix(c(
  0, 15, 1,      # flat / near-flat
  15, 90, 2      # pitched
), ncol = 3, byrow = TRUE))

# Aspect is meaningless on near-flat pixels (noise flips it between bins),
# so give them their own aspect class 0
aspect_bin[slope_roof < 5] <- 0

facet_class <- aspect_bin * 10 + slope_bin

# Majority filter: removes salt-and-pepper noise, keeps real facet boundaries
facet_class_clean <- focal(facet_class, w = 3, fun = "modal", na.rm = TRUE)

##### visualise the facet classes with the palette #####
f <- facet_class_clean
present <- facet_lut[facet_lut$id %in% unique(values(f)), ]
levels(f) <- present[, c("id", "label")]
coltab(f) <- data.frame(value = present$id, col = present$col)

plot(f, main = "Facet classes (aspect × slope)")
plot(vect(buildings_tile), add = TRUE, border = "black", lwd = 0.3)

##### raster -> facet polygons #####
# Bake the building ID into the value so facets never merge across buildings
building_ids <- rasterize(vect(buildings_tile), facet_class_clean, field = "build_id")
composite <- building_ids * 1000 + facet_class_clean

facet_polys <- as.polygons(composite, dissolve = TRUE) |> st_as_sf()
facet_polys <- suppressWarnings(st_cast(facet_polys, "POLYGON"))
facet_polys <- st_set_crs(facet_polys, 25832)
facet_polys$facet_id <- seq_len(nrow(facet_polys))

facet_polys$area_m2 <- as.numeric(st_area(facet_polys))
facet_polys <- facet_polys[facet_polys$area_m2 > min_facet_area, ]


################################################################################
# 14. PER-FACET IRRADIATION AND SUITABILITY
facet_polys$mean_irr_Wm2 <- terra::extract(irr_raster, vect(facet_polys),
                                           fun = mean, na.rm = TRUE)[, 2]
# Alternative, more precise on partial-pixel edges (needs exactextractr):
# library(exactextractr)
# facet_polys$mean_irr_Wm2 <- exact_extract(irr_raster, facet_polys, "mean")

facet_polys$irr_kWh_m2_yr <- facet_polys$mean_irr_Wm2 / 1000
facet_polys$pv_potential_kWh_yr <- facet_polys$irr_kWh_m2_yr *
  facet_polys$area_m2 * panel_efficiency * performance_ratio

summary(facet_polys$irr_kWh_m2_yr)   # check value range before trusting the breaks

facet_polys$suitability <- cut(facet_polys$irr_kWh_m2_yr,
                               breaks = suit_breaks, labels = suit_labels)

# per-facet suitability with building outlines
ggplot() +
  geom_sf(data = buildings_tile, fill = NA, color = "grey40", linewidth = 0.2) +
  geom_sf(data = facet_polys, aes(fill = suitability), color = NA) +
  scale_fill_manual(values = suit_cols, drop = FALSE) +
  theme_minimal() +
  labs(title = "Per-facet Solar Suitability", fill = "Suitability")


################################################################################
# 15. SAVE RESULTS
st_write(buildings_tile, paste0(data_dir, "solar_potential_per_building_", tile_id, ".shp"), delete_layer = TRUE)
st_write(facet_polys,    paste0(data_dir, "solar_potential_per_facet_",    tile_id, ".shp"), delete_layer = TRUE)

save(aoi, dtm, dsm, roof_dsm, buildings_tile, facet_polys,
     slope, aspect, slope_smooth, aspect_smooth,
     irr_raster, irr_raster_noshade, lat, lon,
     file = paste0(data_dir, "workspace_", tile_id, ".RData"))


################################################################################
# 16. 3D VIEWS (rayshader)
dsm_matrix <- raster_to_matrix(dsm)

##### Option A: continuous irradiation draped over the DSM #####
# irr_raster only covers the roof extent, so resample it to the full DSM grid
irr_full   <- resample(irr_raster, dsm)
irr_matrix <- raster_to_matrix(irr_full)

irr_colors <- height_shade(irr_matrix, texture = irr_ramp(256))

dsm_matrix |>
  sphere_shade(texture = "imhof1") |>
  add_overlay(irr_colors, alphalayer = 0.85) |>
  plot_3d(dsm_matrix, zscale = 0.5, fov = 0, theta = 45, phi = 35,
          windowsize = c(1000, 800), zoom = 0.65, background = "white")

render_snapshot(paste0(data_dir, "solar_3d_irradiation.png"))
rgl::close3d()

##### Option B: per-facet suitability classes draped over the DSM #####
suitability_num <- facet_polys
suitability_num$suit_code <- as.numeric(suitability_num$suitability)   # 1 = Low ... 4 = Excellent

suit_raster <- rasterize(vect(suitability_num), dsm, field = "suit_code")
suit_matrix <- raster_to_matrix(suit_raster)

# Map each class code to its colour; NA (non-roof) stays transparent
suit_overlay <- matrix(unname(suit_cols)[suit_matrix], nrow = nrow(suit_matrix))
suit_overlay[is.na(suit_matrix)] <- NA

dsm_matrix |>
  sphere_shade(texture = "imhof1") |>
  add_overlay(suit_overlay, alphalayer = 0.9) |>
  plot_3d(dsm_matrix, zscale = 0.5, fov = 0, theta = 45, phi = 35,
          windowsize = c(1000, 800), zoom = 0.65, background = "white")

render_snapshot(paste0(data_dir, "solar_3d_suitability.png"))

# option c
nr <- nrow(suit_matrix); nc <- ncol(suit_matrix)
ok <- !is.na(suit_matrix)

rgb_vals <- t(col2rgb(unname(suit_cols)[suit_matrix[ok]])) / 255

r <- g <- b <- a <- matrix(0, nr, nc)
r[ok] <- rgb_vals[, 1]; g[ok] <- rgb_vals[, 2]; b[ok] <- rgb_vals[, 3]
a[ok] <- 1                                   # non-roof stays fully transparent

suit_overlay <- array(c(r, g, b, a), dim = c(nr, nc, 4))

dsm_matrix |>
  sphere_shade(texture = "imhof1") |>
  add_overlay(suit_overlay, alphalayer = 0.9) |>
  plot_3d(dsm_matrix, zscale = 0.5, fov = 0, theta = 45, phi = 35,
          windowsize = c(1000, 800), zoom = 0.65, background = "white")

render_snapshot(paste0(data_dir, "solar_3d_suitability.png"))
