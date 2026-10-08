# Data

No data is stored in this repository (the point clouds are large). Download it yourself and place it in this folder:

| File | Dataset | Where to get it |
|---|---|---|
| `651_5478.laz` | Laserdaten (ALS point cloud, 1 km × 1 km tile) | [Bayern OpenData – Laserdaten](https://geodaten.bayern.de/opengeodata/OpenDataDetail.html?pn=laserdaten) |
| `hausumringe.shp` (+ `.dbf`, `.shx`, `.prj`) | Hausumringe (building outlines) | [Bayern OpenData](https://geodaten.bayern.de/opengeodata/) |

Notes:

- The portal offers bulk downloads as a **Metalink** file, a small XML list of download URLs. Open it with a download manager, or read the URLs in R and download only the tiles you need.
- The `.laz` header carries no CRS. Bavarian data is in **ETRS89 / UTM zone 32N (EPSG:25832)**; the script assigns this manually.
- Data accessed on 14.09.2026.
- To use another tile, change `tile_id` at the top of `R/rooftop_solar_potential.R`.
