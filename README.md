# AeroMargin

Safe walking routes on Mars terrain, with an oxygen margin for the way home.
Team Cheese Pizza, NASA Space Apps Challenge 2026, "Interplanetary Survival Guide: Martian Map".

Live site: https://yeapbaz-ui.github.io/aeromargin/

This project was created for the NASA Space Apps Challenge. NASA does not endorse this product.

## What it does

An astronaut on a surface walk cannot wait 4 to 24 minutes for Earth to say "turn back". AeroMargin answers one question on its own: can I walk there and still get home?

- Builds the safest route from A to B over real Mars elevation data and shows the shortest one next to it.
- Colours the ground green, yellow or red, and lets the crew add hazard zones that routes then avoid.
- Draws the point of no return: the edge of the area you can walk to and still reach the airlock before the oxygen reserve.
- Shows what a walk costs in distance, time, steps and oxygen, and shrinks the safe area when the astronaut is tired.
- Takes voice or typed commands in English and Russian, and previews the route as it could look on a helmet display.

## Data sources

Every dataset is listed in [`data/manifest.json`](data/manifest.json) and on the [Data & Attribution page](attribution.html).

| Layer | Source | Product | In this build |
| --- | --- | --- | --- |
| Elevation | USGS Astrogeology (NASA MGS MOLA + ESA Mars Express HRSC) | `Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2` | **Data.** A 55 x 36 sample window, loaded and checked |
| Radiation | NASA JPL / SwRI, MSL RAD (Hassler et al. 2014) | 0.64 mSv/day, Gale crater | One reference figure, scaled by the model |
| Thermal inertia | ASU THEMIS / USGS (Fergason et al. 2006) | not retrieved | **Not loaded.** One reference value is used |
| Oxygen model basis | NASA Human Integration Design Handbook | `NASA/SP-2010-3407/REV1` | Basis of the team's model |
| Images | NASA/JPL-Caltech/UArizona, HiRISE | six image IDs | Pictures of the same kind of ground |

### Why this elevation model

The task offered MOLA (463 m), Mars Trek tiles, or the USGS archive. We chose the USGS HRSC and MOLA blended DEM because:

- it is 200 m per pixel, more than twice as fine as MOLA alone, and a walking map needs every metre of that;
- Mars Trek serves pictures of elevation (colour and shading), not the elevation values a route search needs;
- the file is hosted by USGS with a stable address and a published citation, so the tile can be cut again by anyone.

The region is fixed: the north wall and floor of Melas Chasma, centre 8.83 S, 287.39 E, 10.2 x 6.4 km. The window was found by scanning a 300 x 180 km survey of the same DEM for a place with a steep wall, a flat floor and gentle ground for the airlock.

### How the tile is made and checked

`tools/fetch_terrain.sh` uses GDAL to read only the needed byte ranges of the 11 GB source file over HTTP and writes `data/terrain/melas/dem.f32` (Float32, 200 m, metres above the areoid). The GDAL metadata of the window and of the source file is stored next to it. The same script runs in the "Fetch terrain window" GitHub Action.

In the browser, `AeroData.load` reads the manifest, then checks the tile's size, value range and SHA-256 against it. If anything differs, the map does not start and the page says why. There is no synthetic fallback: the procedural terrain generator of the first version has been removed.

## Model assumptions

Everything computed from the elevation data is the team's model and is marked `MODEL` in the interface.

- **Grid.** The 200 m samples are interpolated (Catmull-Rom) to 256 x 160 cells of 40 m. Slope is taken from that grid.
- **Safety classes.** Red: slope over 30 degrees, slope of 25 to 30 degrees (where seasonal streaks can occur; a rule, not a map of RSL sites), and red zones reported by the crew. Yellow: slope over 18 degrees, reported caution zones and a buffer around red zones.
- **Metabolic power.** `W = 115 + 0.28 m g v + 1.95 m g v max(0, sin(slope)) - 0.45 m g v max(0, -sin(slope)) + 38 (452 / I)^1.65`, with 195 kg suited mass in Mars gravity.
- **Oxygen.** `V = W x 60 / 20200 + 0.035` L/min (20.2 kJ per litre of oxygen, 0.035 L/min leak). In metabolic equivalents, 1 MET = 1.163 W per kg for an 80 kg crew member: standing in the suit is 1.2 MET, a steady walk about 3.4 MET. The form follows the NASA Human Integration Design Handbook; the coefficients are ours. The exact HIDH table has not been checked by the team.
- **Soil.** No per-cell soil data is loaded, so every cell uses 295 TIU (average regolith in our specification, from Curiosity REMS at Gale crater).
- **Routes.** Dijkstra search where a step costs litres of oxygen. The safest route multiplies the cost of yellow cells; red cells are never entered.
- **Point of no return.** A cell is inside when oxygen to get there plus oxygen to get home, with a fatigue factor, still leaves the 45 L reserve.
- **Dose.** 0.64 mSv/day (26.7 microSv/h) times a sky-view factor from the elevation grid.
- **Steps.** Route length divided by a 0.6 m stride.

## How to run

The site is static: no build step, no backend, no keys.

```bash
python3 -m http.server 8080     # or: npm run serve
# open http://localhost:8080
```

Open it through a web server, not straight from disk: browsers block loading the terrain tile from `file://`.

Tests need Node 20 or newer and no packages:

```bash
npm test
```

They validate the manifest against the files on disk and check the model: oxygen flow against metabolic rate, the two units of the RAD dose rate, and that the loader rejects tiles that do not match the manifest. The task asked for Vitest or Jest; we used Node's built-in test runner so that the project keeps zero dependencies.

## How to deploy

Push to `main`. GitHub Pages serves the repository root. `.nojekyll` keeps the `data/` folder as it is. The "Tests" workflow runs on every push and pull request.

To share hazard reports between devices, put the address of a Firebase Realtime Database in `AEROMARGIN_DB_URL` at the end of `index.html`. Without it, reports stay in the browser.

## Attribution

- Elevation: Fergason, R. L., Hare, T. M., & Laura, J. (2018). HRSC and MOLA Blended Digital Elevation Model at 200m v2. Astrogeology PDS Annex, U.S. Geological Survey. MOLA data: NASA, CC0. HRSC data: ESA/DLR/FU Berlin, CC BY-SA 3.0 IGO.
- Radiation: Hassler, D. M., et al. (2014). Science, 343(6169), 1244797. MSL RAD, NASA JPL / SwRI.
- Thermal inertia (referenced, not loaded): Fergason, Christensen & Kieffer (2006), doi:10.1029/2006JE002735. THEMIS, ASU.
- Oxygen model basis: NASA (2014), Human Integration Design Handbook, NASA/SP-2010-3407/REV1.
- Images: HiRISE, NASA/JPL-Caltech/UArizona.

See [`CITATION.cff`](CITATION.cff).

## Known limitations

What is real:

- The elevation of every map cell.

What is modelled by the team:

- Slope, safety classes, routes, distance, time, steps, oxygen, the point of no return, dose rate and the science ranking.

What is simulated:

- The tank, the clock, the walk, the accumulated dose, the fatigue index and the telemetry packet.
- Science targets S1 to S5: the team picked their positions on the real grid; their geology is not mapped.
- Two hazard reports are preloaded as samples.

Other limits:

- The elevation model is 200 m per pixel. Boulders, ripples and small scarps are not in it. CTX (20 m) and HiRISE (1 m) elevation models are the next step.
- No thermal inertia tile is loaded, so loose sand is not found from data. The USGS host of the THEMIS quantitative mosaic did not answer when the tiles were cut.
- The radiation figure is an average for Gale crater near solar maximum, not a measurement in Valles Marineris.
- The oxygen model has not been checked against suit test data. Two identifiers are marked "verify" in the manifest instead of being guessed: the PDS data set of RAD and the HIDH table number.
- Voice input needs Chrome, Edge or Safari. The AR view is a preview in the browser, not a headset build.
- Photographs show the same kind of ground from orbit, not the mapped spot.

## Repository layout

```
index.html            the site
attribution.html      Data & Attribution page
app.js                AeroCore (model), AeroData (tile loader), the interface
data/manifest.json    datasets and tiles, with checksums
data/terrain/melas/   the elevation tile and its GDAL metadata
tools/                fetch_terrain.sh, validate_manifest.mjs
tests/                model and loader tests
img/                  HiRISE images
```
