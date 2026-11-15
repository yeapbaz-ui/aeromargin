#!/usr/bin/env bash
# Cut a small window out of public Mars rasters hosted by USGS Astrogeology, without downloading the full files.
# GDAL reads only the byte ranges it needs over HTTP (/vsicurl/).
#
# Usage: tools/fetch_terrain.sh <name> <center_lon_east_deg> <center_lat_deg> <width_m> <height_m> <dem_res_m> <ti_res_m>
# Output: data/terrain/<name>/{dem,ti}.f32 (little-endian Float32, row-major, north row first),
#         plus GDAL metadata of the windows and of the source files.
set -euo pipefail
NAME=$1; LON=$2; LAT=$3; WM=$4; HM=$5; RES=$6; TIRES=${7:-$6}
OUT=data/terrain/$NAME; mkdir -p "$OUT"
R=3396190
DEM_URL=https://planetarymaps.usgs.gov/mosaic/Mars/HRSC_MOLA_Blend/Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2.tif
TI_BASE=${TI_BASE:-http://astropedia.astrogeology.usgs.gov/download}/Mars/Odyssey/THEMIS-Global-Thermal-Inertia-Mosaic/Quantitative-32-Bit
export GDAL_HTTP_MAX_RETRY=4 GDAL_HTTP_RETRY_DELAY=3 GDAL_DISABLE_READDIR_ON_OPEN=EMPTY_DIR CPL_VSIL_CURL_CACHE_SIZE=200000000

SRS="+proj=eqc +lat_ts=$LAT +lat_0=0 +lon_0=$LON +R=$R +units=m +no_defs"
read -r XMIN YMIN XMAX YMAX < <(python3 -c "
import math
lat=float('$LAT'); wm=float('$WM'); hm=float('$HM'); R=$R
yc=R*math.radians(lat)
print(-wm/2, yc-hm/2, wm/2, yc+hm/2)")
echo "target window (local equirectangular, metres): $XMIN $YMIN $XMAX $YMAX"

if [ "${SKIP_DEM:-0}" != "1" ]; then
echo "== DEM"
gdalinfo -json "/vsicurl/$DEM_URL" > "$OUT/source_dem.gdalinfo.json"
gdalwarp -q -overwrite -t_srs "$SRS" -te "$XMIN" "$YMIN" "$XMAX" "$YMAX" -tr "$RES" "$RES" -r bilinear -ot Float32 -dstnodata -9999 "/vsicurl/$DEM_URL" "$OUT/dem.tif"
gdalinfo -json -stats "$OUT/dem.tif" > "$OUT/dem.gdalinfo.json"
python3 tools/tif_to_f32.py "$OUT/dem.tif" "$OUT/dem.f32"
echo "$DEM_URL" > "$OUT/dem.source.txt"

fi

if [ "${SKIP_TI:-0}" != "1" ]; then
echo "== THEMIS thermal inertia"
echo "-- directory probes"
for D in https://planetarymaps.usgs.gov/mosaic/Mars/ https://planetarymaps.usgs.gov/mosaic/Mars/THEMIS_controlled_mosaics/ https://planetarymaps.usgs.gov/mosaic/; do echo "## $D"; curl -sL -m 40 "$D" | grep -o -i 'href="[^"]*' | grep -i -E "themis|therm|TI_|inertia|TES" | head -40 || true; done
: > "$OUT/ti.source.txt"
for Q in ${TI_QUADS:-30S240E}; do
  URL="$TI_BASE/THEMIS_TI_Mosaic_Quant_${Q}_100mpp.cub"
  echo "-- $URL"
  curl -sIL -m 60 "$URL" | grep -i -E "^HTTP|content-length|accept-ranges|location|content-type" || true
  curl -sL -m 60 -r 0-1023 -o /tmp/head.bin -w "range probe: http %{http_code}, %{size_download} bytes\n" "$URL" || true
  head -c 600 /tmp/head.bin | tr -c '[:print:]\n' '.' | head -20 || true
  if timeout 240 gdalinfo -json "/vsicurl/$URL" > "$OUT/source_ti_$Q.gdalinfo.json" 2> "$OUT/source_ti_$Q.err"; then
    python3 -c "
import json;d=json.load(open('$OUT/source_ti_$Q.gdalinfo.json'));print('size',d['size'],'gt',d.get('geoTransform'),'corners',d.get('wgs84Extent',d.get('cornerCoordinates')))"
    if timeout ${TI_TIMEOUT:-600} gdalwarp -overwrite -t_srs "$SRS" -te "$XMIN" "$YMIN" "$XMAX" "$YMAX" -tr "$TIRES" "$TIRES" -r bilinear -ot Float32 -dstnodata -9999 "/vsicurl/$URL" "$OUT/ti.tif"; then
      gdalinfo -json -stats "$OUT/ti.tif" > "$OUT/ti.gdalinfo.json" || true
      VALID=$(python3 -c "
import json;d=json.load(open('$OUT/ti.gdalinfo.json'));b=d['bands'][0]
print(b.get('metadata',{}).get('',{}).get('STATISTICS_VALID_PERCENT','0'))")
      echo "valid percent: $VALID"
      if python3 -c "import sys; sys.exit(0 if float('$VALID')>${TI_MIN_VALID:-20} else 1)"; then
        python3 tools/tif_to_f32.py "$OUT/ti.tif" "$OUT/ti.f32"; echo "$URL" > "$OUT/ti.source.txt"; break
      fi
    else
      echo "warp failed or timed out for $URL"
    fi
  else
    echo "could not open $URL"; cat "$OUT/source_ti_$Q.err" | head -5 || true; rm -f "$OUT/source_ti_$Q.gdalinfo.json"
  fi
done
fi
rm -f "$OUT"/*.tif "$OUT"/*.aux.xml "$OUT"/*.err
ls -la "$OUT"
