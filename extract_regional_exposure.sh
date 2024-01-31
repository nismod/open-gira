set -e
set -x

year="2030"
rcp="rcp4p5"
gcm="MIROC-ESM-CHEM"
level="1"

rcp=$1
gcm=$2
year=$3
level=$4

echo "inunriver_${rcp}_${gcm}_${year} admin${level}"


exactextract \
  -r "rp00002:hazard-aqueduct-river/exposure_0.5m/inunriver_${rcp}_${gcm}_${year}_rp00002.tif" \
  -r "rp00005:hazard-aqueduct-river/exposure_0.5m/inunriver_${rcp}_${gcm}_${year}_rp00005.tif" \
  -r "rp00010:hazard-aqueduct-river/exposure_0.5m/inunriver_${rcp}_${gcm}_${year}_rp00010.tif" \
  -r "rp00025:hazard-aqueduct-river/exposure_0.5m/inunriver_${rcp}_${gcm}_${year}_rp00025.tif" \
  -r "rp00050:hazard-aqueduct-river/exposure_0.5m/inunriver_${rcp}_${gcm}_${year}_rp00050.tif" \
  -r "rp00100:hazard-aqueduct-river/exposure_0.5m/inunriver_${rcp}_${gcm}_${year}_rp00100.tif" \
  -r "rp00250:hazard-aqueduct-river/exposure_0.5m/inunriver_${rcp}_${gcm}_${year}_rp00250.tif" \
  -r "rp00500:hazard-aqueduct-river/exposure_0.5m/inunriver_${rcp}_${gcm}_${year}_rp00500.tif" \
  -r "rp01000:hazard-aqueduct-river/exposure_0.5m/inunriver_${rcp}_${gcm}_${year}_rp01000.tif" \
  -p "admin-boundaries/gadm36_level${level}.shp" \
  -f "GID_${level}" \
  -s "pop_exposed_inunriver_${rcp}_${gcm}_${year}_rp00002=sum(rp00002)" \
  -s "pop_exposed_inunriver_${rcp}_${gcm}_${year}_rp00005=sum(rp00005)" \
  -s "pop_exposed_inunriver_${rcp}_${gcm}_${year}_rp00010=sum(rp00010)" \
  -s "pop_exposed_inunriver_${rcp}_${gcm}_${year}_rp00025=sum(rp00025)" \
  -s "pop_exposed_inunriver_${rcp}_${gcm}_${year}_rp00050=sum(rp00050)" \
  -s "pop_exposed_inunriver_${rcp}_${gcm}_${year}_rp00100=sum(rp00100)" \
  -s "pop_exposed_inunriver_${rcp}_${gcm}_${year}_rp00250=sum(rp00250)" \
  -s "pop_exposed_inunriver_${rcp}_${gcm}_${year}_rp00500=sum(rp00500)" \
  -s "pop_exposed_inunriver_${rcp}_${gcm}_${year}_rp01000=sum(rp01000)" \
  -o "pop_exposed_inunriver_${rcp}_${gcm}_${year}_admin${level}.csv"
