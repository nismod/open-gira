

for fname in nbs_stack landslide_slope_vegetation mangrove river_basin_afforestation landslide_slope_vegetation_with_EAD mangrove_with_EAD river_basin_afforestation_with_EAD; do
    echo "Checking $fname"
    python scripts/check_missing_slice_dirs.py results/slices/planet-latest_nbs $fname.parquet
done
