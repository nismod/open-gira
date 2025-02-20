
for filter in rail road-tertiary; do
for hazard in river-jrc aqueduct-river aqueduct-coast coastal-deltares; do
for dir in EAD EAD_and_cost_per_RP cost_per_RP fraction_per_RP; do
    echo "Checking $filter $hazard $dir"
    python scripts/check_missing_slices.py results/direct_damages/planet-latest_filter-$filter/hazard-$hazard/$dir
done
done
done

for filter in rail road-tertiary; do
for hazard in landslide-arup; do
for dir in EAD_and_cost_per_trigger split_EAD_and_cost_per_trigger; do
    echo "Checking $filter $hazard $dir"
    python scripts/check_missing_slices.py results/direct_damages/planet-latest_filter-$filter/hazard-$hazard/$dir
done
done
done
