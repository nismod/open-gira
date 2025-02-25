# Count number of EAD files in each direct damages run, in planet-latest... directory
for dirn in results/direct_damages/planet-latest_filter-*/*/EAD_*; do sh -c "echo $dirn; ls $dirn | wc -l"; done
