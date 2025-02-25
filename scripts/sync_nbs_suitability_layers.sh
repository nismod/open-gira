# Sync from project drive to snakemake inputs directory
# - working copy location on SoGE Linux Cluster; to be replaced with Zenodo download rule
source_dir=/ouce-home/projects/LCNR_DR/GLOBAL/300m_9s/RestorationCostBenefit/
target_dir=./results/input/nbs-suitability/raw/

rsync $source_dir/Final\ Mangrove\ Outputs/ManRestorClass_9s.tif $target_dir
rsync $source_dir/Final\ Global\ Outputs/*.tif $target_dir
