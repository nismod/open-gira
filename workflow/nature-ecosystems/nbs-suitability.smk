"""Download Harwood et al. (2024) NbS opportunities
"""

rule download_nbs_suitability:
    """
    Fetch NbS cost, benefit and opportunity layers.

    biodiversity_benefit_tif: Biodiversity benefit (this is a relative number
        derived from the Biodiversity Habitat Index as the slope of the species
        area curve to give a relative prioritisation)

    carbon_benefit_tif: Additional Carbon of Mature forest after nominal 50
        years (t/ha)

    planting_cost_tif: Native planting costs (Dollars per hectare of restoration
        2020)

    regeneration_cost_tif: Natural regeneration costs (Dollars per hectare of
        restoration 2020)

    tree_potential_tif: Areas suitable for tree planting, which are not already
        tree-covered, crops, urban, or permanent water, and not within 3km of
        the coastline. Intended to identify afforestation areas for
        catchment-scale river flood reduction.

    mangrove_potential_tif: Potential mangrove restoration based on Delft method
        with additional consideration of tidal range. Values are categorical: 1
        for accreting (expanding) shorelines (ideal); 2 for static to moderate
        retreating shorelines (sub ideal); 3 for fast retreating shorelines
        (might be worth considering if coastal managment is applied to block
        fill an otherwise appropriate area)

    landslide_slope_potential_tif: These are maps of currently non-woodland
        areas which could be restored to woodland (under the two scenarios
        potential vegetation and treeline) with slope 8-36 degrees. Values are
        categorical: 1 for Crops, 2 for Anything Else 3 for Bare Ground.
    """
    output:
        biodiversity_benefit_tif="{OUTPUT_DIR}/input/nbs-suitability/raw/G_BioBenefit_9s.tif",
        carbon_benefit_tif="{OUTPUT_DIR}/input/nbs-suitability/raw/G_CarbonBenefit_9s.tif",
        planting_cost_tif="{OUTPUT_DIR}/input/nbs-suitability/raw/G_PlantingCost_9s.tif",
        regeneration_cost_tif="{OUTPUT_DIR}/input/nbs-suitability/raw/G_RegenCost_9s.tif",
        tree_potential_tif="{OUTPUT_DIR}/input/nbs-suitability/raw/G_PotentialNonCoastalTreeNBS_9s.tif",
        mangrove_potential_tif="{OUTPUT_DIR}/input/nbs-suitability/raw/ManRestorClass_9s.tif",
        landslide_slope_potential_tif="{OUTPUT_DIR}/input/nbs-suitability/raw/G_LandslideNbS_123_9s.tif",
    shell:
        """
        # pushd $(dirname {output.mangrove})
        #     zenodo_get -w links.txt --record=XXXXXX # TODO
        #     wget -nc -i links.txt
        #     md5sum -c md5sums.txt
        # popd
        """
