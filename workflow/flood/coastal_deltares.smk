"""Deltares Coastal Floods

https://planetarycomputer.microsoft.com/dataset/deltares-floods#overview

Deltares has produced inundation maps of flood depth using a model that takes
into account water level attenuation and is forced by sea level. At the
coastline, the model is forced by extreme water levels containing surge and
tide from GTSMip6. The water level at the coastline is extended landwards to
all areas that are hydrodynamically connected to the coast following a
bathtub-like approach and calculates the flood depth as the difference
between the water level and the topography. Unlike a simple 'bathtub' model,
this model attenuates the water level over land with a maximum attenuation
factor of 0.5m/km-1. The attenuation factor simulates the dampening of the
flood levels due to the roughness over land.
"""

rule download_deltares:
    """Note - URL may need signing as of updates to Planetary Computer storage/access.
    """
    output:
        nc="results/input/hazard-coastal-deltares/raw/raw/GFM_global_{DEM}DEM{RESOLUTION}_{YEAR}slr_{RP}_masked.nc",
    shell:
        """
        base_url="https://deltaresfloodssa.blob.core.windows.net/floods/v2021.06/global/{wildcards.DEM}/{wildcards.RESOLUTION}"
        wget $base_url/GFM_global_{wildcards.DEM}DEM{wildcards.RESOLUTION}_{wildcards.YEAR}slr_{wildcards.RP}_masked.nc
        """

rule convert_deltares:
    input:
        nc="results/input/hazard-coastal-deltares/raw/GFM_global_{DEM}DEM{RESOLUTION}_{YEAR}slr_{RP}_masked.nc",
    output:
        tif="results/input/hazard-coastal-deltares/raw/GFM_global_{DEM}DEM{RESOLUTION}_{YEAR}slr_{RP}_masked.tif",
    shell:
        """
        gdal_translate -co COMPRESS=LZW -a_srs "EPSG:4326" -q {input.nc} {output.tif}
        """

rule prepare_deltares_all:
    input:
        tiffs=expand(
            "results/input/hazard-coastal-deltares/raw/GFM_global_{DEM}DEM{RESOLUTION}_{YEAR}slr_{RP}_masked.tif",
            DEM=["NASA", "MERIT"],
            RESOLUTION=["90m", "1km"],
            YEAR=[2018, 2050],
            RP=[f"rp{rp:04d}" for rp in (0, 2, 5, 10, 25, 50, 100, 250)],
        )
