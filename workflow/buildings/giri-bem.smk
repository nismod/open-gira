"""
Download GIRI Building Exposure Model

Total Building Stock US$ (Resolution: 5km x 5km)

The Building Exposure Model provides information on the building type and the
economic value of the built environment for the non-residential (employment,
health and education) sectors, for each country and territory of the world.

This dataset was produced by UNEP/GRID-Geneva in May 2023.

Source
------
https://giri.unepgrid.ch

Reference
---------
Thomas Piller, Antonio Benvenuti & Andrea De Bono (2023) The GIRI global
building exposure model (BEM)
https://giri.unepgrid.ch/sites/default/files/2023-09/GIRI_BEM_report_UNIGE.pdf
"""

rule download_giri_bem:
    output:
        res="{OUTPUT_DIR}/input/giri/bem_5x5_valfis_res.tif",
        nres="{OUTPUT_DIR}/input/giri/bem_5x5_valfis_nres.tif",
    shell:
        """
        output_dir=$(dirname {output})

        wget -nc https://hazards-data.unepgrid.ch/bem_5x5_valfis_res.tif \
            --directory-prefix=$output_dir
        wget -nc https://hazards-data.unepgrid.ch/bem_5x5_valfis_nres.tif \
            --directory-prefix=$output_dir
        """
