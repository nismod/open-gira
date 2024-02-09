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
        output_dir=$(dirname {output.res})

        wget -nc https://hazards-data.unepgrid.ch/bem_5x5_valfis_res.tif \
            --directory-prefix=$output_dir
        wget -nc https://hazards-data.unepgrid.ch/bem_5x5_valfis_nres.tif \
            --directory-prefix=$output_dir
        """

rule summarise_giri_bem_admin:
    output:
        adm1="{OUTPUT_DIR}/input/giri/bem_5x5_valfis_adm1.csv",
        adm0="{OUTPUT_DIR}/input/giri/bem_5x5_valfis_adm0.csv",
    shell:
        """
        exactextract \
            -p ./results/input/admin-boundaries/adm1.shp \
            -r "res:results/input/giri/bem_5x5_valfis_res.tif" \
            -r "nres:results/input/giri/bem_5x5_valfis_nres.tif" \
            -f GID_1 \
            -o adm1_bem-res.csv \
            -s "sum(res)" \
            -s "sum(nres)"

        exactextract \
            -p ./results/input/admin-boundaries/adm0.shp \
            -r "res:results/input/giri/bem_5x5_valfis_res.tif" \
            -r "nres:results/input/giri/bem_5x5_valfis_nres.tif" \
            -f GID_0 \
            -o adm0_bem-res.csv \
            -s "sum(res)" \
            -s "sum(nres)"
        """
