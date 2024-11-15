"""
Download JRC Global Human Settlement Layer Built-up Surface

GHS-BUILT-S R2023A - GHS built-up surface grid, derived from Sentinel2 composite
and Landsat, multitemporal (1975-2030)

Reference
---------
https://ghsl.jrc.ec.europa.eu/download.php?ds=bu

Dataset:

> Pesaresi M., Politis P. (2023): GHS-BUILT-S R2023A - GHS built-up surface
> grid, derived from Sentinel2 composite and Landsat, multitemporal
> (1975-2030)European Commission, Joint Research Centre (JRC) PID:
> http://data.europa.eu/89h/9f06f36f-4b11-47ec-abb0-4f8b7b1d72ea,
> doi:10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA

Concept & Methodology:

> European Commission GHSL Data Package 2023, Publications Office of the
> European Union, Luxembourg, 2023, JRC133256, ISBN 978-92-68-02341-9
> doi:10.2760/098587

"""

rule download_ghsl_built_s:
    output:
        "{OUTPUT_DIR}/input/ghsl/GHS_{VAR}_E{YEAR}_GLOBE_R2023A_4326_3ss_V1_0.tif"
    params:
        GROUP=lambda wildcards, output: wildcards.VAR.replace("_NRES", "").replace("_ANBH", "")
    shell:
        """
        output_dir=$(dirname {output})

        mkdir -p $output_dir

        wget \
            -nc \
            --directory-prefix=$output_dir \
            https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_{params.GROUP}_GLOBE_R2023A/GHS_{wildcards.VAR}_E{wildcards.YEAR}_GLOBE_R2023A_4326_3ss/V1-0/GHS_{wildcards.VAR}_E{wildcards.YEAR}_GLOBE_R2023A_4326_3ss_V1_0.zip

        unzip -o $output_dir/GHS_{wildcards.VAR}_E{wildcards.YEAR}_GLOBE_R2023A_4326_3ss_V1_0.zip \
            -d $output_dir
        """

rule download_ghsl_built_s_all:
    input:
        expand(
            os.path.join(
                "results",
                "input",
                "ghsl",
                "GHS_BUILT_{VAR}_E{YEAR}_GLOBE_R2023A_4326_3ss_V1_0.tif",
            ),
            VAR=("S", "S_NRES", "V", "V_NRES"),
            YEAR=(2020, )
        ) + [
            os.path.join(
                "results",
                "input",
                "ghsl",
                "GHS_BUILT_H_ANBH_E2018_GLOBE_R2023A_4326_3ss_V1_0.tif",
            )
        ]
