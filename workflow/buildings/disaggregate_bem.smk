"""
Upscale (disaggregate) building exposure layers according to built volume
"""

rule disaggregate_bem:
    input:
        bem_res="{OUTPUT_DIR}/input/giri/{ISO3}/bem_5x5_valfis_res.tif",
        ghsl_res="{OUTPUT_DIR}/input/ghsl/{ISO3}/GHS_BUILT_V_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        bem_nres="{OUTPUT_DIR}/input/giri/{ISO3}/bem_5x5_valfis_nres.tif",
        ghsl_nres="{OUTPUT_DIR}/input/ghsl/{ISO3}/GHS_BUILT_V_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
    output:
        bem_res_3ss="{OUTPUT_DIR}/buildings/{ISO3}/building_exposure_res_3ss.tif",
        bem_nres_3ss="{OUTPUT_DIR}/buildings/{ISO3}/building_exposure_nres_3ss.tif",
    shell:
        """

        """
