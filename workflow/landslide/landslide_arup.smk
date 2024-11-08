
rule download_metadata:
    output:
        json="{OUTPUT_DIR}/input/hazard-landslide-arup/raw/metadata.json",
    shell:
        """
        wget --output-document={output.json} \
            https://datacatalogapi.worldbank.org/ddhxext/DatasetDownload?dataset_unique_id=0037584&version_id=
        """

rule download_landslides:
    input:
        json="{OUTPUT_DIR}/input/hazard-landslide-arup/raw/metadata.json",
    output:
        "{OUTPUT_DIR}/input/hazard-landslide-arup/raw/global-landslide-hazard-map-report.pdf",
        "{OUTPUT_DIR}/input/hazard-landslide-arup/raw/ls_eq_tiled.tif",
        "{OUTPUT_DIR}/input/hazard-landslide-arup/raw/LS_RF_Median_1980-2018_COG.tif",
        # Also available:
        # "{OUTPUT_DIR}/input/hazard-landslide-arup/raw/LS_RF_Mean_1980-2018_COG.tif",
        # "{OUTPUT_DIR}/input/hazard-landslide-arup/raw/LS_TH_COG.tif",
    run:
        import json
        import os
        import zipfile
        from pathlib import Path
        import requests

        out_dir = os.path.dirname(input.json)
        output_fnames = set(os.path.basename(fname) for fname in output)

        with open(input.json, 'r') as fh:
            meta = json.load(fh)

        for file_meta in meta["resources"]:
            fname = file_meta["distribution"]["file_name"]
            url = file_meta["distribution"]["url"]
            out_file = os.path.join(out_dir, fname)

            if Path(out_file).exists() or fname not in output_fnames:
                print("Skipped downloading", fname)
            else:
                print("Downloading", url)
                r = requests.get(url)
                with open(out_file, 'wb') as fd:
                    for chunk in r.iter_content(chunk_size=1024):
                        fd.write(chunk)

            if ".zip" in fname:
                print("Extracting zip", out_file)
                with zipfile.ZipFile(out_file, 'r') as zh:
                    zh.extractall(out_dir)
