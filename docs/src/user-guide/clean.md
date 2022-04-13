# Cleaning the output

The output directory `./results/` is orderly but large.
You can remove all the directories except `./results/input/` using the command:
```shell
snakemake -R clean
```  
All directories _including_ `./results/input/` can be removed using the command:
```shell
snakemake -R clean_all
``` 

Neither of these commands will remove the output, the final
`./results/<dataset>_filter-highway-core_hazard-<hazard>.geoparquet` files.
To remove everything, simply delete the `./results/` directory:
```shell
rm -rf results
```