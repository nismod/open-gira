# Cleaning the output

The output directory `./results/` is orderly but large.

You can remove intermediate files using the command:
```shell
snakemake -c1 -R clean
```

To remove everything, simply delete the `./results/` directory:
```shell
rm -rf results
```
