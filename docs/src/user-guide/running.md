# Running

When we run snakemake we have to tell it how many CPU cores it can use to do its processing.
If you have fewer than 4 cores, or wish to use more, substitute an appropriate number for the 4 below:

```shell
snakemake --cores 4 -- <target_file>
```

You should see a lot of text flashing by, and some loading bars.
Eventually, everything should finish with a report along the lines of:

```text
Finished job 0.
111 of 111 steps (100%) done
Complete log: /mnt/f/OxRSE/open-gira/.snakemake/log/2022-01-24T154611.005270.snakemake.log
```