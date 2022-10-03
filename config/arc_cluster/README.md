# Running open-gira on ARC

## Persistence

Either create or attach to a tmux session:
`tmux attach-session -t <id>`
`tmux`

## Resources

Allocate some nodes for use:
`salloc --ntasks-per-node=<max tasks per node> --nodes=<num nodes> --partition=<short|medium|long> --time=01:00:00 --mem=8000`

## Python environment

We need the snakemake executable to create jobs for us. We could get this from
the module load system, but their version is quite old (6.10.0). N.B. Versions
<7.0.0 may cause the following problem:
https://github.com/snakemake/snakemake/issues/1392

Instead use a conda environment we have created which contains snakemake:
`source activate /data/ouce-gri-jba/envs/snakemake-7.12.1`

## Osmium

N.B. Some open-gira jobs require the use of osmium. This has been compiled on
the cluster with `$DATA/osmium/build_osmium.sh`. The executable (actually a
wrapper) is in ~/bin. N.B. To _run_ osmium, the Boost C++ libraries must be
available. They are provided by the cluster module system as
Boost/1.77.0-GCC-11.2.0, and loaded as a step in the wrapper script.

## Invocation

From the open-gira repo call the command we wish to run, using the cluster
specific profile. For more details on the cluster execution, see the
config.yaml file in the profile directory.
`snakemake --profile config/arc_cluster <target name>`
