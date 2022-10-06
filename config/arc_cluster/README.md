# Running open-gira on ARC

As open-gira is built with snakemake, its use is remarkably similar from a
laptop to a cluster. However there are a few differences in setup and use. They
are discussed here.

## Python environment

### Initialising our shared conda installation

There is a conda install which users may share. This means we don't need to
create many duplicate environments unnecessarily (snakemake will check to see
if an equivalent environment has already been created).

Prior to using the conda install for the first time you must initialise it for
your shell with the following command:
`/data/ouce-gri-jba/anaconda/condabin/conda init`

Your `~/.bashrc` should then contain someting like this:
```
i# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
    __conda_setup="$('/data/ouce-gri-jba/anaconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/ouce-gri-jba/anaconda/etc/profile.d/conda.sh" ]; then
        . "/data/ouce-gri-jba/anaconda/etc/profile.d/conda.sh"
    else
        export PATH="/data/ouce-gri-jba/anaconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```

### Enabling an environment with snakemake

We need the snakemake executable to create jobs for us. We could get this from
the module load system, but their version is quite old (6.10.0). N.B. Versions
<7.0.0 may cause the following problem:
https://github.com/snakemake/snakemake/issues/1392

Instead use a conda environment we have created which contains snakemake:
`conda activate snakemake-7.12.1`

Your prompt should then change to something like:
`(snakemake-7.12.1) [cenv0899@arc-login01 ~]$`

## Osmium

open-gira jobs which filter Open Street Map datasets may require the use of a
tool called osmium. This has been compiled on the cluster (with
`/data/ouce-gri-jba/osmium/build_osmium.sh`). To run osmium, place a symlink
somewhere on your $PATH, pointing to the wrapper script. For example:
```
mkdir -p ~/bin
ln -s /data/ouce-gri-jba/osmium/run_osmium.sh ~/bin/osmium
```

## Session persistence

To persist a terminal over time (and despite dropped SSH connections) consider using `tmux`.

Here's a [friendly guide](https://www.hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/) and the [official documentation](https://github.com/tmux/tmux/wiki/Getting-Started).

`Ctrl+b` is used to preface commands (often written as `C-b` in the docs).
`Ctrl+b, d` to detach (keep running, but leave) a session.
`tmux ls` to list running sessions.
`tmux attach-session -t <session_id>` to reattach to a session.

## Allocate resources

Allocate some nodes for use:
`salloc --ntasks-per-node=<max tasks per node> --nodes=<num nodes> --partition=<short|medium|long> --time=01:00:00 --mem=8000`

## Invoke pipeline

Having allocated resources with `salloc` (see above), you can then invoke
snakemake to dispatch jobs and satisfy your target rule. From the open-gira
repository call the command you wish to run, using the cluster specific
profile. For more details on the cluster execution, see the config.yaml file
in the profile directory.
`snakemake --profile config/arc_cluster <target name>`

To test the pipeline with a short job, try the following:
`snakemake --profile config/arc_cluster results/exposure/tanzania-mini_filter-road/hazard-aqueduct-river/img`

## Interpreting errors

Each submitted job will have its `stdout` logged to file. This is very useful
in diagnosing errors.

For example:

```
Error in rule create_overall_bbox:
jobid: 10
output: results/json/tanzania-mini.json
conda-env: /data/ouce-gri-jba/anaconda/envs/1b06bc7edf94bde3731973d2898f9d82
cluster_jobid: 2682538
Error executing rule create_overall_bbox on cluster (jobid: 10, external: 2682538, jobscript: /data/ouce-gri-jba/mert2014/open-gira/.snakemake/tmp.pvgmp5pz/snakejob.create_overall_bbox.10.sh). For error details see the cluster log and the log files of the involved rule(s).
Exiting because a job execution failed. Look above for error message
Complete log: .snakemake/log/2022-10-05T181802.083378.snakemake.log
```

This reports an error in the rule create_overall_bbox so we should look in
logs/create_overall_bbox:

```
ls logs
ls logs/create_overall_bbox
less logs/create_overall_bbox/create_overall_bbox-DATASET\=tanzania-mini\,OUTPUT_DIR\=results-2682538.out `
```

Which prints the full log for that job, where we can look for a traceback:

```
Traceback (most recent call last):
  File "/data/ouce-gri-jba/mert2014/open-gira/.snakemake/scripts/tmp5ruqzjkn.create_overall_bbox.py", line 24, in <module>
    bboxes = subprocess.check_output(["osmium", "fileinfo", osm_file, "-g", "header.boxes"])
  File "/data/ouce-gri-jba/anaconda/envs/1b06bc7edf94bde3731973d2898f9d82/lib/python3.9/subprocess.py", line 424, in check_output
    return run(*popenargs, stdout=PIPE, timeout=timeout, check=True,
  File "/data/ouce-gri-jba/anaconda/envs/1b06bc7edf94bde3731973d2898f9d82/lib/python3.9/subprocess.py", line 505, in run
    with Popen(*popenargs, **kwargs) as process:
  File "/data/ouce-gri-jba/anaconda/envs/1b06bc7edf94bde3731973d2898f9d82/lib/python3.9/subprocess.py", line 951, in __init__
    self._execute_child(args, executable, preexec_fn, close_fds,
  File "/data/ouce-gri-jba/anaconda/envs/1b06bc7edf94bde3731973d2898f9d82/lib/python3.9/subprocess.py", line 1821, in _execute_child
    raise child_exception_type(errno_num, err_msg, err_filename)
FileNotFoundError: [Errno 2] No such file or directory: 'osmium'
```

In particular, here, the FileNotFoundError says that the job runner couldn't find osmium, which is needed to run this rule.
