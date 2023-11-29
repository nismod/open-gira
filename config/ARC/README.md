# Running open-gira on ARC

As open-gira is built using `snakemake`, its use is fairly similar from a
laptop to a cluster. However there are a few differences, notably using a
`profile` (discussed here).

## Python environment

### Micromamba

I recommend installing
[micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html#install-script)
into your userspace as a package manager for Python packages (and more).

### Creating an execution environment

To create an environment on ARC containing the necessary software to run the workflows:
```
micromamba create -f environment.yml -y
```

To activate this:
```
micromamba activate open-gira
```

Your prompt should then change to something like:
```
(open-gira) [cenv0899@arc-login01 ~]$
```

## Exactextract

There is one dependency of `open-gira` that is not available via the conda
(micromamba) ecosystem, `exactextract`. To install this, see
[here](https://github.com/isciences/exactextract#compiling) and place the
compiled binary in your `PATH`.

To build (and run) `exactextract` on ARC you will need to use the `module`
program to load two dependencies:
```
module load GEOS/3.10.3-GCC-11.3.0
module load GDAL/3.5.0-foss-2022a
```

I suggest placing these lines in your `~/.bashrc` file so they automatically run on login.

## Session persistence

To persist a terminal over time (and despite dropped SSH connections) consider using `tmux`.

Here's a [friendly guide](https://www.hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/) and the [official documentation](https://github.com/tmux/tmux/wiki/Getting-Started).

`tmux` to start a new session.

`Ctrl+b` is used to preface commands (often written as `C-b` in the docs).

`Ctrl+b, d` to detach (keep running, but leave) a session.

`tmux ls` to list running sessions.

`tmux attach-session -t <session_id>` to reattach to a session.

## Invoke workflow

The general pattern to doing work with `open-gira` on `ARC` is to activate the
environment (see above) and issue a request for a target file:
```
snakemake --profile config/ARC <target name>
```

`snakemake` will then identify what work is required, issue job requests to
`SLURM` and monitor the filesystem to watch for completed results.

To test the pipeline with a short job, try the following:
```
snakemake --profile config/ARC results/exposure/tanzania-mini_filter-road/hazard-aqueduct-river/img
```

Resource allocation is defined per rule, with defaults in `config/ARC/config.yaml`.

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

This reports an error in the rule `create_overall_bbox` so we should look in
`logs/create_overall_bbox`:

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

In particular, here, the `FileNotFoundError` says that the job runner couldn't
find `osmium`, which is needed to run this rule.
