# Testing

For each rule used in `workflow/Snakefile` there is a unit test ensuring that
the expected outputs are generated, given inputs based on the test dataset.

Unit tests are adapted from those automatically generated using `snakemake
--generate-unit-tests`, see [automatically generating unit
tests](https://snakemake.readthedocs.io/en/stable/snakefiles/testing.html) in
the Snakemake docs.

## Test dataset

The pipeline is tested on a small sample of OpenStreetMap data; the country of
Djibouti in the Horn of Africa.

Djibouti has both road and rail networks, land borders with 3 countries and a
coastline, but is quite small. Its size means the datasets are lightweight and
the runtimes brief.

## External files

Some of the tests rely on some external files, i.e. files that should not be a
part of the workflow because they are not produced by the workflow but are
instead located elsewhere. The download_* rules, for instance, import files
from either local or remote locations. For testing purposes, we store those
remote files outside the working directory of the test, in
`tests/external_files`.

The contents of this directory are:

| filename | description |
|----------|-------------|
| `djibouti-latest.osm.pbf` | OpenStreetMap test data (see [above](#test-dataset)) |
| `hazard_sources.txt` | Text file that contains the file path to the hazard file below |
| `inunriver_rcp4p5_MIROC-ESM-CHEM_2030_rp01000.tif` | Hazard data raster file |

## Configuration

The tests have their own configuration, located in `tests/config/config.yaml`.
This is largely a copy of the production config file, although there are some
important differences.

The testing config file identifies the files in `tests/external_files` as
remote targets for download rules and network_filter expression files e.g.
`tests/config/road-secondary.txt`.

The slice count is set to 4. This means that splitting and joining logic is
tested, but the number of slices remains small for simplicity.

## Adding new tests

### Manually

When you add a new rule, you should add a new test for the rule.
That is done as follows:

1. Add to the testing config file as necessary to list your infrastructure
   and hazard data sources. Changing what's already there might break extant
   tests.
2. Run a snakemake job that executes your rule and creates the output files you
need. To this end you may find it useful to run a job from open-gira/ with the
testing configuration. This can be achieved as follows: `snakemake --cores all
--configfile tests/config/config.yaml --directory . <desired_output>`
3. Create a new directory `tests/unit/<rule_name>` where `<rule_name>` is
the name of your new rule (without angle brackets).
4. Create a new Python file `tests/unit/test_<rule_name>.py`.
5. Inside your `tests/unit/<rule_name>` directory, create two folders:
   1. `data`, containing any files your rule takes as inputs (or an empty
   `.gitkeep` file if there are no inputs).
   2. `expected` containing any outputs your rule produces
6. Copy the contents of one of the other test files into your
`tests/unit/test_<rule_name>.py` file.
   1. Change the function name to `test_<rule_name>`.
   2. Change the first argument of the `run_snakemake_test` call to `<rule_name>`.
   3. Change the second argument of the `run_snakemake_test` call to be a tuple
   of the desired output paths i.e. the target files or directories.

### Automatically (sort of)

Running `snakemake` with the `--generate-unit-tests` options automatically
create unit tests for each rule in the workflow file.  Input data and expected
outputs are copied into a specific directory under `.tests/unit`.

N.B. Before generating the tests, the pipeline must have been successfully run
at least once, in order for the expected outputs to be there.

The generated tests can then be adapted from `.tests/unit` into our
`tests/unit`. Some modifications are necessary to create the necessary
structure, see [above](#manually).

## Test runner

The tests rely on a runner and common output checking mechanism, located in the
automatically generated `tests/unit/runner.py`.  This has been adapted to
reduce repetition in the test files, by adding in a `run_snakemake_test`
function that takes the test target name and targets and runs the three-stage
test procedure for each.

The test procedure stages are:
1. Setup temporary test environment and copy required files
2. Run the snakemake command being tested
3. Check output against an expected output

The first stage additionally copies across the `tests/config`,
`tests/external_files` and `tests/bundled_data` directories to the temporary
working directory.  The last stage is specifically told to ignore those
directories when testing the output.  This will not be a problem unless a rule
attempts to modify the configuration files, which it should not do.

## Continuous Integration

Tests are run any time a branch is updated, as defined by
`.github/workflows/test.yml`.
