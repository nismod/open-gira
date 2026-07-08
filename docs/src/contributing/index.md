# Contributing to open-gira

Contributions are welcome — bug reports, documentation, new hazards, new asset
types, new sectors. This page describes how to set up a development environment,
how the tests are structured, and the conventions a contribution is expected to
follow. It reflects what the continuous-integration checks actually enforce, so
that "it passes locally" and "it passes in CI" mean the same thing.

## Development setup

open-gira uses [pixi](https://pixi.sh/) to manage its environment from a
lockfile. From a clone of the repository:

```bash
pixi install        # create the environment from the lockfile
pixi shell          # enter a shell with it activated
```

Work inside `pixi shell` (or prefix commands with `pixi run`). The full list of
dependencies is in `pixi.toml`, and the locked versions in `pixi.lock`.

## Before you open a pull request

The CI pipeline (`.github/workflows/test.yml`) runs three checks. Run their
equivalents locally first:

1. **Formatting.** CI runs `ruff format --check`. Format your code before
   committing:

   ```bash
   ruff format .
   ```

2. **Environment consistency.** CI regenerates `environment.yml` from the pixi
   workspace and fails if it differs from what is committed. If you change
   dependencies, keep the two in sync:

   ```bash
   pixi workspace export conda-environment --name open-gira > environment.yml
   ```

3. **Tests.** CI runs the full test suite in parallel:

   ```bash
   python -m pytest -n auto
   ```

   Some tests need `osmium` and `imagemagick` available; see
   [Installation](../user-guide/installation.md) for these.

If all three pass locally, CI should agree.

## How the tests are structured

This is the part of the codebase most often described as "tribal knowledge," so
it is worth spelling out. `tests/README.md` is the fuller account; the essentials
are:

- **One test per rule.** For each rule in `workflow/Snakefile` there is a test
  that runs the rule on a small sample dataset and checks its outputs against a
  stored expectation.
- **The sample dataset is Djibouti** — small, but with road and rail networks, a
  coastline, and land borders, so it exercises real logic quickly.
- **The fixture convention.** Each rule's test lives in
  `tests/integration/<rule_name>/`, alongside a Python file
  `tests/integration/test_<rule_name>.py`. The directory contains two folders:
  - `data/` — the input files the rule consumes (or an empty `.gitkeep` if the
    rule takes no inputs);
  - `expected/` — the output files the rule should produce.

  The test runner copies the fixture into a temporary working directory, runs
  the rule via snakemake, and compares the produced outputs against `expected/`.

### Adding a test for a new rule

When you add a rule, add a test:

1. Extend `tests/config/config.yaml` if your rule needs new data sources listed.
2. Run the rule once (with the test configuration) to generate correct outputs
   to use as the expectation.
3. Create `tests/integration/<rule_name>/` with `data/` and `expected/` folders.
4. Create `tests/integration/test_<rule_name>.py`, copying an existing test and
   changing the function name and the `run_snakemake_test` arguments to your
   rule name and target paths.

`tests/README.md` gives the step-by-step, including the semi-automatic route via
snakemake's `--generate-unit-tests`.

## Extending the workflow

For common extensions there are focused how-to guides:

- [Add an asset type and damage curve](add-asset-type.md)

*(further extension guides — adding a hazard dataset, adding an OSM filter — are
planned; see the documentation outline.)*

## Reporting issues

Bug reports and feature requests are welcome on the
[issue tracker](https://github.com/nismod/open-gira/issues). A good bug report
includes the target you requested, the command you ran, and what happened versus
what you expected. A dry run (`snakemake -n -- <target>`) attached to the report
is often enough to diagnose a problem.
