# Contributor Guide

## Building the documentation locally

The documentation is built with [Sphinx](https://www.sphinx-doc.org) using the
[MyST parser](https://myst-parser.readthedocs.io) so that Markdown files are
rendered directly.

### Set up the environment

A conda environment with all required packages is provided in `docs/environment.yaml`:

```bash
conda env create -f docs/environment.yaml
conda activate hycom-cice-docs
```

If you already have a conda environment for the model, you can install the docs
dependencies into it instead:

```bash
conda activate <your-env>
pip install -r docs/requirements.txt
```

### Build

Run from the `docs/` directory:

```bash
cd docs
make html
```

Then open `docs/_build/html/index.html` in a browser. Other useful targets:

```bash
make clean   # remove the build directory
make latex   # build a PDF via LaTeX
make help    # list all available targets
```

### Previewing documentation changes in a PR

To see how your changes to the documentation render, you have two options:

1. Build the documentation locally — see [Build](#build) above for instructions.
2. After pushing your changes to the PR, once the Read the Docs build has finished,
   click the yellow link in the PR's checks list (as shown below) to preview the
   rendered docs for this PR:

   <img width="926" height="449" alt="Read the Docs PR check with the Details link highlighted" src="https://github.com/user-attachments/assets/bf73471a-0bee-4dd4-a386-7ce50c019566" />

### Adding or editing pages

- All documentation lives in `docs/` as Markdown files.
- The table of contents is defined in `docs/index.rst`.
- To add a new page, create a `.md` file in `docs/` and add its name (without
  extension) to the appropriate `toctree` block in `docs/index.rst`.

## Development environment

To work on the code locally, set up the Python environment from the repo root:

```bash
conda env create -f environment/python.yaml
conda activate hycom-cice
```

Then install the local Python libraries:

```bash
pip install pythonlibs/modeltools pythonlibs/modelgrid pythonlibs/gridxsec pythonlibs/abfile
```

See [installation.md](installation.md) for HPC-specific setup (Betzy, Olivia).

## Running the tests

With the environment active, run the full test suite from the repo root:

```bash
conda activate hycom-cice
pytest
```

To run a single file:

```bash
pytest tests/test_phase3.py
```
