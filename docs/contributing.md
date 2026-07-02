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

### Adding or editing pages

- All documentation lives in `docs/` as Markdown files.
- The table of contents is defined in `docs/index.rst`.
- To add a new page, create a `.md` file in `docs/` and add its name (without
  extension) to the appropriate `toctree` block in `docs/index.rst`.
