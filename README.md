# Computable Model of Chemistry Textbook

## Install

```
uv venv .venv --python 3.12.0
source .venv/bin/activate
uv pip install -e '.[dev]'
uv pip install -e '.[doc]'
```

## Build Sphinx Document

From root:

```
make -C docs clean html
sphinx-autobuild docs/source docs/build/html
```
