# Development memo

## Installation

```shell
git clone git@github.com:spglib/spinspg.git
cd spinspg
conda create -y -n spinspg python=3.10 pip
conda activate spinspg
pip install -e ".[dev,docs]"
pre-commit install
```

## Compile documents

```shell
sphinx-autobuild docs docs_build
# open localhost:8000 in your browser
```

## Release

```shell
# Confirm the version number via `setuptools-scm`
python -m setuptools_scm

# Update changelog here
vim docs/changelog.md

# Push with tag
git tag <next-version>
git push origin main
git push origin <next-version>
```
