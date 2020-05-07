
# Documentation

This is the directory containing files to generate the documentation for galgebra (except for `books/`, which contains books about `galgebra` or Geometric Algebra in general from different sources).

For the generated documentation, please visit https://galgebra.readthedocs.io .

The structure of the `doc` diretory is:

```bash
│                                   # Sphinx Doc Source Files:
│                                   #
├─ index.rst                        # The entry point of the Sphinx doc, it references both
│                                   # galgebra_guide.ipynb and api.rst
├─ galgebra.tex                     # The original LaTeX source that generated books/galgebra.pdf
├─ galgebra_guide.rst               # The convert vertion of the latex document.
│                                   # This was converted to markdown manually, then to rST via nodedown and nbsphinx.
├─ api.rst                          # Use automodule to extract doc from galgebra Python source files
│
│                                   # Configurations:
│                                   #
├─ readthedocs-pip-requirements.txt # PIP requirements.txt for readthedocs
├─ conf.py                          # Configurations for Sphinx
│
│                                   # Scripts:
│                                   #
├─ Makefile                         # So `make html` can generate the doc
├─ make.bat                         # Do the same as above for Windows
│
│                                   # Directories Used By Sphinx:
│                                   #
├─ _static/                         # Configured html_static_path for Sphinx, see conf.py
├─ _templates/                      # Configured templates_path for Sphinx, see conf.py
├─ images/                          # Images used in doc
│
│                                   # Legacy Directories:
│                                   #
├─ python/                          # Some Python scripts referenced in the doc but not all of them works
└─ books/                           # Books about galgebra or Geometric Algebra from different sources
```
