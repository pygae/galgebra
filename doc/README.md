
# Documentation

This is the diretory containing files to generate the documentation for galgebra ( except for `books/` which contains books about galgebra or Geometric Algebra in general from different sources).

For the generated documentation, please visit https://galgebra.readthedocs.io .

The structure of the `doc` diretory is:

```bash
│                                   # Sphinx Doc Source Files:
│                                   #
├─ index.rst                        # The entry point of the Sphinx doc, it references both
│                                   # galgebra_guide.ipynb and api.rst
├─ galgebra.tex                     # The orignal LaTeX source that generated books/galgebra.pdf
├─ galgebra.md                      # The Markdown source semi-auto converted from galgebra.tex
│                                   # please edit this instead of galgebra.tex (deprecated by it)
│                                   # or galgebra_guide.ipynb (auto generated from it)
├─ galgebra_guide.ipynb             # The Jupyter notebook converted from galgebra.md
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
├─ ipython/                         # Some Jupyter notebooks that are now possibly obsolete
└─ books/                           # Books about galgebra or Geometric Algebra from different sources
```
