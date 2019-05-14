# Test Suite for GAlgebra

## `test_test.py`

It's a 1:1 port from `examples/Terminal/terminal_check.py`, manually change the `print()` calls to assertion of string equality.

## `test_mv.py`

It's only some basic test for class `Mv`.

## `examples/ipython`

We are primarily using [nbval](https://github.com/computationalmodelling/nbval) to verify LaTeX or enhanced console output for GAlgebra by the Jupyter Notebooks in `examples/ipython`.

### How it works

The test notebooks are used as a recorder of outputs (plain text, color ANSI console output, LaTeX), the author of the test notebook will need to visually inspect the outputs in Jupyter or [nbviewer](https://nbviewer.jupyter.org/) to confirm that GAlgebra is working as expected.

Then `nbval` will be run as a plugin of `pytest` in the CI environment for various versions of Python or even of dependencies to ensure the behavior is identical by checking the actual outputs with the recorded outputs. See `.circleci/config.yml` for up-to-date instructions to run them in CI.

While `nbval` is an effective way to write and execute tests for a symbolic library that has multiple output formats, it's also sensitive to subtle changes between versions.

When the tests failed due to these subtle changes, if it should be ignored, the maintainer of test notebooks may utilize `test/.nbval_sanitize.cfg` to sanitize outputs before checking, i.e. replacing outputs by Regular Expressions, or to add a comment like `# NBVAL_IGNORE_OUTPUT` to a specific cell to disable checking for the cell and still keep it in the notebook.

Most of the time, the behavior do change after modifying GAlgebra and some of the tests will fail, the maintainer should carefully review the changed outputs and re-execute the notebook with failed tests and commit the changed notebook. The maintainer should not re-execute and commit the notebook regardless of whether the outputs are correct.

Jupyter Notebooks is actually a JSON file, when it's executed in different OS or Python versions, the output may be identical, but output numbers, Python version info, cell meta data, line endings will definitely change drastically. In order to keep the diff of a commit with modifications to a notebook minimal ( and can be reviewed by reviewers ), the maintainer should always do the following before commiting:

- In Jupyter, restart kernel and execute all cells, so the output numbers will stay linear and the same (if no cells are inserted)
- Shutdown Jupyter, so it won't auto-save and disrupt the following step
- Visually inspect the diff and revert irrelevant changes such as Python version info, cell meta data
- Keep the line endings the same as the original notebook, this can be done by replacing `\r\n"` to `\n"` or vise versa

These steps might look tedious, but one may compare them with manually copying the outputs to the expected output strings of assertions in test, or reviewing such changes to tests.

### Test purposes of each Jupyter Notebook

These notebooks simply execute every example `.py` file (which are converted from old examples) in the corresponding subdirectory of `examples`

- `Terminal.ipynb`:
  - These examples demonstrate various functionalities of GAlgebra
  - output format: color ANSI console output
- `LaTeX.ipynb`
  - Same as the above
  - output format: complete and compileable LaTeX source files
- `Old Format.ipynb`
  - These examples calls deprecated interfaces of GAlgebra

These notebooks are to demonstrate individual scenarios of interest and they are added to tests because they help increasing test coverage of such scenarios :

- `simple_ga_test.ipynb`
- `Smith Sphere.ipynb`
- `second_derivative.ipynb`
- `dop.ipynb`
- `st4.ipynb`
- `colored_christoffel_symbols.ipynb`
- `gr_metrics.ipynb`
- `inner_product.ipynb`
- `verify_doc_python.ipynb`
