=========
Changelog
=========

- :release:`0.5.1 <2024.03.31>`

- :bug:`495` ``MatrixFunction`` is broken since SymPy 1.11, which is required by initializing :class:`~galgebra.ga.Ga` with ``gsym`` and ``coords``, this is now fixed with a workaround (:issue:`507`).

- :bug:`494` There is an extra ``\cdot`` in LaTeX output for multiplication between a numeric power and a polynomial with a leading numeric coefficient since SymPy 1.10. This has to be fixed by SymPy, before that, be cautious when seeing unexpected ``\cdot`` in LaTeX output.

- :feature:`493` SymPy 1.9-1.12 are tested and supported.

- :feature:`492` (also :issue:`473`) Support Python 3.8-3.11, drop support for earlier Python versions.

- :support:`503` The README now supports Github dark mode, and is included by ``sphinx_mdinclude`` in Sphinx with improved support for TOC and math formulas.

- :support:`497` The documentation is now built with the latest Sphinx (7.2.6) and MathJax 3.

- :support:`487` CI has been changed to Github Actions.

- :bug:`468` Hyperbolic functions now requires explicit ``collect`` and ``trigsimp`` to simplify unlike SymPy 1.6.1 and before.

- :support:`476` Make tests more robust by directly testing the printer, rather than ``__str__``.

- :bug:`467` Basis vectors are not normalized in :class:`~galgebra.metric.Metric` with ``norm=True`` when derivatives of basis vectors aren't available.

- :support:`459` Various cleanups and fixes to :class:`~galgebra.lt.Lt`:

  * Construction of ``Lt`` from a linear function is fixed
  * The broken construction of ``Lt`` from a tuple is removed
  * The code of ``Dictionary_to_Matrix`` and ``Matrix_to_dictionary`` is improved
  * Properties of ``Lt`` are tidied, documented, and better tested (:issue:`460`)

- :support:`458` The output of :class:`~galgebra.lt.Lt` objects now use :math:`\mapsto` instead of naming the transformation.

- :support:`457` Line-wrapping of the plaintext printing are disabled by default, as it makes reviewing the output easier, to re-enable it, call ``sympy.init_printing`` with ``wrap_line=True``.

- :support:`455` The static method ``format`` and the attribute ``Lt.mat_fmt`` of :class:`galgebra.lt.Lt` are removed as they were never used.

- :support:`454` The ``l.coords`` and ``l.X`` attributes of :class:`galgebra.lt.Lt` objects are deprecated in favor of using ``l.Ga.coords`` and ``l.Ga.coord_vec``.

- :feature:`454` :class:`galgebra.lt.Lt` objects can now be used with :class:`~galgebra.ga.Ga`\ s that are not constructed with a ``coords`` argument.

- :feature:`450` ANSI printing can now be turned off with ``Eprint(deriv=None, fct=None, base=None)``.

- :support:`450` Many undocumented parts of :func:`galgebra.printer.Eprint` have been removed:

  * The ``debug`` argument to ``Eprint``, which just printed out the other arguments.
  * The ``on`` argument to ``Eprint``, which if ``False`` just made the call do nothing at all.
    Use ``if cond: Eprint(...)`` instead of ``Eprint(..., on=cond)``.
  * The overrideable default settings ``Eprint.defaults``.
  * The constant dictionaries ``Eprint.ColorCode`` and ``Eprint.InvColorCode``.
  * The class methods ``Eprint.Base``, ``Eprint.Fct``, ``Eprint.Deriv``, and ``Eprint.Strip``.
  * The class attributes ``Eprint.base``, ``Eprint.fct``, ``Eprint.deriv``, and ``Eprint.normal``.

  If you were relying on these details to implement your own ANSI printing, it is recommended that you use a package like ``colorama`` instead.

- :bug:`448` Calling :func:`galgebra.printer.Eprint` no longer causes ANSI escape codes to appear in the names of the coefficients of ``ga.mv('A', 'vector')``.

- :feature:`436` The ``Ga`` constructor now allows algebras with only a single basis vector, via a trailing comma in the list of bases.
  This enables algebras like the complex (``Ga('i,', g=[-1])``) and dual (``Ga('delta,', g=[0])``) numbers to be used.

- :release:`0.5.0 <2020.06.05>`

- :bug:`431` Left and right contraction are no longer swapped on scalar :class:`~sympy.core.expr.Expr` instances.
- :feature:`232` Python 3.8 is now supported and tested.
- :feature:`212` Python 2.7 is no longer supported, allowing more features from Python 3.x to be used.

- :feature:`0` Generally improved documentation, including:

  * The re-inclusion of missing code examples that were already referenced from the text (:issue:`293`).
  * Fixed footnote and citation references (:issue:`287`, :issue:`290`, :issue:`291`)
  * Faster and smaller pages with less MathJax loading time (:issue:`292`, :issue:`296`).
  * An adapted version of a notebook tutorial from Alan Bromborsky is now present in the docs, :doc:`tutorials/algebra` (:issue:`375`).
  * Type annotations on many functions in the :doc:`API docs <api>` (:issue:`305`).
  * Correction of curly quotes that led to invalid python syntax (:issue:`294`).

- :support:`393` ``galgebra.dop.Sdop.str_mode`` has been removed, as it was always ``False``.

- :bug:`391` Subtracting a scalar from a :class:`~galgebra.mv.Mv` instance no longer fails.

- :bug:`387` :func:`galgebra.printer.Print_Function` no longer emits invalid LaTeX with an unwanted ``\\``

- :support:`386` The ``title`` attribute of :class:`~galgebra.mv.Dop` has been removed. as it was always ``None``.

- :support:`385` ``galgebra.printer.print_replace`` has been removed.

- :support:`384` The ``dop`` argument to :func:`~galgebra.printer.Format`, along with the corresponding static member of :class:`~galgebra.printer.LatexPrinter`, has been removed, as had no effect.

- :bug:`382` The result of calling :func:`galgebra.printer.latex` on a multi-vector when the global ``galgebra_mv_fmt`` setting is not 1 is now valid to use within math mode.

- :bug:`378` (also :issue:`379`) The post-processing performed by :func:`galgebra.printer.tex` is now better-behaved:

  * It no longer translates ``grad`` into :math:`\nabla` if it is part of a longer word like ``gradual``.
  * It no longer mangles occurrences of ``@@``.
  * It no longer throws away ``#`` and ``%`` symbols that appear in the middle of lines.
  * It no longer injects stale information from previous lines into subsequent lines.
  * It now strips all leading and trailing whitespace before wrapping in ``equation`` or ``align``.

- :support:`376` (also :issue:`258`, :issue:`371`) The ``repr()`` and ``str()`` of both :mod:`galgebra` and :mod:`sympy` objects will no longer ever return latex strings, instead always returning plaintext strings. If you want latex strings, use :func:`galgebra.printer.latex`.
  Note that this does not affect support for ``print(sympy_object)``, which continues to print latex representations provided :func:`galgebra.printer.Format` has been called.

- :bug:`372` For scalar multivectors, the printed result is now mathematically equivalent in plaintext and latex mode.

- :bug:`369` (also :issue:`380`) The ``Fmt`` method of :meth:`Mv <galgebra.mv.Mv.Fmt>`, :meth:`Dop <galgebra.mv.Dop.Fmt>`, :meth:`Lt <galgebra.lt.Lt.Fmt>`, :meth:`Mv <galgebra.lt.Mlt.Fmt>` now works properly in both IPython (giving plaintext output) and Jupyter (giving LaTeX output).
- :bug:`369` The ``fmt`` argument fo the ``Fmt`` method of :meth:`Mv <galgebra.mv.Mv.Fmt>`, :meth:`Dop <galgebra.mv.Dop.Fmt>`, :meth:`Lt <galgebra.lt.Lt.Fmt>`, :meth:`Mlt <galgebra.lt.Mlt.Fmt>` no longer has side effects on subsequent ``print`` statements.

- :support:`367` The ``fmt_dop`` argument to and ``dop_fmt`` attribute of :class:`galgebra.mv.Dop` have been removed, as they had no effect.
- :support:`367` (also :issue:`364`, :issue:`369`) The following properties of :class:`~galgebra.printer.GaLatexPrinter` and :class:`~galgebra.printer.GaPrinter` have been removed:

  * ``fmt``
  * ``prev_fmt``
  * ``dop_fmt``
  * ``prev_dop_fmt``
  * ``lt_fmt``
  * ``prev_lt_fmt``
  * ``Dmode``
  * ``Fmode``
  * ``inv_trig_style``
  * ``dop``

- :feature:`364` The galgebra printer now respects the global ``inv_trig_style`` sympy print setting.
- :feature:`364` The galgebra print customizations ``Dmode``, ``Fmode``, ``fmt``, and ``Mlt.lcnt`` can now be set via::
    
    sympy.init_printing(
        omit_function_args=True,  # Fmode
        omit_partial_derivative_fraction=True,  # Dmode
        galgebra_mv_fmt=1,  # fmt
        galgebra_mlt_lcnt=1,  # lcnt
    )

  These names are provisional, and may change in future.

- :support:`358` ``galgebra.printer.find_executable`` has been removed, as the functionality is provided by :func:`shutil.which`.

- :bug:`354` :func:`galgebra.printer.oprint` no longer strips the last ``)`` from :class:`sympy.Matrix` objects within lists.

- :support:`349` :func:`galgebra.printer.Get_Program` is deprecated. This
  function used to have to be called before :func:`galgebra.printer.Print_Function`,
  but now all it does is act as an enable/disable flag for the latter, such that ``Get_Program(True); Print_Function()`` was a no-op.

- :bug:`348` :func:`galgebra.printer.Print_Function` no longer emits code after the function body. Previously, examples had a ``def dummy(): pass`` function to work around this.

- :bug:`345` :attr:`~galgebra.ga.Ga.e_sq` no longer contains a floating point term.

- :feature:`336` The scalar product :math:`a * b` is now available via :meth:`galgebra.ga.Ga.scalar_product`.
  Note that there is no operator overload on multivectors yet.

- :bug:`323` (also :issue:`377`, :issue:`340`) Many non-``Mv`` objects now render better in the default sympy printer. This bug manifested most often when using ``Mv.obj``, and would result in misrenderings like :math:`e^e_{xy}` when :math:`e_x \wedge e_y` was intended, or :math:`e_{x.y}` when :math:`e_x \cdot e_y` was intended.

- :bug:`320` :class:`galgebra.lt.Mlt` no longer crashes at construction, arithmetic, or multiplication.
- :support:`320` The following attributes of :class:`galgebra.ga.Ga` have been removed:

  - `a`

- :bug:`319` :meth:`galgebra.mv.Mv.get_coefs` now returns ``0`` in the place of empty coefficients.
- :bug:`319` :meth:`galgebra.ga.Ga.make_grad` no longer has a broken cache that ignores ``cmpflg`` if both the left and right gradient operator are requested for the same vector
- :bug:`319` :meth:`galgebra.ga.Ga.make_grad` no longer crashes when called on an algebra with no :attr:`~galgebra.ga.Ga.coords`.

- :support:`286` The :mod:`galgebra.deprecated` module, which provides old-style APIs, now emits :exc:`DeprecationWarning`\ s.

- :support:`284` (also :issue:`317`) :meth:`galgebra.lt.Lt.setup` has been deprecated, along with the attributes it populated :attr:`galgebra.ga.Ga.lt_x` (the same as :attr:`~galgebra.ga.Ga.coords_vec`) and :attr:`galgebra.ga.Ga.lt_coords` (the same as :attr:`~galgebra.ga.Ga.coords`).

- :bug:`264` :class:`~galgebra.mv.Dop` no longer emits latex strings when printed in non-latex mode.

- :support:`259` ``galgebra.printer.oprint`` now aligns results in columns

- :bug:`258` The result of simplifying sympy expressions is no longer dependent on whether :func:`galgebra.printer.Format` has been called.

- :support:`252` (also :issue:`310`) The ``inverse_metric()`` and ``derivatives_of_g()`` methods of :class:`~galgebra.metric.Metric` are deprecated, as the properties they computed (:attr:`~galgebra.metric.Metric.g_inv` and :attr:`~galgebra.metric.Metric.dg`) are now computed automatically.

- :feature:`252` (also :issue:`310`) Many attributes of :class:`~galgebra.ga.Ga` instances are now _always_ present, rather than being conditionally present depending on the type of algebra. These include:

  * :attr:`~galgebra.ga.Ga.r_basis`
  * :attr:`~galgebra.ga.Ga.r_basis_mv`
  * :attr:`~galgebra.ga.Ga.r_basis_dict`
  * :attr:`~galgebra.metric.Metric.g_inv`
  * :attr:`~galgebra.metric.Metric.g_adj`

  Other properties that are still only meaningful for some algebras now exist but raise clearer error messages:

  * :attr:`~galgebra.ga.Ga.coord_vec`
  * :attr:`~galgebra.ga.Ga.bases`
  * :attr:`~galgebra.ga.Ga.basic_mul`

  These lists are not exhaustive.

- :feature:`243` :meth:`galgebra.mv.Mv.subs` now accepts all the same arguments as :func:`sympy.subs`.

- :support:`216` ``galgebra.metric.test_init_slots`` has been removed. The functionality this provided is superseded by the language feature of keyword-only arguments.

- :support:`202` :class:`~galgebra.dop.Pdop` and :class:`~galgebra.dop.Sdop` instance are no longer associated with a Ga. As a result, their ``.Ga`` attribute no longer exists, and the :meth:`~galgebra.ga.Ga.pdop` and :meth:`~galgebra.ga.Ga.sdop` methods of :class:`~galgebra.ga.Ga` are deprecated in favor of calling the constructors directly.
  For Ga-aware operators, continue to use :class:`~galgebra.mv.Dop`.


- :bug:`53` Calling :func:`sympy.sympify` (or any other sympy function) on a :class:`~galgebra.mv.Mv` instance no longer raises :exc:`RecursionError`, and instead raises :exc:`TypeError` with a helpful message.

- :bug:`344` :func:`galgebra.metric.collect` no longer discards terms that were not requested.

- :support:`338` The undocumented and mispelt static method :meth:`galgebra.lt.Mlt.extact_basis_indexes` (which just computed a value equivalent to :attr:`galgebra.ga.Ga.basis_super_scripts`) has been deprecated.

- :feature:`326` The :meth:`galgebra.ga.Ga.make_grad` function now accepts multivectors, not just vectors.

- :support:`335` The undocumented method :meth:`galgebra.ga.Ga.X` (the same as :attr:`~galgebra.ga.Ga.coords_vec`) has been deprecated.


- :support:`280` The ``galgebra.utils`` module, which provided Python 2 compatibility helpers, has been deprecated.

- :support:`245` (also :issue:`248`, :issue:`334`) The following attributes have been deprecated to reduce the number of similar members in :class:`~galgebra.ga.Ga`.

  * Unified into the :attr:`~galgebra.ga.ProductFunction.table_dict` attribute of :class:`~galgebra.ga.ProductFunction` instances:

    * :attr:`galgebra.ga.Ga.wedge_table_dict` |rarr| ``wedge.table_dict``
    * :attr:`galgebra.ga.Ga.left_contract_table_dict` |rarr| ``left_contract.table_dict``
    * :attr:`galgebra.ga.Ga.right_contract_table_dict` |rarr| ``right_contract.table_dict``
    * :attr:`galgebra.ga.Ga.dot_table_dict` |rarr| ``hestenes_dot.table_dict``
    * :attr:`galgebra.ga.Ga.mul_table_dict` |rarr| ``mul.table_dict``
    * :attr:`galgebra.ga.Ga.basic_mul_table_dict` |rarr| ``basic_mul.table_dict``

  * Unified into methods of :class:`~galgebra.ga.ProductFunction` instances:

    * :attr:`galgebra.ga.Ga.non_orthogonal_bases_products` |rarr| ``basic_mul.of_basis_bases``
    * :attr:`galgebra.ga.Ga.wedge_product_basis_blades` |rarr| ``wedge.of_basis_blades``
    * :attr:`galgebra.ga.Ga.geometric_product_basis_blades` |rarr| ``mul.of_basis_blades``
    * :attr:`galgebra.ga.Ga.dot_product_basis_blades` and :attr:`galgebra.ga.Ga.non_orthogonal_dot_product_basis_blades` |rarr| ``<some_dot_type>.of_basis_blades``.
      In future it will no longer be possible to call the non-ortho method for an orthogonal ga, or vice-versa.

  * Unified into :attr:`galgebra.ga.Ga.blade_expansion_dict`:

    * :attr:`galgebra.ga.Ga.blade_expansion` |rarr| ``blade_expansion_dict.items()``

  * Unified into :attr:`galgebra.ga.Ga.base_expansion_dict`:

    * :attr:`galgebra.ga.Ga.base_expansion` |rarr| ``base_expansion_dict.items()``

  * Unified into :attr:`galgebra.ga.Ga.basic_mul_table_dict`:

    * :attr:`galgebra.ga.Ga.basic_mul_table` |rarr| ``basic_mul_table_dict.items()``
    * :attr:`galgebra.ga.Ga.basic_mul_keys` |rarr| ``basic_mul_table_dict.keys()``
    * :attr:`galgebra.ga.Ga.basic_mul_values` |rarr| ``basic_mul_table_dict.values()``

  * Unified into :attr:`galgebra.ga.Ga.indexes_to_bases_dict`:

    * :attr:`galgebra.ga.Ga.indexes_to_bases` |rarr| ``indexes_to_bases_dict.items()``
    * :attr:`galgebra.ga.Ga.bases_to_indexes` |rarr| ``indexes_to_bases_dict.inverse.items()``
    * :attr:`galgebra.ga.Ga.bases_to_indexes_dict` |rarr| ``indexes_to_bases_dict.inverse``

  * Unified into :attr:`galgebra.ga.Ga.indexes_to_bases_dict`:

    * :attr:`galgebra.ga.Ga.indexes_to_blades` |rarr| ``indexes_to_blades_dict.items()``
    * :attr:`galgebra.ga.Ga.blades_to_indexes` |rarr| ``indexes_to_blades_dict.inverse.items()``
    * :attr:`galgebra.ga.Ga.blades_to_indexes_dict` |rarr| ``indexes_to_blades_dict.inverse``

- :support:`64` The following attributes have been deprecated in favor of using the new :attr:`~galgebra.ga.GradedTuple.flat` attribute:

  * :attr:`galgebra.ga.Ga.indexes_lst` |rarr| :attr:`galgebra.ga.Ga.indexes` ``.flat``
  * :attr:`galgebra.ga.Ga.bases_lst` |rarr| :attr:`galgebra.ga.Ga.bases` ``.flat``
  * :attr:`galgebra.ga.Ga.blades_lst` |rarr| :attr:`galgebra.ga.Ga.blades` ``.flat``
  * :attr:`galgebra.ga.Ga.mv_blades_lst` |rarr| :attr:`galgebra.ga.Ga.mv_blades` ``.flat``

  Since the replacement properties also include the scalar element, ``ga.blades_lst`` is equivalent to ``ga.blades.flat[1:]`` (and likewise for the other properties).

- :release:`0.4.5 <2019.12.31>`
- :support:`83` This is the last release to support Python 2.7.
- :feature:`195` Deprecate confusing parts of the public API:

  * :meth:`galgebra.ga.Ga.mv_I`, which was a bad way to spell :meth:`galgebra.ga.Ga.E`.
  * :meth:`galgebra.ga.Ga.mv_x`, which was a bad way to create a vector with a weird name.
  * ``Pdop(None)``, which is better spelt ``Pdop({})``, the latter of which has always worked (:issue:`187`).
  * :meth:`galgebra.ga.Ga.Pdop_identity`, which is better spelt ``ga.pdop({})`` (:issue:`194`).
  * :meth:`galgebra.ga.Ga.Pdiffs`, which is better spelt ``ga.pdop(x)`` (:issue:`194`).
  * :meth:`galgebra.ga.Ga.sPds`, which is better spelt ``ga.sdop(x)`` (:issue:`194`).

- :bug:`188` :meth:`~galgebra.ga.Ga.pdop` and :meth:`~galgebra.ga.Ga.sdop` are now transparent aliases to :class:`~galgebra.mv.Pdop` and :class:`~galgebra.mv.Sdop`, respectively.
  Previously, they would crash.
- :bug:`188` Passing a ``ga`` keyword-argument to :meth:`~galgebra.ga.Ga.sm`, :meth:`~galgebra.ga.Ga.mv`, and :meth:`~galgebra.ga.Ga.lt` is now an error, previously it was silently ignored.
- :bug:`187` ``Sdop(x)`` now functions as intended as an alias of ``Sdop([(S(1), Pdop({x: 1}))])``, which previously crashed.
- :bug:`187` :class:`~galgebra.mv.Pdop` no longer silently accepts extra arguments.
- :bug:`183` ``op1 == op2`` now works correctly for :class:`~galgebra.mv.Sdop` and :class:`~galgebra.mv.Dop` instances.
- :bug:`157` Adding an :class:`~galgebra.mv.Sdop` instance to a scalar now works as intended, rather than crashing.
- :bug:`151` ``Dop([], ga=ga)`` and ``Sdop([], ga=ga)`` now evaluate to multiplication by zero, not by one.
  Multiplication by one can as always be spelt ``Dop([(S(1), Pdop({}, ga=ga}))], ga=ga)``.
- :bug:`177` :class:`~galgebra.mv.Dop` objects that evaluate to ``0`` no longer raise cryptic ``ValueError``\ s when operated on.
- :feature:`172` :data:`galgebra.__version__` has been added, which contains the version string.
- :feature:`164` (and :issue:`169`, :issue:`170`) Sympy 1.5 is officially supported and tested.
- :support:`167` Python 3.4 is no longer supported.
- :bug:`165` :func:`galgebra.metric.linear_expand` no longer accepts a mode argument, as this did not work properly.
  For the old behavior of ``linear_expand(x, mode=True)``, use ``linear_expand_terms(x)`` instead.
- :bug:`151` (also :issue:`150`) :class:`~galgebra.mv.Dop`, :class:`~galgebra.mv.Sdop`, and :class:`~galgebra.mv.Pdop` no longer have mutating methods.
  This fixed issues where not only would the laplacian be sometimes calculated incorrectly, but its correctness would vary depending on whether it had been printed!
- :bug:`134` :attr:`~galgebra.ga.Ga.dot_table_dict` now contains correct values (zero) for scalar keys
- :bug:`90` :attr:`galgebra.ga.Ga.blades`, :attr:`galgebra.ga.Ga.bases`, and :attr:`galgebra.ga.Ga.indices` now reference the scalar ``S(0)`` as the single grade-0 object. Previously they listed no grade 0 objects.
- :bug:`81` (also :issue:`180`) Passing coefficients as ``Mv(coefs, 'odd', ga=ga)`` is forbidden.
- :bug:`80` (also :issue:`57`, :issue:`58`, :issue:`97`) The :class:`galgebra.mv.Mv` constructor no longer silently accepts illegal arguments, and produces better error messages.
- :feature:`78` :meth:`~galgebra.ga.Ga.grads` now raises a better error when it fails, and is faster.
- :support:`72` Other internal cleanup.
- :feature:`66` Remove unused or unusable code in the public API:

  * ``Ga.mul_table``, ``Ga.wedge_table``, ``Ga.dot_table``, ``Ga.left_contract_table``,
    and ``Ga.right_contract_table``, all of which were the empty list, ``[]``.
  * ``galgebra.mv.modules``, a string which served no purpose (:issue:`71`).
  * ``__add_ab__``, ``__sub_ab__``, ``__mul_ab__``, and ``__div_ab__``, none of are real magic method names (:issue:`67`).
    No code should be calling these directly anyway.
  * ``Dop.flatten_one_level`` which is better spelt ``itertools.chain.from_iterable`` (:issue:`175`).
  * ``Dop.basic``, which was a non-working version of :meth:`~galgebra.ga.Ga.grads()` (:issue:`185`).
  * ``Pdop.setGa``, ``Sdop.setGa``, ``Dop.setGa``, all of which rely on dangerous global state - use the ``ga`` keyword argument of the constructors of these types instead (:issue:`163`).

- :feature:`66` The :attr:`~galgebra.ga.Ga.mul_table_dict` table, and the equivalent tables for the other products, are now computed lazily when indexed. These are now all documented too.
- :bug:`61` Make contraction and Hestenes dot products thread-safe.
  Previously these relied on the :attr:`~galgebra.ga.Ga.dot_mode` setting not being changed mid-operation.
  The :meth:`~galgebra.ga.Ga.dot` method still respects this setting, but is no longer used internally.
- :bug:`60` (also :issue:`141`) Make the following operations on :class:`galgebra.mv.Mv` non-mutating:

  * :meth:`~galgebra.mv.Mv.blade_rep`
  * :meth:`~galgebra.mv.Mv.base_rep`
  * :meth:`~galgebra.mv.Mv.diff`
  * :meth:`~galgebra.mv.Mv.simplify`
  * :meth:`~galgebra.mv.Mv.expand`
  * :meth:`~galgebra.mv.Mv.collect`
  * ``print(mv)``

  Any code relying on this behavior will need to change from ``x.method()`` to ``x = x.method()``.
  Note that the latter syntax was always supported even before this change.

- :support:`59` (also :issue:`65`) Make internal attributes and helper functions private.
  All of the following have been made private by being renamed:

  * ``Mv.make_grad`` (:issue:`59`).
  * ``Mv.make_scalar`` (:issue:`59`).
  * ``Mv.make_vector`` (:issue:`59`).
  * ``Mv.make_bivector`` (:issue:`59`).
  * ``Mv.make_pseudo_scalar`` (:issue:`59`).
  * ``Mv.make_multivector`` (:issue:`59`).
  * ``Mv.make_spinor`` (:issue:`59`).
  * ``Mv.make_odd`` (:issue:`59`).
  * ``Ga.build_bases()`` (:issue:`65`).
  * ``Ga.basis_product_tables()`` (:issue:`65`).
  * ``Ga.build_reciprocal_basis()`` (:issue:`65`).
  * ``Ga.build_connection()`` (:issue:`65`).
  * ``Ga.non_orthogonal_mul_table()`` (:issue:`65`).
  * ``Ga.base_blade_conversions()`` (:issue:`65`).
  * ``Ga.init_connect_flg()`` (:issue:`65`).
  * ``Ga.derivatives_of_basis()`` (:issue:`65`).
  * ``Ga.lt_flg`` (:issue:`65`).
  * ``Ga.agrads`` (:issue:`65`).
  * ``Ga.dbases`` (:issue:`65`).
  * ``Ga.XOX`` (:issue:`195`).

- :support:`55` Rename ``*kargs`` to ``*args`` internally, to match convention.
  This has no effect on callers, but makes the docs and source easier to read.
- :feature:`50` (also :issue:`51`, :issue:`56`) Improve documentation formatting:

  * LaTeX and code samples are now appropriately formatted
  * Attributes of classes now have permalinks

- :support:`46` (also :issue:`69`) Remove unnecessary executable bit from importable python files, and the corresponding no-op code that would be run.

- :release:`0.4.4 <2019.09.30>`
- :feature:`17` Fix examples under both Python 2 & 3

  * Fix `examples/*` and verify them in CI using `nbval`
  * Test coverage increased from 48.89% to 66.52%
  * Add README for `test` and `examples`

- :feature:`9` Documentation now available at https://galgebra.readthedocs.io/

  * Convert doc to Sphinx with the help of `pandoc`, `notedown <https://github.com/aaren/notedown>`_  and `nbsphinx <https://nbsphinx.readthedocs.io/en/0.3.5/>`_
  * Add `Getting Started` guide to README
  * Update installation instructions in README
  * Add migration guide from `sympy.galgebra` and `brombo/galgebra`
  * Add Changelog
  * Add doc for examples, tests and bundled resources

- :bug:`15` Fix printing of some products and inverses of multivectors
- :bug:`18` Fix TypeError of unicode string, use `BytesIO` instead of `StringIO`
- :bug:`26` Fix calculation of the Christoffel symbols
- :bug:`27` Fix broken class MV
- :bug:`29` Fix that sometimes plain text is printed with or instead of LaTeX in Jupyter Notebooks
- :bug:`30` Fix bugs of using LaTeX as symbol name
- :bug:`32` Fix calculation of reciprocal basis for non-orthogonal basis
- :bug:`31` Freeze the depended version of SymPy to 1.3
- :support:`17` Setup Circle CI to build more Python versions faster

  * TravisCI build for PRs is now removed

- :release:`0.4.3 <2018.02.18>`
- :feature:`2` Support Python 3

  * Convert code to be Python 3 compatible
  * Pass tests under both Python 2 & 3
  * Support installing from PyPI instead of setting `pth`
  * Support importing with `from galgebra.<package name> import *`

- :support:`2` Setup Travis CI
- :support:`8` Add test coverage in CI using using `pytest <https://pytest.org/>`_ and `CodeCov <https://codecov.io/gh/pygae/galgebra>`_
- :support:`8` Validate existing Jupyter notebooks using `nbval <https://github.com/computationalmodelling/nbval>`_
- :support:`8` Remove NumPy dependency
- :support:`2` Remove .pyc files and add standard .gitignore for python
- :support:`2` Clean up obsolete code in setup.py and make it simple
- :bug:`2` Fixes `brombo/galgebra#19 <https://github.com/brombo/galgebra/issues/19>`_: Failures in `test_noneuclidian_distance_calculation`
- :bug:`2` Fixes `brombo/galgebra#20 <https://github.com/brombo/galgebra/issues/20>`_: Incorrect LaTeX output in `test_derivatives_in_spherical_coordinates`

.. include:: <isonum.txt>
