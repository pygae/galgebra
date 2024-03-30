Representations in *sympy*
==========================

.. math::
   \newcommand{\bm}[1]{\boldsymbol{#1}}
   \newcommand{\ubh}{\bm{\hat{u}}}
   \newcommand{\ebh}{\bm{\hat{e}}}
   \newcommand{\ebf}{\bm{e}}
   \newcommand{\mat}[1]{\left [ {#1} \right ]}
   \newcommand{\bra}[1]{{#1}_{\mathcal{G}}}
   \newcommand{\ket}[1]{{#1}_{\mathcal{D}}}
   \newcommand{\ds}{\displaystyle}
   \newcommand{\bfrac}[2]{\displaystyle\frac{#1}{#2}}
   \newcommand{\lp}{\left (}
   \newcommand{\rp}{\right )}
   \newcommand{\half}{\frac{1}{2}}
   \newcommand{\llt}{\left <}
   \newcommand{\rgt}{\right >}
   \newcommand{\abs}[1]{\left |{#1}\right |}
   \newcommand{\pdiff}[2]{\bfrac{\partial {#1}}{\partial {#2}}}
   \newcommand{\pdifftwo}[3]{\bfrac{\partial^{2} {#1}}{\partial {#2}\partial {#3}}}
   \newcommand{\lbrc}{\left \{}
   \newcommand{\rbrc}{\right \}}
   \newcommand{\set}[1]{\lbrc {#1} \rbrc}
   \newcommand{\W}{\wedge}
   \newcommand{\R}{\dagger}
   \newcommand{\lbrk}{\left [}
   \newcommand{\rbrk}{\right ]}
   \newcommand{\com}[1]{\lbrk {#1} \rbrk}
   \newcommand{\proj}[2]{\llt {#1} \rgt_{#2}}
   %\newcommand{\bm}{\boldsymbol}
   \newcommand{\braces}[1]{\left \{ {#1} \right \}}
   \newcommand{\grade}[1]{\left < {#1} \right >}
   \newcommand{\f}[2]{{#1}\lp {#2} \rp }
   \newcommand{\paren}[1]{\lp {#1} \rp }
   \newcommand{\eval}[2]{\left . {#1} \right |_{#2}}
   \newcommand{\prm}[1]{{#1}'}
   \newcommand{\ddt}[1]{\bfrac{d{#1}}{dt}}
   \newcommand{\deriv}[3]{\bfrac{d^{#3}#1}{d{#2}^{#3}}}
   \newcommand{\be}{\begin{equation}}
   \newcommand{\ee}{\end{equation}}
   \newcommand{\eb}{\bm{e}}
   \newcommand{\ehb}{\bm{\hat{e}}}
   \newcommand{\Tn}[2]{\f{\mathcal{T}_{#2}}{#1}}
   \newcommand{\tr}{\mbox{tr}}
   \newcommand{\T}[1]{\texttt{#1}}
   \newcommand{\grd}{\bm{\nabla}}
   \newcommand{\indices}[1]{#1}
   \newcommand{\xRightarrow}[1]{\overset{#1}{\Rightarrow}}

Representation of Multivectors
------------------------------

The *sympy* python module offers a simple way of representing multivectors using linear combinations of commutative expressions (expressions consisting only of commuting *sympy* objects) and non-commutative symbols. We start by defining :math:`n` non-commutative *sympy* symbols as a basis for the vector space

``e_1, ..., e_n = symbols('e_1,...,e_n', commutative=False, real=True)``

Several software packages for numerical geometric algebra calculations are available from Doran-Lasenby group and the Dorst group. Symbolic packages for Clifford algebra using orthogonal bases such as :math:`{{\eb}}_{i}{{\eb}}_{j}+{{\eb}}_{j}{{\eb}}_{i} = 2\eta_{ij}`, where :math:`\eta_{ij}` is a numeric array are available in Maple and Mathematica. The symbolic algebra module, *ga*, developed for python does not depend on an orthogonal basis representation, but rather is generated from a set of
:math:`n` arbitrary symbolic vectors :math:`{{\eb}}_{1},{{\eb}}_{2},\dots,{{\eb}}_{n}` and a symbolic metric tensor :math:`g_{ij} = {{\eb}}_{i}\cdot {{\eb}}_{j}` (the symbolic metric can be symbolic constants or symbolic function in the case of a manifold).

In order not to reinvent the wheel all scalar symbolic algebra is handled by the python module *sympy* and the abstract basis vectors are encoded as non-commuting *sympy* symbols.

The basic geometric algebra operations will be implemented in python by defining a geometric algebra class, *Ga*, that performs all required geometric algebra an calculus operations on *sympy* expressions of the form (Einstein summation convention)

.. math:: \be F +\sum_{r=1}^{n}F^{i_{1}\dots i_{r}}\eb_{i_{1}}\dots\eb_{i_{r}} \ee

where the :math:`F`\ ’s are *sympy* symbolic constants or functions of the coordinates and a multivector class, *Mv*, that wraps *Ga* and overloads the python operators to provide all the needed multivector operations as shown in Table :ref:`ops` where :math:`A` and :math:`B` are any two multivectors (In the case of :math:`+`, :math:`-`, :math:`*`, :math:`{\wedge}`, :math:`|`, :math:`<`, and :math:`>` the operation is also defined if :math:`A` or :math:`B` is a *sympy* symbol or a *sympy* real
number).

.. _ops:

.. table:: Operators

   ================== =================================
   :math:`A+B`        sum of multivectors
   :math:`A-B`        difference of multivectors
   :math:`A*B`        geometric product of multivectors
   :math:`A{\wedge}B` outer product of multivectors
   :math:`A{\vert}B`  inner product of multivectors
   :math:`A{<}B`      left contraction of multivectors
   :math:`A{>}B`      right contraction of multivectors
   :math:`A{/}B`      division of multivectors
   ================== =================================

Multivector operations for GA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since ``<`` and ``>`` have no r-forms (in python for the ``<`` and ``>`` operators there are no ``__rlt__()`` and ``__rgt__()`` member functions to overload) we can only have mixed modes (sympy scalars and multivectors) if the first operand is a multivector.

Except for ``<`` and ``>`` all the multivector operators have r-forms so that as long as one of the operands, left or right, is a multivector the other can be a multivector or a scalar (*sympy* symbol or number).

Operator Precedence
~~~~~~~~~~~~~~~~~~~

**Note that** the operator order precedence is determined by python and is not necessarily that used by geometric algebra. It is *absolutely essential* to use parenthesis in multivector expressions containing ``^``, ``|``, ``<``, and/or ``>``. As an example let ``A`` and ``B`` be any two multivectors. Then ``A + A*B = A +(A*B)``, but ``A+A^B = (2*A)^B`` since in python the ``^`` operator has a lower precedence than the ``+`` operator. In geometric algebra the outer and inner products and the
left and right contractions have a higher precedence than the geometric product and the geometric product has a higher precedence than addition and subtraction. In python the ``^``, ``|``, ``>``, and ``<`` all have a lower precedence than ``+`` and ``-`` while ``*`` has a higher precedence than ``+`` and ``-``.

**Additional care has to be used** when using the operators ``!=`` and ``==`` with the operators ``<`` and ``>``. All these operators have the same precedence and are evaluated chained from left to right. To be completely safe for expressions such as ``A == B`` or ``A != B`` always user ``(A) == (B)`` and ``(A) != (B)`` if ``A`` or ``B`` contains a left, ``<``, or right, ``>``, contraction.

For those users who wish to define a default operator precedence the functions ``def_prec()`` and ``GAeval()`` are available in the module printer.

.. autofunction:: galgebra.printer.def_prec
   :noindex:

.. autofunction:: galgebra.printer.GAeval
   :noindex:

.. _BasisMetric:

Vector Basis and Metric
-----------------------

The two structures that define the ``metric`` class (inherited by the geometric algebra class) are the symbolic basis vectors and the symbolic metric. The symbolic basis vectors are input as a string with the symbol name separated by spaces. For example if we are calculating the geometric algebra of a system with three vectors that we wish to denote as ``a0``, ``a1``, and ``a2`` we would define the string variable:

.. code:: python

   basis = 'a0 a1 a2'

that would be input into the geometric algebra class instantiation function, ``Ga()``. The next step would be to define the symbolic metric for the geometric algebra of the basis we have defined. The default metric is the most general and is the matrix of the following symbols

.. math::

   \begin{equation}\label{metric}
     g = \lbrk
     \begin{array}{ccc}
       (a0.a0)   & (a0.a1)  & (a0.a2) \\
       (a0.a1) & (a1.a1)  & (a1.a2) \\
       (a0.a2) & (a1.a2) & (a2.a2) \\
     \end{array}
     \rbrk
     \end{equation}

where each of the :math:`g_{ij}` is a symbol representing all of the dot products of the basis vectors. Note that the symbols are named so that :math:`g_{ij} = g_{ji}` since for the symbol function :math:`(a0.a1) \ne (a1.a0)`.

Note that the strings shown in the above equation are only used when the values of :math:`g_{ij}` are output (printed). In the ga module (library) the :math:`g_{ij}` symbols are stored in a member of the geometric algebra instance so that if ``o3d`` is a geometric algebra then ``o3d.g`` is the metric tensor ( :math:`g_{ij} =` ``o3d.g[i, j]``) for that algebra.

The default definition of :math:`g` can be overwritten by specifying a string that will define :math:`g`. As an example consider a symbolic representation for conformal geometry. Define for a basis

.. code:: python

   basis = 'a0 a1 a2 n nbar'

and for a metric

.. code:: python

   g = '# # # 0 0, # # # 0 0, # # # 0 0, 0 0 0 0 2, 0 0 0 2 0'

then calling ``cf3d = Ga(basis, g=g)`` would initialize the metric tensor

.. math::

   \be g = \lbrk\begin{array}{ccccc}
       (a0.a0) & (a0.a1)  & (a0.a2) & 0 & 0\\
       (a0.a1) & (a1.a1)  & (a1.a2) & 0 & 0\\
       (a0.a2) & (a1.a2)  & (a2.a2) & 0 & 0 \\
       0 & 0 & 0 & 0 & 2 \\
       0 & 0 & 0 & 2 & 0
     \end{array}
     \rbrk \ee

for the ``cf3d`` (conformal 3-d) geometric algebra.

Here we have specified that ``n`` and ``nbar`` are orthogonal to all the ``a``\ ’s, ``(n.n) = (nbar.nbar) = 0``, and ``(n.nbar) = 2``. Using ``#`` in the metric definition string just tells the program to use the default symbol for that value.

When ``Ga`` is called multivector representations of the basis local to the program are instantiated. For the case of an orthogonal 3-d vector space that means the symbolic vectors named ``a0``, ``a1``, and ``a2`` are created. We can instantiate the geometric algebra and obtain the basis vectors with -

.. code:: python

   o3d = Ga('a_1 a_2 a_3', g=[1, 1, 1])
   a_1, a_2, a_3 = o3d.mv()

or use the ``Ga.build()`` function -

.. code:: python

   o3d, a_1, a_2, a_3 = Ga.build('a_1 a_2 a_3', g=[1, 1, 1])

Note that the python variable name for a basis vector does not have to correspond to the name give in ``Ga()`` or ``Ga.build()``, one may wish to use a shortened python variable name to reduce programming (typing) errors, for example one could use -

.. code:: python

   o3d, a1, a2, a3 = Ga.build('a_1 a_2 a_3', g=[1, 1, 1])

or

.. code:: python

   st4d, g0, g1, g2, g3 = Ga.build('gamma_0 gamma_1 gamma_2 gamma_3',
                                   g=[1, -1, -1, -1])

for Minkowski space time.

If the latex printer is used ``e1`` would print as :math:`{\boldsymbol{e_{1}}}` and ``g1`` as :math:`{\boldsymbol{\gamma_{1}}}`.

Representation and Reduction of Multivector Bases
-------------------------------------------------

In our symbolic geometric algebra all multivectors can be obtained from the symbolic basis vectors we have input, via the different operations available to geometric algebra. The first problem we have is representing the general multivector in terms terms of the basis vectors. To do this we form the ordered geometric products of the basis vectors and develop an internal representation of these products in terms of python classes. The ordered geometric products are all multivectors of the form
:math:`a_{i_{1}}a_{i_{2}}\dots a_{i_{r}}` where :math:`i_{1}<i_{2}<\dots <i_{r}` and :math:`r \le n`. We call these multivectors bases and represent them internally with non-commutative symbols so for example :math:`a_{1}a_{2}a_{3}` is represented by

.. code:: python

   Symbol('a_1*a_2*a_3', commutative=False)

In the simplest case of two basis vectors ``a_1`` and ``a_2`` we have a list of bases

.. code:: python

   self.bases = ((Integer(1),)
                 (Symbol('a_1', commutative=False, real=True),
                  Symbol('a_2', commutative=False, real=True)),
                 (Symbol('a_1*a_2', commutative=False, real=True),))

For the case of the basis blades we have

.. code:: python

   self.blades = ((Integer(1),)
                  (Symbol('a_1', commutative=False, real=True),
                   Symbol('a_2', commutative=False, real=True)),
                  (Symbol('a_1^a_2', commutative=False, real=True)))

The index tuples for the bases of each pseudo grade and each grade for the case of dimension 3 is

.. code:: python

   self.indexes = (((),),
                   ((0,), (1,), (2,)),
                   ((0, 1), (0, 2), (1, 2)),
                   ((0, 1, 2),))

Then the non-commutative symbol representing each base is constructed from each index tuple. For example for ``self.indexes[1][1]`` the symbol is ``Symbol('a_1*a_3', commutative=False)``.

Base Representation of Multivectors
-----------------------------------

In terms of the bases defined as non-commutative *sympy* symbols the general multivector is a linear combination (scalar *sympy* coefficients) of bases so that for the case of two bases the most general multivector is given by -

.. code:: python

   A = A_0+A__1*self.bases[1][0]+A__2*self.bases[1][1]+\
       A__12*self.bases[2][0]

If we have another multivector ``B`` to multiply with ``A`` we can calculate the product in terms of a linear combination of bases if we have a multiplication table for the bases.

Blade Representation of Multivectors
------------------------------------

Since we can now calculate the symbolic geometric product of any two multivectors we can also calculate the blades corresponding to the product of the symbolic basis vectors using the formula

.. math:: \be A_{r}{\wedge}b = {\frac{1}{2}}\lp A_{r}b+\lp -1 \rp ^{r}bA_{r} \rp , \ee

where :math:`A_{r}` is a multivector of grade :math:`r` and :math:`b` is a vector. For our example basis the result is shown in Table :ref:`bladexpand`.

.. code-block:: python
   :caption: Blade expansions
   :name: bladexpand

   1 = 1
   a0 = a0
   a1 = a1
   a2 = a2
   a0^a1 = {-(a0.a1)}1+a0a1
   a0^a2 = {-(a0.a2)}1+a0a2
   a1^a2 = {-(a1.a2)}1+a1a2
   a0^a1^a2 = {-(a1.a2)}a0+{(a0.a2)}a1+{-(a0.a1)}a2+a0a1a2

The important thing to notice about Table :ref:`bladexpand` is that it is a triagonal (lower triangular) system of equations so that using a simple back substitution algorithm we can solve for the pseudo bases in terms of the blades giving Table :ref:`baseexpand`.

.. code-block:: python
   :caption: Base expansions
   :name: baseexpand

   1 = 1
   a0 = a0
   a1 = a1
   a2 = a2
   a0a1 = {(a0.a1)}1+a0^a1
   a0a2 = {(a0.a2)}1+a0^a2
   a1a2 = {(a1.a2)}1+a1^a2
   a0a1a2 = {(a1.a2)}a0+{-(a0.a2)}a1+{(a0.a1)}a2+a0^a1^a2

Using Table :ref:`baseexpand` and simple substitution we can convert from a base multivector representation to a blade representation. Likewise, using Table :ref:`bladexpand` we can convert from blades to bases.

Using the blade representation it becomes simple to program functions that will calculate the grade projection, reverse, even, and odd multivector functions.

Note that in the multivector class ``Mv`` there is a class variable for each instantiation, ``self.is_blade_rep``, that is set to ``False`` for a base representation and ``True`` for a blade representation. One needs to keep track of which representation is in use since various multivector operations require conversion from one representation to the other.