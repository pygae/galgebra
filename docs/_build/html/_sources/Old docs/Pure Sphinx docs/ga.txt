.. raw:: html

    <script type="text/javascript" >
        MathJax.Hub.Config({
            TeX: { equationNumbers: { autoNumber: "AMS" } }
        });
    </script>

.. role:: red
   :class: color:red

.. role:: raw-html(raw)
   :format: html
   
*****************
Geometric Algebra
*****************

:Author: Alan Bromborsky

.. |release| replace:: 0.10

.. % Complete documentation on the extended LaTeX markup used for Python
.. % documentation is available in ``Documenting Python'', which is part
.. % of the standard documentation for Python.  It may be found online
.. % at:
.. %
.. % http://www.python.org/doc/current/doc/doc.html
.. % \lstset{language=Python}
.. % \input{macros}
.. % This is a template for short or medium-size Python-related documents,
.. % mostly notably the series of HOWTOs, but it can be used for any
.. % document you like.
.. % The title should be descriptive enough for people to be able to find
.. % the relevant document.

.. % Increment the release number whenever significant changes are made.
.. % The author and/or editor can define 'significant' however they like.

.. % At minimum, give your name and an email address.  You can include a
.. % snail-mail address if you like.

.. % This makes the Abstract go on a separate page in the HTML version;
.. % if a copyright notice is used, it should go immediately after this.
.. %
.. % \ifhtml
.. % \chapter*{Front Matter\label{front}}
.. % \fi
.. % Copyright statement should go here, if needed.
.. % ...
.. % The abstract should be a paragraph or two long, and describe the

.. % scope of the document.

.. topic:: Abstract

   This document describes the implementation, installation and use of a
   geometric algebra module written in
   python that utilizes the sympy symbolic algebra library.  The python
   module ga has been developed for coordinate free calculations using
   the operations (geometric, outer, and inner products etc.) of geometric algebra.
   The operations can be defined using a completely arbitrary metric defined
   by the inner products of a set of arbitrary vectors or the metric can be
   restricted to enforce orthogonality and signature constraints on the set of
   vectors.  Additionally, a metric that is a function of a coordinate set can
   be defined so that a geometric algebra over a manifold can be implemented.
   Geometric algebras over submanifolds of the base manifold are also supported as
   well as linear multivector differential operators and linear transformations.
   In addition the module includes the geometric, outer (curl) and inner
   (div) derivatives. Tensors are included in the module as multilinear
   functions of vectors with contraction and covariant differentiation
   defined without the need of component indices.  For latex output a
   latex distribution must be installed.  A more detail description of the
   module and the mathematics behind it is at `GA[pg 12] <../../../LaTeX_docs/GA.pdf#page=13>`_.

.. module:: sympy.galgebra.ga


What is Geometric Algebra?
==========================

    .. math::

        \newcommand{\bm}[1]{\boldsymbol{#1}}
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
        \newcommand{\braces}[1]{\left \{ {#1} \right \}}
        \newcommand{\grade}[1]{\left < {#1} \right >}
        \newcommand{\f}[2]{{#1}\lp {#2} \rp}
        \newcommand{\paren}[1]{\lp {#1} \rp}
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

    Geometric algebra is the Clifford algebra of a real finite dimensional vector
    space or the algebra that results when the vector space
    is extended with a product of vectors (geometric product) that is associative,
    left and right distributive, and yields a real number for the square (geometric
    product) of any vector [Hestenes]_, [Doran]_.  The elements of the geometric
    algebra are called multivectors and consist of the linear combination of
    scalars, vectors, and the geometric product of two or more vectors. The
    additional axioms for the geometric algebra are that for any vectors :math:`a`,
    :math:`b`, and :math:`c` in the base vector space ([Doran]_,p85):

    .. math::
      :nowrap:

      \begin{equation*}
          \begin{array}{c}
          a\left ( bc \right ) = \left ( ab \right ) c \\
          a\left ( b+c \right ) = ab+ac \\
          \left ( a + b \right ) c = ac+bc \\
          aa = a^{2} \in \Re
          \end{array}
      \end{equation*}

    The dot product of two vectors is defined by ([Doran]_,p86)

    .. math::
      :nowrap:

      \begin{equation}
         a\cdot b \equiv (ab+ba)/2
      \end{equation}

    Then consider

    .. math::
      :nowrap:

      \begin{align}
         c &= a+b \\
         c^{2} &= (a+b)^{2} \\
         c^{2} &= a^{2}+ab+ba+b^{2} \\
         a\cdot b &= (c^{2}-a^{2}-b^{2})/2 \in \Re
      \end{align}

    Thus :math:`a\cdot b`  is real.  The objects generated from linear combinations
    of the geometric products of vectors are called multivectors.  If a basis for
    the underlying vector space is the set of vectors formed from :math:`\boldsymbol{e}_{1},\dots,\boldsymbol{e}_{n}` (we use
    boldface :math:`\boldsymbol{e}`'s to denote basis vectors)
    a complete basis for the geometric algebra is given by the scalar :math:`1`, the vectors :math:`\boldsymbol{e}_{1},\dots,\boldsymbol{e}_{n}`
    and all geometric products of basis vectors

    .. math::
      :nowrap:

      \begin{equation}
          \boldsymbol{e}_{i_{1}}\boldsymbol{e}_{i_{2}}\dots \boldsymbol{e}_{i_{r}} \mbox{ where } 0\le r \le n,\;0 \le i_{j}
                \le n \mbox{ and } i_{1}\lt i_{2}\lt \dots\lt i_{r}
      \end{equation}

    Each base of the complete basis is represented by a noncommutative symbol (except for the scalar 1)
    with name :math:`\boldsymbol{e}_{i_{1}}\dots \boldsymbol{e}_{i_{r}}` so that the general multivector :math:`\boldsymbol{A}` is represented by
    (:math:`A` is the scalar part of the multivector and the :math:`A^{i_{1},\dots,i_{r}}` are scalars)

    .. math::
      :nowrap:

      \begin{equation}
          \boldsymbol{A} = A + \sum_{r=1}^{n}\sum_{i_{1},\dots,i_{r}}^{0\le i_{j}\lt i_{j+1} \le n}
          A^{i_{1},\dots,i_{r}}\boldsymbol{e}_{i_{1}}\boldsymbol{e}_{i_{2}}\dots \boldsymbol{e}_{i_{r}}
      \end{equation}

    The critical operation in setting up the geometric algebra is reducing
    the geometric product of any two bases to a linear combination of bases so that
    we can calculate a multiplication table for the bases.  Since the geometric
    product is associative we can use the operation (by definition for two vectors :math:`a\cdot b \equiv (ab+ba)/2`  which is a scalar)

    .. math::
      :nowrap:

       \begin{equation}\label{reduce}
          \boldsymbol{e}_{i_{j+1}}\boldsymbol{e}_{i_{j}} = 2\boldsymbol{e}_{i_{j+1}}\cdot \boldsymbol{e}_{i_{j}} - \boldsymbol{e}_{i_{j}}\boldsymbol{e}_{i_{j+1}}
       \end{equation}

    These processes are repeated untill every basis list in :math:`\boldsymbol{A}` is in normal
    (ascending) order with no repeated elements. As an example consider the
    following

    .. math::
      :nowrap:

       \begin{align*}
          \boldsymbol{e}_{3}\boldsymbol{e}_{2}\boldsymbol{e}_{1} &= (2(\boldsymbol{e}_{2}\cdot \boldsymbol{e}_{3}) - \boldsymbol{e}_{2}\boldsymbol{e}_{3})\boldsymbol{e}_{1} \\
                          &= 2\left ( \boldsymbol{e}_{2}\cdot \boldsymbol{e}_{3}\right )\boldsymbol{e}_{1} -
                             \boldsymbol{e}_{2}\boldsymbol{e}_{3}\boldsymbol{e}_{1} \\
                          &= 2\left ( \boldsymbol{e}_{2}\cdot \boldsymbol{e}_{3}\right )\boldsymbol{e}_{1} -
                            \boldsymbol{e}_{2}\left ( 2\left ( \boldsymbol{e}_{1}\cdot \boldsymbol{e}_{3}\right ) -
                            \boldsymbol{e}_{1}\boldsymbol{e}_{3}\right ) \\
                          &= 2\left ( \left ( \boldsymbol{e}_{2}\cdot \boldsymbol{e}_{3}\right )\boldsymbol{e}_{1} -
                             \left ( \boldsymbol{e}_{1}\cdot \boldsymbol{e}_{3}\right )\boldsymbol{e}_{2}\right )+\boldsymbol{e}_{2}\boldsymbol{e}_{1}\boldsymbol{e}_{3} \\
                          &= 2\left ( \left ( {\boldsymbol{e}_{2}\cdot \boldsymbol{e}_{3}}\right )\boldsymbol{e}_{1} -
                          \left ( {\boldsymbol{e}_{1}\cdot \boldsymbol{e}_{3}}\right )\boldsymbol{e}_{2}+
                             \left ( \boldsymbol{e}_{1}\cdot \boldsymbol{e}_{2}\right )\boldsymbol{e}_{3}\right )-\boldsymbol{e}_{1}\boldsymbol{e}_{2}\boldsymbol{e}_{3}
       \end{align*}

    which results from repeated application of eq. (:math:`\ref{reduce}`).  If the product of basis vectors contains repeated factors
    eq. (:math:`\ref{reduce}`) can be used to bring the repeated factors next to one another so that if :math:`\boldsymbol{e}_{i_{j}} = \boldsymbol{e}_{i_{j+1}}`
    then :math:`\boldsymbol{e}_{i_{j}}\boldsymbol{e}_{i_{j+1}} = \boldsymbol{e}_{i_{j}}\cdot \boldsymbol{e}_{i_{j+1}}` which is a scalar that commutes with all the terms in the product
    and can be brought to the front of the product.  Since every repeated pair of vectors in a geometric product of :math:`r` factors
    reduces the number of noncommutative factors in the product by :math:`r-2`. The number of bases in the multivector algebra is :math:`2^{n}`
    and the number containing :math:`r` factors is :math:`{n\choose r}` which is the number of combinations or :math:`n` things
    taken :math:`r` at a time (binominal coefficient).

    The other construction required for formulating the geometric algebra is the outer or wedge product (symbol :math:`\wedge`) of :math:`r`
    vectors denoted by :math:`a_{1}\wedge\dots\wedge a_{r}`.  The wedge product of :math:`r` vectors is called an :math:`r`-blade and is defined
    by ([Doran]_,p86)

    .. math::
      :nowrap:

       \begin{equation}
          a_{1}\wedge\dots\wedge a_{r} \equiv \sum_{i_{j_{1}}\dots i_{j_{r}}} \epsilon^{i_{j_{1}}\dots i_{j_{r}}}a_{i_{j_{1}}}\dots a_{i_{j_{1}}}
       \end{equation}

    where :math:`\epsilon^{i_{j_{1}}\dots i_{j_{r}}}` is the contravariant permutation symbol which is :math:`+1` for an even permutation of the
    superscripts, :math:`0` if any superscripts are repeated, and :math:`-1` for an odd permutation of the superscripts. From the definition
    :math:`a_{1}\wedge\dots\wedge a_{r}` is antisymmetric in all its arguments and the following relation for the wedge product of a vector :math:`a` and an
    :math:`r`-blade :math:`\boldsymbol{B}_{r}` can be derived

    .. math::
      :nowrap:

       \begin{equation}\label{wedge}
          a\wedge \boldsymbol{B}_{r} = (a\boldsymbol{B}_{r}+(-1)^{r}\boldsymbol{B}_{r}a)/2
       \end{equation}

    Using eq. (:math:`\ref{wedge}`) one can represent the wedge product of all the basis vectors
    in terms of the geometric product of all the basis vectors so that one can solve (the system
    of equations is lower diagonal) for the geometric product of all the basis vectors in terms of
    the wedge product of all the basis vectors.  Thus a general multivector :math:`\boldsymbol{B}` can be
    represented as a linear combination of a scalar and the basis blades.

    .. math::
      :nowrap:

       \begin{equation}
          \boldsymbol{B} = B + \sum_{r=1}^{n}\sum_{i_{1},\dots,i_{r},\;\forall\; 0\le i_{j} \le n} B^{i_{1},\dots,i_{r}}\boldsymbol{e}_{i_{1}}\wedge \boldsymbol{e}_{i_{2}}\wedge\dots\wedge \boldsymbol{e}_{r}
       \end{equation}

    Using the blades :math:`\boldsymbol{e}_{i_{1}}\wedge \boldsymbol{e}_{i_{2}}\wedge\dots\wedge \boldsymbol{e}_{r}` creates a graded
    algebra where :math:`r` is the grade of the basis blades.  The grade-:math:`r`
    part of :math:`\boldsymbol{B}` is the linear combination of all terms with
    grade :math:`r` basis blades. The scalar part of :math:`\boldsymbol{B}` is defined to
    be grade-:math:`0`.  Now that the blade expansion of :math:`\boldsymbol{B}` is defined
    we can also define the grade projection operator :math:`\left < \boldsymbol{B}\right >_{r}` by

    .. math::
      :nowrap:

       \begin{equation}
          \left < \boldsymbol{B}\right >_{r} = \sum_{i_{1},\dots,i_{r},\;\forall\; 0\le i_{j} \le n} B^{i_{1},\dots,i_{r}}\boldsymbol{e}_{i_{1}}\wedge \boldsymbol{e}_{i_{2}}\wedge\dots\wedge \boldsymbol{e}_{r}
       \end{equation}

    and

    .. math::
      :nowrap:

       \begin{equation}
          \left < \boldsymbol{B}\right >_{} \equiv \left < \boldsymbol{B}\right >_{0} = B
       \end{equation}

    Then if :math:`\boldsymbol{A}_{r}` is an :math:`r`-grade multivector and :math:`\boldsymbol{B}_{s}` is an :math:`s`-grade multivector we have

    .. math::
      :nowrap:

       \begin{equation}
          \boldsymbol{A}_{r}\boldsymbol{B}_{s} = \left < \boldsymbol{A}_{r}\boldsymbol{B}_{s}\right >_{\left | {r-s}\right |}+
                                 \left < \boldsymbol{A}_{r}\boldsymbol{B}_{s}\right >_{\left | {r-s}+2\right |}+\cdots
                                 \left < \boldsymbol{A}_{r}\boldsymbol{B}_{s}\right >_{r+s}
       \end{equation}

    and define ([Hestenes]_,p6)

    .. math::
      :nowrap:

       \begin{align}
          \boldsymbol{A}_{r}\wedge\boldsymbol{B}_{s} &\equiv \left < \boldsymbol{A}_{r}\boldsymbol{B}_{s}\right >_{r+s} \\
          \boldsymbol{A}_{r}\cdot\boldsymbol{B}_{s} &\equiv \left \{ \begin{array}{cc}
          r\mbox{ and }s \ne 0: & \left < \boldsymbol{A}_{r}\boldsymbol{B}_{s}\right >_{\left | {r-s}\right |}  \\
          r\mbox{ or }s = 0: & 0 \end{array} \right \}
       \end{align}

    where :math:`\boldsymbol{A}_{r}\cdot\boldsymbol{B}_{s}` is called the dot or inner product of
    two pure grade multivectors.  For the case of two non-pure grade multivectors

    .. math::
      :nowrap:

       \begin{align}
          \boldsymbol{A}\wedge\boldsymbol{B} &= \sum_{r,s}\left < \boldsymbol{A}\right >_{r}\wedge\left < \boldsymbol{B}\right >_{{s}} \\
          \boldsymbol{A}\cdot\boldsymbol{B} &= \sum_{r,s\ne 0}\left < \boldsymbol{A}\right >_{r}\cdot\left < \boldsymbol{B}\right >_{{s}}
       \end{align}

    Two other products, the right (:math:`\rfloor`) and left (:math:`\lfloor`) contractions, are defined by

    .. math::
      :nowrap:

       \begin{align}
          \boldsymbol{A}\lfloor\boldsymbol{B} &\equiv \sum_{r,s}\left \{\begin{array}{cc} \left < \boldsymbol{A}_r\boldsymbol{B}_{s}\right >_{r-s} & r \ge s \\
                                                      0                                               & r < s \end{array}\right \}  \\
          \boldsymbol{A}\rfloor\boldsymbol{B} &\equiv \sum_{r,s}\left \{\begin{array}{cc} \left < \boldsymbol{A}_{r}\boldsymbol{B}_{s}\right >_{s-r} & s \ge r \\
                                                      0                                               & s < r\end{array}\right \}
       \end{align}

    A final operation for multivectors is the reverse.  If a multivector :math:`\boldsymbol{A}` is the geometric product of :math:`r` vectors (versor)
    so that :math:`\boldsymbol{A} = a_{1}\dots a_{r}` the reverse is defined by

    .. math::
      :nowrap:

       \begin{align}
          \boldsymbol{A}^{\dagger} \equiv a_{r}\dots a_{1}
       \end{align}

    where for a general multivector we have (the the sum of the reverse of versors)

    .. math::
      :nowrap:

       \begin{equation}
          \boldsymbol{A}^{\dagger} = A + \sum_{r=1}^{n}(-1)^{r(r-1)/2}\sum_{i_{1},\dots,i_{r},\;\forall\; 0\le i_{j} \le n} A^{i_{1},\dots,i_{r}}\boldsymbol{e}_{i_{1}}\wedge \boldsymbol{e}_{i_{2}}\wedge\dots\wedge \boldsymbol{e}_{r}
       \end{equation}

    note that if :math:`\boldsymbol{A}` is a versor then :math:`\boldsymbol{A}\boldsymbol{A}^{\dagger}\in\Re` and if
    :math:`\boldsymbol{AA}^{\dagger} \ne 0` then

    .. math::
      :nowrap:

       \begin{equation}
          \boldsymbol{A}^{-1} = \frac{\boldsymbol{A}^{\dagger}}{\boldsymbol{AA}^{\dagger}}.
       \end{equation}


Representation of Multivectors in *sympy*
=========================================

    The sympy python module offers a simple way of representing multivectors using linear
    combinations of commutative expressions (expressions consisting only of commuting sympy objects)
    and noncommutative symbols. We start by defining :math:`n` noncommutative sympy symbols as a basis for
    the vector space

    .. code-block:: python

       (e_1,...,e_n) = symbols('e_1,...,e_n',commutative=False)


    Several software packages for numerical geometric algebra calculations are
    available from Doran-Lasenby group and the Dorst group. Symbolic packages for
    Clifford algebra using orthongonal bases such as
    :math:`\boldsymbol{e}_{i}\boldsymbol{e}_{j}+\boldsymbol{e}_{j}\boldsymbol{e}_{i} = 2\eta_{ij}`, where :math:`\eta_{ij}` is a numeric
    array are available in Maple and Mathematica. The symbolic algebra module,
    *ga*, developed for python does not depend on an orthogonal basis
    representation, but rather is generated from a set of :math:`n` arbitrary
    symbolic vectors :math:`\boldsymbol{e}_{1},\boldsymbol{e}_{2},\dots,\boldsymbol{e}_{n}` and a symbolic metric
    tensor :math:`g_{ij} = \boldsymbol{e}_{i}\cdot \boldsymbol{e}_{j}` (the symbolic metric can be symbolic constants
    or symbolic functions in the case of a manifold).

    All scalar symbolic algebra is handled by the
    python module sympy and the abstract basis vectors are encoded as
    noncommuting sympy symbols.

    The basic geometic algebra operations will be implemented in python by defining
    a geometric algebra class, *Ga*, that performs all required geometric algebra an
    calculus operations on sympy expressions of the form (Einstein summation convention)

    .. math::
      :nowrap:

        \begin{equation}
           F +\sum_{r=1}^{n}F^{i_{1}\dots i_{r}}\boldsymbol{e}_{i_{1}}\dots\boldsymbol{e}_{i_{r}}
        \end{equation}

    where the :math:`F`'s are sympy symbolic constants or functions of the
    coordinates and a multivector class, *Mv*, that wraps *Ga* and overloads the python operators to provide
    all the needed multivector operations as shown in the table of multivector operations
    where *A* and *B*  are any two multivectors (In the case of
    *+*, *-*, *\**, *^*, *|*, *<*, and *>* the operation is also defined if *A* or
    *B* is a sympy symbol or a sympy real number).

    .. image:: tabels/mvops.png
        :width: 350px
        :align: center


    Since *<* and *>* have no r-forms (in python for the *<* and *>* operators there are no *__rlt__()* and
    *__rgt__()* member functions to overload)
    we can only have mixed modes (scalars and multivectors) if the first operand is a multivector.


        Except for *<* and *>* all the multivector operators have r-forms so that as long as one of the
        operands, left or right, is a multivector the other can be a multivector or a scalar (sympy symbol or integer).

        Note that the operator order precedence is determined by python and is not
        necessarily that used by geometric algebra. It is **absolutely essential** to
        use parenthesis in multivector
        expressions containing *\^*, *|*, *<*, and/or *>*.  As an example let
        *A* and *B* be any two multivectors. Then *A + A*B = A +(A\*B)*, but
        *A+A^B = (2\*A)^B* since in python the *\^* operator has a lower precedence
        than the *+* operator.  In geometric algebra the outer and inner products and
        the left and right contractions have a higher precedence than the geometric
        product and the geometric product has a higher precedence than addition and
        subtraction.  In python the *\^*, *|*, *>*, and *<* all have a lower
        precedence than *+* and *-* while *\** has a higher precedence than
        *+* and *-*.

    For those users who wish to define a default operator precedence the functions
    *def\_prec()* and *GAeval()* are available in the module printer.

       *def_prec(gd,op_ord='<>|,^,\*')*

           Define the precedence of the multivector operations.  The function
           *def\_prec()* must be called from the main program and the
           first argument *gd* must be set to *globals()*.  The second argument
           *op_ord* determines the operator precedence for expressions input to
           the function *GAeval()*. The default value of *op_ord* is *<>|,^,\**.
           For the default value the *<*, *>*, and *|* operations have equal
           precedence followed by *^*, and *^* is followed by *\**.


       *GAeval(s,pstr=False)*

           The function *GAeval()* returns a multivector expression defined by the
           string *s* where the operations in the string are parsed according to
           the precedences defined by *define_precedence()*. *pstr* is a flag
           to print the input and output of *GAeval()* for debugging purposes.
           *GAeval()* works by adding parenthesis to the input string *s* with the
           precedence defined by *op_ord='<>|,^,\*'*.  Then the parsed string is
           converted to a sympy expression using the python *eval()* function.
           For example consider where *X*, *Y*, *Z*, and *W* are multivectors

           .. code-block:: python

              def_prec(globals())
              V = GAeval('X|Y^Z*W')

       The sympy variable *V* would evaluate to *((X|Y)^Z)\*W*.


..  _vbm:

Vector Basis and Metric
-----------------------

    The two structures that define the *metric* class (inherited by the
    geometric algebra class) are the
    symbolic basis vectors and the symbolic metric.  The symbolic basis
    vectors are input as a string with the symbol name separated by spaces.  For
    example if we are calculating the geometric algebra of a system with three
    vectors that we wish to denote as *a0*, *a1*, and *a2* we would define the
    string variable:

    .. code-block:: python

       basis = 'a0 a1 a2'

    that would be input into the function which instantiates the geometric
    algebra.  The next step would be
    to define the symbolic metric for the geometric algebra of the basis we
    have defined. The default metric is the most general and is the matrix of
    the following symbols

        .. math::
          :nowrap:

          \begin{equation}\label{metric}
          g = \left [
          \begin{array}{ccc}
            (a0.a0) & (a0.a1)  & (a0.a2) \\
            (a0.a1) & (a1.a1)  & (a1.a2) \\
            (a0.a2) & (a1.a2) & (a2.a2)
          \end{array}
          \right ]
          \end{equation}

    where each of the :math:`g_{ij}` is a symbol representing all of the dot
    products of the basis vectors. Note that the symbols are named so that
    :math:`g_{ij} = g_{ji}` since for the sympy symbols :math:`(a0.a1) \ne (a1.a0)`.
    Note that the strings shown in eq :math:`\ref{metric}` are only used when the values
    of :math:`g_{ij}` are output (printed).   In the *ga* module (library)
    the :math:`g_{ij}` symbols are stored in a member of the geometric algebra
    instance so that if  *o3d* is a geometric algebra then *o3d.g* is
    the metric tensor (:math:`g_{ij} =` *o3d.g[i,j]*) for that algebra.

    The default definition of :math:`g` can be overwritten by specifying a string
    that will define :math:`g`. As an example consider a symbolic representation
    for conformal geometry. Define a basis

       .. code-block:: python

          basis = 'a0 a1 a2 n nbar'

    and a metric

       .. code-block:: python

          g = '# # # 0 0, # # # 0 0, # # # 0 0, 0 0 0 0 2, 0 0 0 2 0'


    then calling *cf3d = Ga(basis,g=g)* would initialize the metric tensor


    .. math::
      :nowrap:

      \begin{equation}
      g = \left [
      \begin{array}{ccccc}
        (a0.a0) & (a0.a1)  & (a0.a2) & 0 & 0\\
        (a0.a1) & (a1.a1)  & (a1.a2) & 0 & 0\\
        (a0.a2) & (a1.a2)  & (a2.a2) & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 \\
        0 & 0 & 0 & 2 & 0
      \end{array}
      \right [
      \end{equation}

    for the  *cf3d* (conformal 3-d) geometric algebra.

    Here we have specified that *n* and *nbar* are orthogonal to all the
    *a*'s, *(n.n) = (nbar.nbar) = 0*, and *(n.nbar) = 2*. Using
    *\#* in the metric definition string just tells the program to use the
    default symbol for that value.

    When *Ga* is called multivector representations of the basis local to
    the program are instantiated.  For the case of an orthogonal 3-d vector
    space that means the
    symbolic vectors named *a0*, *a1*, and *a2* are created. We can
    instantiate the geometric algebra and obtain the basis vectors with -

       .. code-block:: python

          o3d = Ga('a_1 a_2 a_3',g=[1,1,1])
          (a_1,a_2,a_3) = o3d.mv()

    or use the *Ga.build()* function -

       .. code-block:: python

          (o3d,a_1,a_2,a_3) = Ga.build('a_1 a_2 a_3',g=[1,1,1])

    Note that the python variable name for a basis vector does not have to
    correspond to the name give in *Ga()* or *Ga.build()*, one may wish to use a
    shortened python variable name to reduce programming (typing) errors, for
    example one could use -

       .. code-block:: python

          (o3d,a1,a2,a3) = Ga.build('a_1 a_2 a_3',g=[1,1,1])

    or

       .. code-block:: python

          (st4d,g0,g1,g2,g3) = Ga.build('gamma_0 gamma_1 gamma_2 gamma_3',g=[1,-1,-1,-1])

    for Minkowski spacetime.

    If the latex printer is used *e1* would print as :math:`\boldsymbol{e_{1}}`
    and *g1* as :math:`\boldsymbol{\gamma_{1}}`.


          Additionally *Ga()* and *Ga.build()* has simpified options for naming a set of basis vectors and for
          inputing an othogonal basis.


          If one wishes to name the basis vectors :math:`\boldsymbol{e}_{x}`, :math:`\boldsymbol{e}_{y}`, and
          :math:`\boldsymbol{e}_{z}` then set *basis='e*x|y|z'* or to name :math:`\boldsymbol{\gamma}_{t}`,
          :math:`\boldsymbol{\gamma}_{x}`, :math:`\boldsymbol{\gamma}_{y}`, and :math:`\boldsymbol{\gamma}_{z}` then set
          *basis=\'gamma\*t|x|y|z\'*.
          For the case of an othogonal basis if the signature of the
          vector space is :math:`(1,1,1)` (Euclidian 3-space) set *g=[1,1,1]* or if it
          is :math:`(1,-1,-1,-1)` (Minkowsi 4-space) set *g=[1,-1,-1,-1]*. If *g* is a
          function of position then *g* can be entered as a sympy matrix with sympy
          functions as the entries of the matrix or as a list of functions for the
          case of a orthogonal metric.  In the case of spherical coordinates we have
          *g=[1,r**2,r**2*sin(th)**2]*.

Representation and Reduction of Multivector Bases
-------------------------------------------------

    In our symbolic geometric algebra all multivectors
    can be obtained from the symbolic basis vectors we have input, via the
    different operations available to geometric algebra. The first problem we have
    is representing the general multivector in terms terms of the basis vectors.  To
    do this we form the ordered geometric products of the basis vectors and develop
    an internal representation of these products in terms of python classes.  The
    ordered geometric products are all multivectors of the form
    :math:`a_{i_{1}}a_{i_{2}}\dots a_{i_{r}}` where :math:`i_{1}<i_{2}<\dots <i_{r}`
    and :math:`r \le n`. We call these multivectors bases and represent them
    internally with noncommutative symbols so for example :math:`a_{1}a_{2}a_{3}`
    is represented by

        .. code-block: python

            Symbol('a_1*a_2*a_3',commutative=False)

    In the simplest case of two basis vectors *a_1* and *a_2* we have a list of
    bases

        .. code-block: python

            self.bases = [[Symbol('a_1',commutative=False),\
                         Symbol('a_2',commutative=False)],\
                         [Symbol('a_1*a_2',commutative=False)]]

    For the case of the basis blades we have

        .. code-block: python

            self.blades = [[Symbol('a_1',commutative=False),\
                          Symbol('a_2',commutative=False)],\
                          [Symbol('a_1^a_2',commutative=False)]]

    .. note:

      For all grades/pseudo-grades greater than one (vectors) the *\** in the name of the base symbol is
      replaced with a *\^* in the name of the blade symbol so that for all basis bases and
      blades of grade/pseudo-grade greater than one there are different symbols for the corresponding
      bases and blades.

    The index tupels for the bases of each pseudo grade and each grade for the case of dimension 3 is

        .. code-block: python

            self.indexes = (((0,),(1,),(2,)),((0,1),(0,2),(1,2)),((0,1,2)))

    Then the noncommutative symbol representing each base is constructed from each index tuple.
    For example for *self.indexes[1][1]* the symbol is *Symbol('a_1\*a_3',commutative=False)*.

    .. note:

        In the case that the metric tensor is diagonal (orthogonal basis vectors) both base and blade
        bases are identical and fewer arrays and dictionaries need to be constructed.


Base Representation of Multivectors
-----------------------------------

    In terms of the bases defined as noncommutative *sympy* symbols the general multivector
    is a linear combination (scalar *sympy* coefficients) of bases so that for the case
    of two bases the most general multivector is given by -

        *A = A_0+A__1\*self.bases[1][0]+A__2\*self.bases[1][1]+A__12\*self.bases[2][0]*

    If we have another multivector *B* to multiply with *A* we can calculate the product in
    terms of a linear combination of bases if we have a multiplication table for the bases.

Blade Representation of Multivectors
------------------------------------

    Since we can now calculate the symbolic geometric product of any two
    multivectors we can also calculate the blades corresponding to the product of
    the symbolic basis vectors using the formula

        .. math::
          :nowrap:

          \begin{equation}
            A_{r}\wedge b = \frac{1}{2}\left ( A_{r}b-\left ( -1 \right )^{r}bA_{r} \right ),
          \end{equation}

    where :math:`A_{r}` is a multivector of grade :math:`r` and :math:`b` is a
    vector.  For our example basis the result is the table

        .. code-block:: python

           1 = 1
           a0 = a0
           a1 = a1
           a2 = a2
           a0^a1 = {-(a0.a1)}1+a0a1
           a0^a2 = {-(a0.a2)}1+a0a2
           a1^a2 = {-(a1.a2)}1+a1a2
           a0^a1^a2 = {-(a1.a2)}a0+{(a0.a2)}a1+{-(a0.a1)}a2+a0a1a2

    which gives the bases blades in terms of bases.

    The important thing to notice about this expansion is that it is a
    triagonal (lower triangular) system of equations so that using a simple back
    substitution algorithm we can solve for the pseudo bases in terms of the blades
    giving the table

        .. code-block:: python

           1 = 1
           a0 = a0
           a1 = a1
           a2 = a2
           a0a1 = {(a0.a1)}1+a0^a1
           a0a2 = {(a0.a2)}1+a0^a2
           a1a2 = {(a1.a2)}1+a1^a2
           a0a1a2 = {(a1.a2)}a0+{-(a0.a2)}a1+{(a0.a1)}a2+a0^a1^a2

    of the bases in terms of the basis blades.

    Using these tables and simple substitution we can convert from a base
    multivector representation to a blade representation and vice versa.

    Using the blade representation it becomes simple to program functions that will
    calculate the grade projection, reverse, even, and odd multivector functions.

    Note that in the multivector class *Mv* there is a class variable for each
    instantiation, *self.is_blade_rep*, that is set to *False* for a base representation
    and *True* for a blade representation.  One needs to keep track of which
    representation is in use since various multivector operations require conversion
    from one representation to the other.

    .. note:

        When the geometric product of two multivectors is calculated the module looks to
        see if either multivector is in blade representation.  If either is the result of
        the geometric product is converted to a blade representation.  One result of this
        is that if either of the multivectors is a simple vector (which is automatically a
        blade) the result will be in a blade representation.  If *a* and *b* are vectors
        then the result *a\*b* will be *(a.b)+a^b* or simply *a^b* if *(a.b) = 0*

Outer and Inner Products, Left and Right Contractions
=====================================================

In geometric algebra any general multivector :math:`A` can be decomposed into
pure grade multivectors (a linear combination of blades of all the same order)
so that in a :math:`n`-dimensional vector space

.. math::
  :nowrap:

  \begin{equation}
  A = \sum_{r = 0}^{n}A_{r}
  \end{equation}


The geometric product of two pure grade multivectors :math:`A_{r}` and
:math:`B_{s}` has the form

.. math::
  :nowrap:

  \begin{equation}
  A_{r}B_{s} = \left < {A_{r}B_{s}} \right >_{\left |{{r-s}}\right |}+\left < {A_{r}B_{s}} \right >_{\left |{{r-s}}\right |+2}+\cdots+\left < {A_{r}B_{s}} \right >_{r+s}
  \end{equation}


where :math:`\left < { } \right >_{t}` projects the :math:`t` grade components of the
multivector argument.  The inner and outer products of :math:`A_{r}` and
:math:`B_{s}` are then defined to be

.. math::
  :nowrap:

  \begin{equation}
  A_{r}\cdot B_{s} = \left < {A_{r}B_{s}} \right >_{\left |{{r-s}}\right |}
  \end{equation}


.. math::
  :nowrap:

  \begin{equation}
  A_{r}\wedge B_{s} = \left < {A_{r}B_{s}} \right >_{r+s}
  \end{equation}


and

.. math::
  :nowrap:

  \begin{equation}
  A\cdot B = \sum_{r,s > 0}A_{r}\cdot B_{s}
  \end{equation}



.. math::
  :nowrap:

  \begin{equation}
  A\wedge B = \sum_{r,s}A_{r}\wedge B_{s}
  \end{equation}


Likewise the right (:math:`\lfloor`) and left (:math:`\rfloor`) contractions are defined as


.. math::
  :nowrap:

  \begin{equation}
  A_{r}\lfloor B_{s} = \left \{ \begin{array}{cc}
     \left < {A_{r}B_{s}} \right >_{r-s} &  r \ge s \\
               0            &  r < s \end{array} \right \}
  \end{equation}


.. math::
  :nowrap:

  \begin{equation}
  A_{r}\rfloor B_{s} = \left \{ \begin{array}{cc}
     \left < {A_{r}B_{s}} \right >_{s-r} &  s \ge r \\
               0            &  s < r \end{array} \right \}
  \end{equation}


and

.. math::
  :nowrap:

  \begin{equation}
  A\lfloor B = \sum_{r,s}A_{r}\lfloor B_{s}
  \end{equation}

.. math::
  :nowrap:

  \begin{equation}
  A\rfloor B = \sum_{r,s}A_{r}\rfloor B_{s}
  \end{equation}

.. warning::

    In the  *MV* class we have overloaded the *^* operator to represent the outer
    product so that instead of calling the outer product function we can write *mv1^ mv2*.
    Due to the precedence rules for python it is **absolutely essential** to enclose outer products
    in parenthesis.

.. warning::

    In the *MV* class we have overloaded the *|* operator for the inner product,
    *>* operator for the right contraction, and *<* operator for the left contraction.
    Instead of calling the inner product function we can write *mv1|mv2*, *mv1>mv2*, or
    *mv1<mv2* respectively for the inner product, right contraction, or left contraction.
    Again, due to the precedence rules for python it is **absolutely essential** to enclose inner
    products and/or contractions in parenthesis.


.. _reverse:

Reverse of Multivector
======================

    If :math:`A` is the geometric product of :math:`r` vectors

    .. math::
      :nowrap:

      \begin{equation}
        A = a_{1}\dots a_{r}
      \end{equation}


    then the reverse of :math:`A` designated :math:`A^{\dagger}` is defined by

    .. math::
      :nowrap:

      \begin{equation}
        A^{\dagger} \equiv a_{r}\dots a_{1}.
      \end{equation}


    The reverse is simply the product with the order of terms reversed.  The reverse
    of a sum of products is defined as the sum of the reverses so that for a general
    multivector A we have

    .. math::
      :nowrap:

      \begin{equation}
        A^{\dagger} = \sum_{i=0}^{N} {\left < {A} \right >_{i}}^{\dagger}
      \end{equation}


    but

    .. math::
      :nowrap:

      \begin{equation}
        {\left < {A} \right >_{i}}^{\dagger} = \left ( -1\right )^{\frac{i\left ( i-1\right )}{2}}\left < {A} \right >_{i}
      \end{equation}


    which is proved by expanding the blade bases in terms of orthogonal vectors and
    showing that eq. :math:`\ref{eq_4}` holds for the geometric product of orthogonal
    vectors.

    The reverse is important in the theory of rotations in :math:`n`-dimensions.  If
    :math:`R` is the product of an even number of vectors and :math:`RR^{\dagger} = 1`
    then :math:`RaR^{\dagger}` is a composition of rotations of the vector :math:`a`.
    If :math:`R` is the product of two vectors then the plane that :math:`R` defines
    is the plane of the rotation.  That is to say that :math:`RaR^{\dagger}` rotates the
    component of :math:`a` that is projected into the plane defined by :math:`a` and
    :math:`b` where :math:`R=ab`.  :math:`R` may be written
    :math:`R = e^{\frac{\theta}{2}U}`, where :math:`\theta` is the angle of rotation
    and :math:`u` is a unit blade :math:`\left ( u^{2} = \pm 1\right )` that defines the
    plane of rotation.


.. _recframe:

Reciprocal Frames
=================

    If we have :math:`M` linearly independent vectors (a frame),
    :math:`a_{1},\dots,a_{M}`, then the reciprocal frame is
    :math:`a^{1},\dots,a^{M}` where :math:`a_{i}\cdot a^{j} = \delta_{i}^{j}`,
    :math:`\delta_{i}^{j}` is the Kronecker delta (zero if :math:`i \ne j` and one
    if :math:`i = j`). The reciprocal frame is constructed as follows:

    .. math::
      :nowrap:

      \begin{equation}
        E_{M} = a_{1}\wedge\dots\wedge a_{M}
      \end{equation}

    .. math::
      :nowrap:

      \begin{equation}
        E_{M}^{-1} = \frac{E_{M}}{E_{M}^{2}}
      \end{equation}

    Then

    .. math::
      :nowrap:

      \begin{equation}
        a^{i} = \left ( -1\right )^{i-1}\left ( a_{1}\wedge\dots\wedge \breve{a}_{i} \wedge\dots\wedge a_{M}\right ) E_{M}^{-1}
      \end{equation}

    where :math:`\breve{a}_{i}` indicates that :math:`a_{i}` is to be deleted from
    the product.  In the standard notation if a vector is denoted with a subscript
    the reciprocal vector is denoted with a superscript. The set of reciprocal vectors
    will be calculated if a coordinate set is given when a geometric algebra is instantiated since
    they are required for geometric differentiation.

.. _manifold:

Manifolds and Submanifolds
==========================

    A :math:`m`-dimensional vector manifold (By the manifold embedding theorem any :math:`m`-dimensional
    manifold is isomorphic to a :math:`m`-dimensional vector manifold), :math:`\mathcal{M}`, is defined by a
    coordinate tuple (tuples are indicated by the vector accent :math:`\vec{x}`)

    .. math::
      :nowrap:

        \begin{equation}
            \vec{x} = \left ( x^{1},\dots,x^{m} \right ),
        \end{equation}

    and the differentiable mapping (:math:`U^{m}` is an :math:`m`-dimensional subset of :math:`\Re^{m}`)

    .. math::
      :nowrap:

        \begin{equation}
            \boldsymbol{e}^{\mathcal{M}}(\vec{x})\colon U^{m}\subseteq\Re^{m}\rightarrow \mathcal{V},
        \end{equation}

    where :math:`\mathcal{V}` is a vector space with an inner product (:math:`\cdot`) and is of
    :math:`\dim (\mathcal{V}) \ge m`.

    Then a set of basis vectors for the tangent space of :math:`\mathcal{M}` at :math:`\vec{x}`,
    :math:`\mathcal{T}_{\vec{x}}\left( \mathcal{M} \right )`, are

    .. math::
      :nowrap:

        \begin{equation}
            \boldsymbol{e}_{i}^{\mathcal{M}} = \partial_{x^{i}}\boldsymbol{e}^{\mathcal{M}}
        \end{equation}

    and

    .. math::
      :nowrap:

        \begin{equation}
            g_{ij}^{\mathcal{M}}\left (\vec{x}\right ) = \boldsymbol{e}_{i}^{\mathcal{M}}\cdot \boldsymbol{e}_{j}^{\mathcal{M}}.
        \end{equation}

    A :math:`n`-dimensional (:math:`n\le m`) submanifold :math:`\mathcal{N}` of :math:`\mathcal{M}` is defined by
    a coordinate tuple

    .. math::
      :nowrap:

        \begin{equation}
            \vec{u} = \left (u^{1},\dots,u^{n} \right ),
        \end{equation}

    and a differentiable mapping

    .. math::
      :nowrap:

        \begin{equation}\label{eq_79}
            \vec{x}(\vec{u})\colon U^{n}\subseteq\Re^{n}\rightarrow U^{m}\subseteq\Re^{m},
        \end{equation}

    which induces a mapping

    .. math::
      :nowrap:

        \begin{equation}
            \boldsymbol{e}^{\mathcal{M}}\left (\vec{x}(\vec{u})\right )\colon U^{n}\subseteq\Re^{n}\rightarrow \mathcal{V}.
        \end{equation}

    Then the basis vectors for the tangent space :math:`\mathcal{T}_{\vec{u}}\left( \mathcal{N} \right )` are
    (using :math:`\boldsymbol{e}^{\mathcal{N}}(\vec{u}) = \boldsymbol{e}^{\mathcal{M}}(\vec{x}(\vec{u}))` and the chain rule)

    .. math::
      :nowrap:

        \begin{equation}
            \boldsymbol{e}_{i}^{\mathcal{N}}(\vec{u}) = \partial_{u^{i}}\boldsymbol{e}^{\mathcal{N}}(\vec{u})
                                        = \partial_{x^{j}}\boldsymbol{e}^{\mathcal{M}}(\vec{x})\partial_{u^{i}}x^{j}
                                        = \boldsymbol{e}_{j}^{\mathcal{M}}(\vec{x}(\vec{u}))\partial_{u^{i}}x^{j},
        \end{equation}

    and

    .. math::
      :nowrap:

        \begin{equation}\label{eq_53}
            g_{ij}^{\mathcal{N}}(\vec{u}) = \partial_{u^{i}}{x^{k}}\partial_{u^{j}}{x^{l}}
                                                    g_{kl}^{\mathcal{M}}\vec{x}(\vec{u}).
        \end{equation}

    Going back to the base manifold, :math:`\mathcal{M}`, note that the mapping
    :math:`\boldsymbol{e}^{\mathcal{M}}(\vec{x})\colon U^{n}\subseteq\Re^{n}\rightarrow \mathcal{V}` allows us to calculate an unormalized pseudo-scalar
    for :math:`\mathcal{T}_{\vec{x}}(\mathcal{M})`,

    .. math::
      :nowrap:

        \begin{equation}
            I^{\mathcal{M}}(\vec{x}) = \boldsymbol{e}_{1}^{\mathcal{M}}(\vec{x})
                                               \wedge\dots\wedge \boldsymbol{e}_{m}^{\mathcal{M}}(\vec{x}).
        \end{equation}

    With the pseudo-scalar we can define a projection operator from :math:`\mathcal{V}`
    to the tangent space of :math:`\mathcal{M}` by

    .. math::
      :nowrap:

        \begin{equation}
            P_{\vec{x}}(v) = \left (v\cdot I^{\mathcal{M}}(\vec{x})\right)
                                      \left (I^{\mathcal{M}}(\vec{x})\right )^{-1} \;\forall\; v\in\mathcal{V}.
        \end{equation}

    In fact for each tangent space :math:`\mathcal{T}_{\vec{x}}(\mathcal{M})` we can define a geometric algebra
    :math:`\mathcal{G}\left (\mathcal{T}_{\vec{x}}(\mathcal{M})\right )` with pseudo-scalar :math:`I^{\mathcal{M}}` so that if
    :math:`A \in \mathcal{G}(\mathcal{V})` then

    .. math::
      :nowrap:

        \begin{equation}
            P_{\vec{x}}(A) = \left (A\cdot I^{\mathcal{M}}(\vec{x})\right )
                                 \left (I^{\mathcal{M}}(\vec{x})\right )^{-1}
                                 \in \mathcal{G}\left (\mathcal{T}_{\vec{x}}(\mathcal{M})\right )\;\forall\;
                                 A \in \mathcal{G}(\mathcal{V})
        \end{equation}

    and similarly for the submanifold :math:`\mathcal{N}`.

    If the embedding :math:`\boldsymbol{e}^{\mathcal{M}}(\vec{x})\colon U^{n}\subseteq\Re^{n}\rightarrow \mathcal{V}` is not given,
    but the metric tensor :math:`g_{ij}^{\mathcal{M}}(\vec{x})` is given, the geometric algebra of the
    tangent space can be constructed.  Also the derivatives of the basis vectors of the tangent space can
    be calculated from the metric tensor using the Christoffel symbols, :math:`\Gamma_{ij}^{k}(\vec{x})`, where
    the derivatives of the basis vectors are given by

    .. math::
      :nowrap:

        \begin{equation}\label{eq_79a}
            \partial_{x^{i}}\boldsymbol{e}_{j}^{\mathcal{M}} = \Gamma_{ij}^{k}(\vec{x})\boldsymbol{e}_{k}^{\mathcal{M}}.
        \end{equation}

    If we have a submanifold, :math:`\mathcal{N}`, defined by eq. (:math:`\ref{eq_79}`) we can calculate the metric of
    :math:`\mathcal{N}` from eq. (:math:`\ref{eq_53}`) and hence construct the geometric algebra and calculus of the
    tangent space, :math:`\mathcal{T}_{\vec{u}}(\mathcal{N})\subseteq \mathcal{T}_{\vec{x}(\vec{u})}{\mathcal{M}}`.

    If the base manifold is normalized (use the hat symbol to denote normalized tangent vectors,
    :math:`\boldsymbol{\hat{e}}_{i}^{\mathcal{M}}`, and the resulting metric tensor, :math:`\hat{g}_{ij}^{\mathcal{M}}` we have
    :math:`\boldsymbol{\hat{e}}_{i}^{\mathcal{M}}\cdot\boldsymbol{\hat{e}}_{i}^{\mathcal{M}} = \pm 1` and :math:`\hat{g}_{ij}^{\mathcal{M}}` does
    not posess enough
    information to calculate :math:`g_{ij}^{\mathcal{N}}`.  In that case we need to know :math:`g_{ij}^{\mathcal{M}}`, the
    metric tensor of the base manifold before normalization.  Likewise, for the case of a vector
    manifold unless the mapping, :math:`\boldsymbol{e}^{\mathcal{M}}(\vec{x})\colon U^{m}\subseteq\Re^{m}\rightarrow \mathcal{V}`, is
    constant the tangent vectors and metric tensor can only be normalized after the fact (one cannot have a
    mapping that automatically normalizes all the tangent vectors).

.. _deriv:

Geometric Derivative
====================

    The directional derivative of a multivector field :math:`F(x)` is defined by (:math:`a` is a vector and :math:`h` is a scalar)

    .. math::
      :nowrap:

      \begin{equation}\label{eq_50}
         (a\cdot\nabla_{x})F \equiv \lim_{h\rightarrow 0}\frac{F(x+ah)-F(x)}{h}.
      \end{equation}

    Note that :math:`a\cdot\nabla_{x}` is a scalar operator.  It will give a result containing only those grades
    that are already in :math:`F`.  :math:`(a\cdot\nabla_{x})F` is the best linear approximation of :math:`F(x)`
    in the direction :math:`a`.  Equation (:math:`\ref{eq_50}`) also defines the operator :math:`\nabla_{x}` which for a set of
    basis vectors, :math:`\left \{ \boldsymbol{e}_{i}\right \}`, has the representation (note that the :math:`\boldsymbol{e}^{j}` are reciprocal
    basis vectors)

    .. math::
      :nowrap:

      \begin{equation}
          \nabla_{x} F = \boldsymbol{e}^{j}\frac{\partial F}{\partial x^{j}}
      \end{equation}

    If :math:`F_{r}` is a :math:`r`-grade multivector (if the independent vector, :math:`x`, is obvious we suppress it in the
    notation and just write :math:`\nabla`) and
    :math:`F_{r} = F_{r}^{i_{1}\dots i_{r}}\boldsymbol{e}_{i_{1}}\wedge\dots\wedge \boldsymbol{e}_{i_{r}}`
    then

    .. math::
      :nowrap:

      \begin{equation}
        \nabla F_{r} = \frac{\partial F_{r}^{i_{1}\dots i_{r}}}{\partial x^{j}}\boldsymbol{e}^{j}\left ( \boldsymbol{e}_{i_{1}}\wedge
                     \dots\wedge \boldsymbol{e}_{i_{r}} \right )
      \end{equation}

    Note that
    :math:`\boldsymbol{e}^{j}\left (\boldsymbol{e}_{i_{1}}\wedge\dots\wedge \boldsymbol{e}_{i_{r}} \right )`
    can only contain grades :math:`r-1` and :math:`r+1` so that :math:`\nabla F_{r}`
    also can only contain those grades. For a grade-:math:`r` multivector
    :math:`F_{r}` the inner (div) and outer (curl) derivatives are

    .. math::
      :nowrap:

      \begin{equation}
      \nabla\cdot F_{r} = \left < \nabla F_{r}\right >_{r-1} = \boldsymbol{e}^{j}\cdot \partial_{x^{j}}F_{r}
      \end{equation}

    and

    .. math::
      :nowrap:

      \begin{equation}
      \nabla\wedge F_{r} = \left < \nabla F_{r}\right >_{r+1} = \boldsymbol{e}^{j}\wedge \partial_{x^{j}}F_{r}
      \end{equation}

    For a general multivector function :math:`F` the inner and outer derivatives are
    just the sum of the inner and outer dervatives of each grade of the multivector
    function.

Geometric Derivative on a Manifold
----------------------------------

    In the case of a manifold the derivatives of the :math:`\boldsymbol{e}_{i}`'s are functions of the coordinates,
    :math:`\left \{x^{i}\right \}`, so that the geometric derivative of a :math:`r`-grade multivector field is (Einstein summation
    convention)

    .. math::
      :nowrap:

      \begin{align}
            \nabla F_{r} &= \boldsymbol{e}^{i}\partial_{x^{i}}F_{r} = \boldsymbol{e}^{i}\partial_{x^{i}}
                           \left ( F_{r}^{i_{1}\dots i_{r}} \boldsymbol{e}_{i_{1}}\wedge\dots\wedge \boldsymbol{e}_{i_{r}}\right ) \nonumber \\
                         &= \partial_{x^{i}} F_{r}^{i_{1}\dots i_{r}} \boldsymbol{e}^{i}\left ( \boldsymbol{e}_{i_{1}}\wedge\dots\wedge \boldsymbol{e}_{i_{r}} \right )
                            +F_{r}^{i_{1}\dots i_{r}}\boldsymbol{e}^{i}\partial_{x^{i}}\left ( \boldsymbol{e}_{i_{1}}\wedge\dots\wedge \boldsymbol{e}_{i_{r}} \right )
      \end{align}

    where the multivector functions :math:`\boldsymbol{e}^{i}\partial_{x^{i}}\left (\boldsymbol{e}_{i_{1}}\wedge\dots\wedge \boldsymbol{e}_{i_{r}}\right )` are the
    connection for the manifold.  We use the Christoffel symbols of the first kind
    to calculate the derivatives of the basis vectors and the product rule to
    calculate the derivatives of the basis blades where

    .. math::
      :nowrap:

      \begin{equation}
        \Gamma_{ijk} = \frac{1}{2} \left ( \partial_{x^{i}}{g_{jk}}+\partial_{x^{j}}{g_{ik}}-\partial_{x^{k}}{g_{ij}}\right ),
      \end{equation}

    and

    .. math::
      :nowrap:

      \begin{equation}
        \partial_{x^{i}}{ \boldsymbol{e}_{j}} = \Gamma_{ijk} \boldsymbol{e}^{k}.
      \end{equation}


    The Christoffel symbols of the second kind,

    .. math::
      :nowrap:

      \begin{equation}
        \Gamma_{ij}^{k} = \frac{1}{2} g^{kl}\left ( \partial_{x^{j}}{g_{li}}+\partial_{x^{i}}{g_{lj}}-\partial_{x^{l}}{g_{ij}}\right ),
      \end{equation}

    could also be used to calculate the derivatives in term of the original basis vectors, but since we need to calculate the
    reciprocal basis vectors for the geometric derivative
    it is more efficient to use the symbols of the first kind.}

    The directional (material/convective) derivative, :math:`(v\cdot\nabla)F_{r}` is given by

    .. math::
      :nowrap:

      \begin{align}
            \left ( v\cdot\nabla \right ) F_{r} &= v^{i}\partial_{x^{i}}{F_{r}} = v^{i}\partial_{x^{i}}
                           \left ( F_{r}^{i_{1}\dots i_{r}}\boldsymbol{e}_{i_{1}}\wedge\dots\wedge \boldsymbol{e}_{i_{r}}\right ) \nonumber \\
                         &= v^{i}\partial_{x^{i}}{F_{r}^{i_{1}\dots i_{r}}}\left ( \boldsymbol{e}_{i_{1}}\wedge\dots\wedge \boldsymbol{e}_{i_{r}} \right )
                            +v^{i}F_{r}^{i_{1}\dots i_{r}}\partial_{x^{i}}\left ( \boldsymbol{e}_{i_{1}}\wedge\dots\wedge \boldsymbol{e}_{i_{r}}\right ),
      \end{align}

    so that the multivector connection functions for the directional derivative are
    :math:`\partial_{x^{i}}\left ( \boldsymbol{e}_{i_{1}}\wedge\dots\wedge \boldsymbol{e}_{i_{r}}\right )`. Be careful and note that
    :math:`(v\cdot\nabla) F_{r} \ne v\cdot \left (\nabla F_{r}\right )` since the dot and geometric products are
    not associative with respect to one another (:math:`v\cdot\nabla` is a scalar operator).

Normalizing Basis for Derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    The basis vector set, :math:`\left \{ \boldsymbol{e}_{i}\right \}`, is not in general normalized.  We define a normalized set of basis
    vectors, :math:`\left \{ \boldsymbol{\hat{e}}_{i}\right \}`, and reciprocal basis vectors, :math:`\left \{ \boldsymbol{\hat{e}}^{i}\right \}`, by

    .. math::
      :nowrap:

      \begin{align}
          \boldsymbol{\hat{e}}_{i} &= \frac{\boldsymbol{e}_{i}}{\sqrt{\left | \left ( \boldsymbol{e}_{i}\right )^{2} \right |}} = \frac{\boldsymbol{e}_{i}}{\left | \boldsymbol{e}_{i}\right |}, \\
          \boldsymbol{\hat{e}}^{i} &= \frac{\boldsymbol{e}^{i}}{\sqrt{\left | \left (\boldsymbol{e}^{i}\right )^{2} \right |}} = \frac{\boldsymbol{e}^{i}}{\left | \boldsymbol{e}^{i}\right |}.
      \end{align}

    This works for all :math:`\boldsymbol{e}_{i}^{2} \neq 0`.  Note that  :math:`\boldsymbol{\hat{e}}_{i}^{2} = \pm 1` and :math:`\left (\boldsymbol{\hat{e}}^{i}\right )^{2} = \pm 1`.
    Using the definition of reciprocal vectors we obtain the relationship between :math:`\left |\boldsymbol{e}_{i}\right |` and :math:`\left | \boldsymbol{e}^{i}\right |`,

    .. math::
      :nowrap:

      \begin{align}
            \boldsymbol{e}^{i}\cdot \boldsymbol{e}_{j} &= \delta^{i}_{j} \nonumber \\
            \left | \boldsymbol{e}^{i}\right |\boldsymbol{\hat{e}}^{i}\cdot \left | \boldsymbol{e}_{j}\right |\boldsymbol{\hat{e}}_{j} &= \delta^{i}_{j} \nonumber \\
            \left | \boldsymbol{e}^{i}\right |\left | \boldsymbol{e}_{i}\right | &= 1 \nonumber \\
            \left | \boldsymbol{e}^{i}\right | &= \frac{1}{\left | \boldsymbol{e}_{i}\right |}.
      \end{align}

    Thus the geometric derivative for a set of normalized basis vectors is (we assume that
    :math:`F_{r} = F_{r}^{i_{1}\dots i_{r}} \boldsymbol{\hat{e}}_{i_{1}}\wedge\dots\wedge \boldsymbol{\hat{e}}_{i_{r}}`)

    .. math::
      :nowrap:

      \begin{equation}
            \nabla F_{r} =  \boldsymbol{e}^{i}\partial_{x^{i}}{F_{r}} = \frac{\boldsymbol{\hat{e}}^{i}}{\left |  \boldsymbol{e}_{i}\right |}\partial_{x^{i}}{F_{r}}
                         =\partial_{x^{i}}{F_{r}^{i_{1}\dots i_{r}}}\frac{\boldsymbol{\hat{e}}^{i}}{\left | \boldsymbol{e}_{i}\right |}
                           \left ( \boldsymbol{\hat{e}}_{i_{1}}\wedge\dots\wedge\boldsymbol{\hat{e}}_{i_{r}}\right )
                            +F_{r}^{i_{1}\dots i_{r}}\frac{\boldsymbol{\hat{e}}^{i}}{\left | \boldsymbol{e}_{i}\right |}\partial_{x^{i}}
                            \left ( \boldsymbol{\hat{e}}_{i_{1}}\wedge\dots\wedge\boldsymbol{\hat{e}}_{i_{r}}\right ).
      \end{equation}

    Additionally, one can calculate the connection of the normalized basis as follows


    .. math::
      :nowrap:

      \begin{align}
            \partial_{x^{j}}\boldsymbol{e}_{i} =& \partial_{x^{j}}\left ( \left | \boldsymbol{e}_{i}\right |\boldsymbol{\hat{e}}_{i}\right ) = \Gamma_{jik}\boldsymbol{e}^{k}, \nonumber \\
            \partial_{x^{j}}\left | \boldsymbol{e}_{i}\right |\boldsymbol{\hat{e}}_{i}
                                              +\left | \boldsymbol{e}_{i}\right |\partial_{x^{j}}\boldsymbol{\hat{e}}_{i} =& \Gamma_{jik}\boldsymbol{e}^{k}, \nonumber \\
            \partial_{x^{j}}\left | \boldsymbol{e}_{i}\right |\boldsymbol{\hat{e}}_{i}
                                              +\left | \boldsymbol{e}_{i}\right |\partial_{x^{j}}\boldsymbol{\hat{e}}_{i} =& \frac{1}{\left | \boldsymbol{e}_{k}\right |}\Gamma_{jik}\boldsymbol{\hat{e}}^{k}, \nonumber \\
            \partial_{x^{j}}\boldsymbol{\hat{e}}_{i} =& \frac{1}{\left | \boldsymbol{e}_{i}\right |}\left ( \frac{1}{\left | \boldsymbol{e}_{k}\right |}\Gamma_{jik}\boldsymbol{\hat{e}}_{k}
                                               -\partial_{x^{j}}\left | \boldsymbol{e}_{i}\right |\boldsymbol{\hat{e}}_{i}\right ), \nonumber \\
                                            =& \frac{1}{\left | \boldsymbol{e}_{i}\right |\left | \boldsymbol{e}_{k}\right |}\Gamma_{jik}\boldsymbol{\hat{e}}_{k}
                                               -\frac{1}{\left | \boldsymbol{e}_{i}\right |}\partial_{x^{j}}\left | \boldsymbol{e}_{i}\right |\boldsymbol{\hat{e}}_{i},  \nonumber \\
                                            =& \frac{1}{\left | \boldsymbol{e}_{i}\right |\left | \boldsymbol{e}_{k}\right |}\Gamma_{jik}\boldsymbol{\hat{e}}_{k}
                                               -\frac{1}{2g_{ii}}\partial_{x^{j}}g_{ii}\boldsymbol{\hat{e}}_{i}.
      \end{align}


.. _DOPS:

Linear Differential Operators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    First a note on partial derivative notation.  We shall use the following notation for a partial derivative where
    the manifold coordinates are :math:`x_{1},\dots,x_{n}`:

    .. math::
      :nowrap:

        \begin{equation}\label{eq_66a}
            \frac{\partial^{j_{1}+\cdots+j_{n}}}{\partial x_{1}^{j_{1}}\dots\partial x_{n}^{j_{n}}} = \partial_{j_{1}\dots j_{n}}.
        \end{equation}

    If :math:`j_{k}=0` the partial derivative with respect to the :math:`k^{th}` coordinate is not taken.  If the :math:`j_{k} = 0` for all
    :math:`1 \le k \le n` then the partial derivative operator is the scalar one.  If we consider a partial derivative where the :math:`x`'s are
    not in normal order such as

    .. math::
      :nowrap:

        \begin{equation}
            \frac{\partial^{j_{1}+\cdots+j_{n}}}{\partial x_{i_{1}}^{j_{1}}\dots\partial x_{i_{n}}^{j_{n}}},
        \end{equation}

    and the :math:`i_{k}`'s are not in ascending order.  The derivative can always be  put in the form in eq (:math:`\ref{eq_66a}`) since the order
    of differentiation does not change the value of the partial derivative (for the smooth functions we are considering).
    Additionally, using our notation the product of two partial derivative operations is given by

    .. math::
      :nowrap:

        \begin{equation}
            \partial_{i_{1}\dots i_{n}}\partial_{j_{1}\dots j_{n}} = \partial_{i_{1}+j_{1},\dots, i_{n}+j_{n}}.
        \end{equation}

    A general general multivector linear differential operator is a linear combination of multivectors and partial derivative operators
    denoted by (in all of this section we will use the Einstein summation convention)

    .. math::
      :nowrap:

        \begin{equation}\label{eq_66b}
            D \equiv D^{i_{1}\dots i_{n}}\partial_{i_{1}\dots i_{n}}.
        \end{equation}

    Equation (:math:`\ref{eq_66b}`) is the normal form of the differential operator in that the partial derivative operators are written to the right
    of the multivector coefficients and do not operate upon the multivector coefficients.
    The operator of eq (:math:`\ref{eq_66b}`) can operate on mulitvector functions, returning a multivector function via the following definitions.


    :math:`F` as (Einstein summation convention)

    .. math::
      :nowrap:

        \begin{equation}\label{eq_67a}
            D\circ F = D^{j_{1}\dots j_{n}}\circ\partial_{j_{1}\dots j_{n}}F,
        \end{equation}

    or

    .. math::
      :nowrap:

        \begin{equation}\label{eq_68a}
            F\circ D = \partial_{j_{1}\dots j_{n}}F\circ D^{j_{1}\dots j_{n}},
        \end{equation}

    where the :math:`D^{j_{1}\dots j_{n}}` are multivector functions and :math:`\circ` is any of the multivector multiplicative operations.

    Equations (:math:`\ref{eq_67a}`) and (:math:`\ref{eq_68a}`) are not the most general multivector linear differential operators, the most general would be

    .. math::
      :nowrap:

        \begin{equation}
            D(F) = D^{j_{1}\dots j_{n}}\left (\partial_{j_{1}\dots j_{n}}F\right ),
        \end{equation}

    where :math:`D^{j_{1}\dots j_{n}}\left (\right )` are linear multivector functionals.

    The definition of the sum of two differential operators is obvious since any multivector operator, :math:`\circ`, is a bilinear operator
    :math:`\left (\left ( D_{A}+D_{B}\right )\circ F = D_{A}\circ F+D_{B}\circ F\right )`, the product of two differential operators :math:`D_{A}` and :math:`D_{B}`
    operating on a multivector function :math:`F` is defined to be (:math:`\circ_{1}` and :math:`\circ_{2}` are any two multivector multiplicative operations)

    .. math::
      :nowrap:

        \begin{align}
            \left ( D_{A}\circ_{1}D_{B}\right )\circ_{2}F &\equiv \left ( D_{A}^{i_{1}\dots i_{n}}\circ_{1}
                                                          \partial_{i_{1}\dots i_{n}}\left ( D_{B}^{j_{1}\dots j_{n}}
                                                          \partial_{j_{1}\dots j_{n}}\right )\right )\circ_{2}F \nonumber \\
                                                  &= \left (D_{A}^{i_{1}\dots i_{n}}\circ_{1}
                                                     \left ( \left (\partial_{i_{1}\dots i_{n}}D_{B}^{j_{1}\dots j_{n}}\right )
                                                     \partial_{j_{1}\dots j_{n}}+
                                                     D_{B}^{j_{1}\dots j_{n}}\right )
                                                     \partial_{i_{1}+j_{1},\dots, i_{n}+j_{n}}\right )\circ_{2}F \nonumber \\
                                                  &= \left ( D_{A}^{i_{1}\dots i_{n}}\circ_{1}\left (\partial_{i_{1}\dots i_{n}}D_{B}^{j_{1}\dots j_{n}}\right )
                                                     \right )\circ_{2}\partial_{j_{1}\dots j_{n}}F+
                                                     \left (D_{A}^{i_{1}\dots i_{n}}\circ_{1}D_{B}^{j_{1}\dots j_{n}}\right )
                                                     \circ_{2}\partial_{i_{1}+j_{1},\dots, i_{n}+j_{n}}F,
        \end{align}

    where we have used the fact that the :math:`\partial` operator is a scalar operator and commutes with :math:`\circ_{1}` and :math:`\circ_{2}`.

    Thus for a pure operator product :math:`D_{A}\circ D_{B}` we have

    .. math::
      :nowrap:

        \begin{equation}\label{eq_71a}
            D_{A}\circ D_{B} = \left ( D_{A}^{i_{1}\dots i_{n}}\circ\left ( \partial_{i_{1}\dots i_{n}}D_{B}^{j_{1}\dots j_{n}}\right )\right )
                                                     \partial_{j_{1}\dots j_{n}}+
                                                     \left (D_{A}^{i_{1}\dots i_{n}}\circ_{1}D_{B}^{j_{1}\dots j_{n}}\right )
                                                     \partial_{i_{1}+j_{1},\dots, i_{n}+j_{n}}
        \end{equation}

    and the form of eq (:math:`\ref{eq_71a}`) is the same as eq(:math:`\ref{eq_67a}`).  The basis of eq (:math:`\ref{eq_71a}`) is that the :math:`\partial` operator
    operates on all object to the right of it as products so that the product rule must be used in all differentiations.  Since eq (:math:`\ref{eq_71a}`)
    puts the product of two differential operators in standard form we also evaluate :math:`F\circ_{2}\left (D_{A}\circ_{1}D_{B}\right )`.

    We now must distinguish between the following cases.  If :math:`D` is a differential operator and :math:`F` a multivector function should :math:`D\circ F` and
    :math:`F\circ D` return a differential operator or a multivector. In order to be consistent with the standard vector analysis we have :math:`D\circ F`
    return a multivector and :math:`F\circ D` return a differential operator.  The we define the complementary differential operator :math:`\bar{D}` which
    is identical to :math:`D` except that :math:`\bar{D}\circ F` returns a differential operator according to eq (:math:`\ref{eq_71a}`) [#f0]_ and :math:`F\circ\bar{D}` returns a multivector according to eq (:math:`\ref{eq_68a}`).

    A general differential operator is built from repeated applications of the basic operator building blocks :math:`\left )\bar{\nabla}\circ A\right )`,
    :math:`(A\circ\bar{\nabla})`, :math:`(\bar{\nabla}\circ\bar{\nabla})`, and :math:`(A\pm \bar{\nabla})`.  Both :math:`\nabla` and
    :math:`\bar{\nabla}` are represented by the operator

    .. math::
      :nowrap:

        \begin{equation}
            \nabla = \bar{\nabla} = \boldsymbol{e}^{i}\partial_{x^{i}},
        \end{equation}

    but are flagged to produce the appropriate result.

    In the our notation the directional derivative operator is :math:`a\cdot\nabla`, the Laplacian
    :math:`\nabla\cdot\nabla` and the expression for the Riemann tensor, :math:`R^{i}_{jkl}`, is

    .. math::
      :nowrap:

        \begin{equation}
            {(\nabla\wedge\nabla}) \boldsymbol{e}^{i} = \frac{1}{2}R^{i}_{jkl}\left (  \boldsymbol{e}^{j}\wedge \boldsymbol{e}^{k}\right ) \boldsymbol{e}^{l}.
        \end{equation}

    We would use the complement if we wish a quantum mechanical type commutator defining

    .. math::
      :nowrap:

        \begin{equation}
            [x,\nabla] \equiv x\nabla - \bar{\nabla}x,
        \end{equation}

    or if we wish to simulate the dot notation (Doran and Lasenby)

    .. math::
      :nowrap:

        \begin{equation}
            \dot{F}\dot{\nabla} = F\bar{\nabla}.
        \end{equation}

.. _Ltrans:

Linear Transformations
======================

    In the tangent space of a manifold, :math:`\mathcal{M}`, (which is a vector space) a linear transformation is the mapping
    :math:`\underline{T}\colon\mathcal{T}_{\vec{x}}\left (\mathcal{M}\right )\rightarrow\mathcal{T}_{\vec{x}}\left (\mathcal{M}\right )`
    (we use an underline to indicate
    a linear transformation) where for all :math:`x,y\in \mathcal{T}_{\vec{x}}\left (\mathcal{M}\right )` and :math:`\alpha\in\Re` we have

    .. math::
      :nowrap:

        \begin{align}
            \underline{T}(x+y) =& \underline{T}(x) + \underline{T}(y) \\
            \underline{T}(\alpha x) =& \alpha\underline{T}(x)
        \end{align}

    The outermorphism induced by :math:`\underline{T}` is defined for :math:`x_{1},\dots,x_{r}\in\mathcal{T}_{\vec{x}}\left (\mathcal{M}\right )` where
    :math:`r\le\dim\left (\mathcal{T}_{\vec{x}}\left (\mathcal{M}\right )\right )`

    .. math::
      :nowrap:

        \begin{equation}
            \underline{T}\left (x_{1}\wedge\dots\wedge x_{r}\right ) \equiv \underline{T}\left (x_{1}\right )\wedge\dots\wedge\underline{T}\left (x_{r} \right )
        \end{equation}

    If :math:`I` is the pseudo scalar for :math:`\mathcal{T}_{\vec{x}}\left (\mathcal{M}\right )` we also have the following
    definitions for determinate, trace, and adjoint (:math:`\overline{T}`) of :math:`\underline{T}` [#f1]_, [#f2]_

    .. math::
      :nowrap:

        \begin{align}
            \underline{T}(I) \equiv&\; \det\left (\underline{T}\right )I \\
            \mbox{tr}\left (\underline{T}\right ) \equiv&\; \nabla_{y}\cdot\underline{T}(y) \\
            x\cdot \overline{T}(y) \equiv&\; y\cdot \underline{T}(x).
        \end{align}

    If :math:`\left \{ \boldsymbol{e}_{i}\right \}` is a basis for :math:`\mathcal{T}_{\vec{x}}\left (\mathcal{M}\right )` then we can represent :math:`\underline{T}` with the matrix :math:`\underline{T}_{i}^{j}` used
    as follows (Einstein summation convention as usual) -

    .. math::
      :nowrap:

        \begin{equation}\label{eq_97}
            \underline{T}\left( \boldsymbol{e}_{i}\right ) = \underline{T}_{i}^{j}\boldsymbol{e}_{j}.
        \end{equation}

    In eq. (:math:`\ref{eq_97}`) the matrix, :math:`\underline{T}_{i}^{j}`, only has it's usual meaning if the :math:`\left \{ \boldsymbol{e}_{i}\right \}` form an orthonormal Euclidan
    basis (Minkowski spaces not allowed). Equations (:math:`\ref{eq_98}`) through (:math:`\ref{eq_100}`) become

    .. math::
      :nowrap:

        \begin{align}
            \det\left (\underline{T}\right ) =&\; \underline{T}\left (\boldsymbol{e}_{1}\wedge\dots\wedge\boldsymbol{e}_{n}\right )
                                 \left (\boldsymbol{e}_{1}\wedge\dots\wedge\boldsymbol{e}_{n}\right )^{-1}, \label{eq_98}\\
            \mbox{tr}\left (\underline{T} \right ) =&\; \underline{T}_{i}^{i},\\
            \overline{T}_{j}^{i} =&\;  g^{il}g_{jp}\underline{T}_{l}^{p}.\label{eq_100}
        \end{align}

.. _MLtrans:

Multi-Linear Transformations (Tensors)
======================================

A multivector multilinear function is a
multivector function :math:`T \paren{ A_{1},\dots,A_{r}}` that is linear in each of it arguments (it could be implicitly
non-linearly dependent on a set of additional arguments such as the postion coordinates, but we only consider the linear arguments).
:math:`T` is a *tensor* of degree :math:`r` if each variable :math:`A_{j}` is restricted to the vector space :math:`\mathcal{V}_{n}`.
More generally if each :math:`A_{j}\in \mathcal{G}\left (\mathcal{V}_{n}\right )` (the geometric algebra of :math:`\mathcal{V}_{n}`),
we call :math:`T` an *extensor* of degree-:math:`r` on :math:`\mathcal{G}\left ( \mathcal{V}_{n}\right )`.

If the values of :math:`T\left (a_{1},\dots,a_{r}\right )` :math:`\left ( a_{j}\in\mathcal{V}_{n}\;\forall\; 1\le j \le r \right )`
are :math:`s`-vectors (pure grade :math:`s` multivectors in :math:`\mathcal{G}\paren{\mathcal{V}_{n}}` we say that
:math:`T` has grade :math:`s` and rank :math:`r+s`.  A tensor of grade zero is called a *multilinear form*.

In the normal definition of tensors as multilinear functions the tensor is defined as a mapping

    .. math::
      :nowrap:

        \begin{equation}

            T:\bigotimes_{i=1}^{r}\mathcal{V}_{i}\rightarrow\Re,

        \end{equation}

so that the standard tensor definition is an example of a grade zero
degree/rank :math:`r` tensor in our definition.

Algebraic Operations
--------------------

The properties of tensors are (:math:`\alpha\in\Re`, :math:`a_{j},b\in\mathcal{V}_{n}`, :math:`T` and :math:`S` are tensors of
rank :math:`r`, and :math:`\circ` is any multivector multiplicative operation)

    .. math::
      :nowrap:

        \begin{align}
            T\left (a_{1},\dots,\alpha a_{j},\dots,a_{r}\right ) =& \alpha T\left (a_{1},\dots,a_{j},\dots,a_{r}\right ), \\
            T\left (a_{1},\dots,a_{j}+b,\dots,a_{r}\right ) =& T\left (a_{1},\dots,a_{j},\dots,a_{r}\right ) +
                                                             T\left (a_{1},\dots,a_{j-1},b,a_{j+1},\dots,a_{r}\right ), \\
            \left ( T\pm S\right )\left (a_{1},\dots,a_{r}\right ) \equiv& T\left (a_{1},\dots,a_{r} \right )
                         \pm S\left (a_{1},\dots,a_{r}\right).
        \end{align}


Now let :math:`T` be of rank :math:`r` and :math:`S` of rank :math:`s` then the product of the two tensors is

    .. math::
      :nowrap:

        \begin{equation}

            \left ( T\circ S\right )\left (a_{1},\dots,a_{r+s}\right ) \equiv T\left (a_{1},\dots,a_{r}\right )
                       \circ S\left (a_{r+1},\dots,a_{r+s}\right ),

        \end{equation}

where ":math:`\circ`" is any multivector multiplicative operation.

Covariant, Contravariant, and Mixed Representations
---------------------------------------------------

The arguments (vectors) of the multilinear fuction can be represented in terms of the basis vectors or the reciprocal basis vectors

    .. math::
      :nowrap:

        \begin{align}
            a_{j} =& a^{i_{j}}\eb_{i_{j}}, \\
                  =& a_{i_{j}}\eb^{i_{j}}.
        \end{align}

These equations gives :math:`a_{j}` in terms of the basis vectors or the reciprocal basis vectors. The index
:math:`j` refers to the argument slot and the indices :math:`i_{j}` the components of the vector in terms of the basis.  The Einstein summation
convention is used throughout.  The covariant representation of the tensor is defined by

    .. math::
      :nowrap:

        \begin{align}
            T_{i_{1}\dots i_{r}} \equiv& \f{T}{\eb_{i_{1}},\dots,\eb_{i_{r}}} \\
            \f{T}{a_{1},\dots,a_{r}} =& \f{T}{a^{i_{1}}\eb_{i_{1}},\dots,a^{i_{r}}\eb_{i_{r}}} \nonumber \\
                                     =& \f{T}{\eb_{i_{1}},\dots,\eb_{i_{r}}}a^{i_{1}}\dots a^{i_{r}} \nonumber \\
                                     =& T_{i_{1}\dots i_{r}}a^{i_{1}}\dots a^{i_{r}}.
        \end{align}

Likewise for the contravariant representation

    .. math::
      :nowrap:

        \begin{align}
            T^{i_{1}\dots i_{r}} \equiv& \f{T}{\eb^{i_{1}},\dots,\eb^{i_{r}}} \\
            \f{T}{a_{1},\dots,a_{r}} =& \f{T}{a_{i_{1}}\eb^{i_{1}},\dots,a_{i_{r}}\eb^{i_{r}}} \nonumber \\
                                     =& \f{T}{\eb^{i_{1}},\dots,\eb^{i_{r}}}a_{i_{1}}\dots a_{i_{r}} \nonumber \\
                                     =& T^{i_{1}\dots i_{r}}a_{i_{1}}\dots a_{i_{r}}.
        \end{align}

One could also have a mixed representation

    .. math::
      :nowrap:

        \begin{align}
            T_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{r}} \equiv& \f{T}{\eb_{i_{1}},\dots,\eb_{i_{s}},\eb^{i_{s+1}}\dots\eb^{i_{r}}} \\
            \f{T}{a_{1},\dots,a_{r}} =& \f{T}{a^{i_{1}}\eb_{i_{1}},\dots,a^{i_{s}}\eb_{i_{s}},
                                        a_{i_{s+1}}\eb^{i_{s}}\dots,a_{i_{r}}\eb^{i_{r}}} \nonumber \\
                                     =& \f{T}{\eb_{i_{1}},\dots,\eb_{i_{s}},\eb^{i_{s+1}},\dots,\eb^{i_{r}}}
                                        a^{i_{1}}\dots a^{i_{s}},a_{i_{s+1}},\dots a_{i_{r}} \nonumber \\
                                     =& T_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{r}}a^{i_{1}}\dots a^{i_{s}},a_{i_{s+1}},\dots a_{i_{r}}.
        \end{align}

In the representation of :math:`T` one could have any combination of covariant (lower) and contravariant (upper) indices.

To convert a covariant index to a contravariant index simply consider

    .. math::
      :nowrap:

        \begin{align}
            \f{T}{\eb_{i_{1}},\dots,\eb^{i_{j}},\dots,\eb_{i_{r}}} =& \f{T}{\eb_{i_{1}},\dots,g^{i_{j}k_{j}}\eb_{k_{j}},\dots,\eb_{i_{r}}} \nonumber \\
                                                                   =& g^{i_{j}k_{j}}\f{T}{\eb_{i_{1}},\dots,\eb_{k_{j}},\dots,\eb_{i_{r}}} \nonumber \\
            T_{i_{1}\dots i_{j-1}i_{j+1}\dots i_{r}}^{i_{j}} =& g^{i_{j}k_{j}}T_{i_{1}\dots i_{j}\dots i_{r}}.
        \end{align}

Similarly one could lower an upper index with :math:`g_{i_{j}k_{j}}`.

Contraction and Differentiation
-------------------------------

The contraction of a tensor between the :math:`j^{th}` and :math:`k^{th}` variables (slots) is

    .. math::
      :nowrap:

        \be
            \f{T}{a_{i},\dots,a_{j-1},\nabla_{a_{k}},a_{j+1},\dots,a_{r}} = \nabla_{a_{j}}\cdot\lp\nabla_{a_{k}}\f{T}{a_{1},\dots,a_{r}}\rp.
        \ee

This operation reduces the rank of the tensor by two.  This definition gives the standard results for *metric contraction* which is
proved as follows for a rank :math:`r` grade zero tensor (the circumflex ":math:`\breve{\:\:}`" indicates that a term is to be deleted from the product).

    .. math::
      :nowrap:

        \begin{align}
            \f{T}{a_{1},\dots,a_{r}} =& a^{i_{1}}\dots a^{i_{r}}T_{i_{1}\dots i_{r}} \\
            \nabla_{a_{j}}T =& \eb^{l_{j}} a^{i_{1}}\dots\lp\partial_{a^{l_j}}a^{i_{j}}\rp\dots a_{i_{r}}T_{i_{1}\dots i_{r}} \nonumber \\
            =& \eb^{l_{j}}\delta_{l_{j}}^{i_{j}} a^{i_{1}}\dots \breve{a}^{i_{j}}\dots a^{i_{r}}T_{i_{1}\dots i_{r}} \\
            \nabla_{a_{m}}\cdot\lp\nabla_{a_{j}}T\rp =& \eb^{k_{m}}\cdot\eb^{l_{j}}\delta_{l_{j}}^{i_{j}}
                                                      a^{i_{1}}\dots \breve{a}^{i_{j}}\dots\lp\partial_{a^{k_m}}a^{i_{m}}\rp
                                                      \dots a^{i_{r}}T_{i_{1}\dots i_{r}} \nonumber \\
                                                     =& g^{k_{m}l_{j}}\delta_{l_{j}}^{i_{j}}\delta_{k_{m}}^{i_{m}}
                                                      a^{i_{1}}\dots \breve{a}^{i_{j}}\dots\breve{a}^{i_{m}}
                                                      \dots a^{i_{r}}T_{i_{1}\dots i_{r}} \nonumber \\
                                                     =& g^{i_{m}i_{j}}a^{i_{1}}\dots \breve{a}^{i_{j}}\dots\breve{a}^{i_{m}}
                                                      \dots a^{i_{r}}T_{i_{1}\dots i_{j}\dots i_{m}\dots i_{r}} \nonumber \\
                                                     =& g^{i_{j}i_{m}}a^{i_{1}}\dots \breve{a}^{i_{j}}\dots\breve{a}^{i_{m}}
                                                      \dots a^{i_{r}}T_{i_{1}\dots i_{j}\dots i_{m}\dots i_{r}}  \nonumber \\
                                                     =& \lp g^{i_{j}i_{m}}T_{i_{1}\dots i_{j}\dots i_{m}\dots i_{r}}\rp a^{i_{1}}\dots
                                                      \breve{a}^{i_{j}}\dots\breve{a}^{i_{m}}\dots a^{i_{r}}\label{eq108}
        \end{align}

Equation :math:`\eqref{eq108}` is the correct formula for the metric contraction of a tensor.

Finally if :math:`\f{T}{a_{1},\dots,a_{r}}` is a tensor field (implicitly a function of position) the tensor derivative is defined as

    .. math::
      :nowrap:

        \begin{align}
            \f{T}{a_{1},\dots,a_{r};a_{r+1}} \equiv \lp a_{r+1}\cdot\nabla\rp\f{T}{a_{1},\dots,a_{r}},
        \end{align}

assuming the :math:`a^{i_{j}}` coefficients are not a function of the coordinates.

This gives for a grade zero rank :math:`r` tensor

    .. math::
      :nowrap:

        \begin{align}
            \lp a_{r+1}\cdot\nabla\rp\f{T}{a_{1},\dots,a_{r}} =& a^{i_{r+1}}\partial_{x^{i_{r+1}}}a^{i_{1}}\dots a^{i_{r}}
                                                                T_{i_{1}\dots i_{r}}, \nonumber \\
                                                             =& a^{i_{1}}\dots a^{i_{r}}a^{i_{r+1}}
                                                                \partial_{x^{i_{r+1}}}T_{i_{1}\dots i_{r}}.
        \end{align}

Covariant Deriviatives
----------------------

The component free form of the covariant derivative (the one used to calculate it in the code) is

    .. math::
      :nowrap:

        \begin{equation}\label{cderiv}
            \mathcal{D}_{a_{r+1}} \f{T}{a_{1},\dots,a_{r};x} \equiv \nabla T
                - \sum_{k=1}^{r}\f{T}{a_{1},\dots,\paren{a_{r+1}\cdot\nabla} a_{k},\dots,a_{r};x}.
        \end{equation}

The effect of :math:`\paren{a_{r+1}\cdot\nabla} a_{k}` in equation :math:`\eqref{cderiv}` is to
parallel transport  :math:`a_{k}` in the direction of :math:`a_{r+1}` which gives the standard definition of the covariant derivative
of a tensor field.  Note that  :math:`a_{k} = a_{k}^{j}\bm{e}_{k}` where :math:`a_{k}^{j}` is not
a function of the coordinates, but in general :math:`\bm{e}_{k}` is a function of the coordinates.


Numpy, LaTeX, and Ansicon Installation
======================================

To install the geometric algebra module on windows,linux, or OSX perform the following operations

    #. Install sympy.  *galgebra* is included in sympy.

    #. To install texlive in linux or windows

        #. Go to <http://www.tug.org/texlive/acquire-netinstall.html> and click on "install-tl.zip" to download
        #. Unzip "install-tl.zip" anywhere on your machine
        #. Open the file "readme.en.html" in the "readme-html.dir" directory.  This file contains the information needed to install texlive.
        #. Open a terminal (console) in the "install-tl-XXXXXX" directory
        #. Follow the instructions in "readme.en.html" file to run the install-tl.bat file in windows or the install-tl script file in linux.

    #. For OSX install mactex from <http://tug.org/mactex/>.

    #. Install python-nympy if you want to calculate numerical matrix functons (determinant, inverse, eigenvalues, etc.).
       For windows go to <http://sourceforge.net/projects/numpy/files/NumPy/1.6.2/> and install the distribution of numpy
       appropriate for your system.  For OSX go to <http://sourceforge.net/projects/numpy/files/NumPy/1.6.1/>.
    #. It is strongly suggested that you go to <http://www.geany.org/Download/Releases> and install the version of the "geany" editor appropriate for your system.
    #. If you wish to use "enhance_print" on windows -

        #. Go to <https://github.com/adoxa/ansicon/downloads> and download "ansicon"
        #. In the Edit -> Preferences -> Tools menu of "geany" enter into the Terminal input the full path of "ansicon.exe"

In addition to the code shown in the examples section of this document there are more examples in the Examples directory under the
*galgebra* directory.

Module Components
=================

Instantiating a Geometric Algebra
---------------------------------

A geometric algebra is instantiated with

   *sympy.galgebra.GA(basis, g=None, coords=None, norm=False, debug=False, X=None)*

   The *basis* and *g* parameters were described in section :ref:`vbm`.

   If *debug=True* the data structure required to initialize the *Ga* class
   are printed out.

   *coords* is a tuple of *sympy* symbols equal in length to
   the number of basis vectors.  These symbols are used as the arguments of a
   multivector field as a function of position and for calculating the derivatives
   of a multivector field. Additionally, *Ga()* calculates the pseudo scalar,
   :math:`I` and makes them available to the programmer as *MV.I* and *MV.Iinv*.
   For the case of instantiating a 3-d geometeric algebra in spherical coordinates we could use

       .. code-block:: python

          (r, th, phi) = coords = symbols('r,theta,phi', real=True)
          basis = 'e_r e_theta e_phi'
          g = [1, r**2, r**2*sin(th)**2]
          sp3d = Ga(basis,g=g,coords=coords,norm=True)

   The input :math:`X` allows the metric to be input as a vector manifold. :math:`X`
   is a list of functions of *coords* dimension, :math:`m`, equal to or greater than
   the number of coordinates. If *g=None* it is assumed that *X* is a vector in an
   :math:`m`-dimensional orthonormal Euclidian vector space.  If it is wished
   the embedding vector space to be non-Euclidian that condition is specified with
   *g*.  For example if we wish the embedding space to be a 5-dimensional Minkowski
   space then *g=[-1,1,1,1,1]*.  Then the *Ga* class uses *X* to calculate the
   manifold basis vectors as a function of the coordinates and from them the metric
   tensor.

   If *norm=True* the basis vectors of the manifold are normalized so that the
   absolute values of the squares of the basis vectors are one. `It is suggested
   that one only use this option for diagonal metric tensors, and even there due so
   with caution, due to the possible
   problems with taking the square root of a general` *sympy* `expression (one that has an
   unknown sign).`

   In addition to the basis vectors, if coordinates are defined for the geometric algebra, the
   left and right geometric derivative operators are calculated and accessed with the *Ga*
   member function *grads()*.

   *Ga.grads()*

       *Ga.grads()* returns a tuple with the left and right geometric derivative operators. A
       typical usage would be

    .. code-block:: python

        (grad,rgrad) = sp3d.grads()

   for the spherical 3-d geometric algebra. The left derivative *grad* :math:`= \nabla` and the
   right derivative *rgrad* :math:`= \bar{\nabla}` have been explained in section :ref:`DOPS`. Again
   the names *grad* and *rgrad* are whatever the user chooses them to be.

   an alternative instantiation method is

   *Ga.build(basis, g=None, coords=None, X=None, norm=False, debug=False)*

     The input parameters for *Ga.build()* are the same as for *Ga()*.  The difference is
     that in addition to returning the geometric algebra *Ga.build()* returns the basis vectors
     at the same time. Using *Ga.build()* in the previous example gives

     .. code-block:: python

       (r, th, phi) = coords = symbols('r,theta,phi', real=True)
       basis = 'e_r e_theta e_phi'
       g = [1, r**2, r**2*sin(th)**2]
       (sp3d,er,eth,ephi) = Ga.build(basis,g=g,coord=coords,norm=True)

   To access the pseudo scalar of the geometric algebra us the member function *I()*.

   *Ga.I()*

       *Ga.I()* returns the normalized pseudo scalar (:math:`\left | {I^{2}}\right |=1`) for the
       geometric algebra. For example :math:`I` = *o3d.I()* for the *o3d* geometric
       algebra.

   In general we have defined member fuctions of the *Ga* class that will instantiate objects
   of other classes since the objects of the other classes are all associated with a particular
   geometric algebra object.  Thus we have

    .. image:: tabels/class_objs.png
        :width: 400px
        :align: center

   for the instantiation of various objects from the *Ga* class.  This means that in order to
   instantiate any of these objects we need only to import *Ga* into our program.


Instantiating a Multivector
---------------------------

    Since we need to associate each multivector with the geometric algebra that contains it
    we use a member function of *Ga* to instantiate every multivector (There is a
    multivector class, *Mv*, but in order the insure that every multivector is associated
    with the correct geometric algebra we always use the member function *Ga.mv* to instantiate
    the multivector.)  The multivector is instantiated with:

    *Ga.mv(name, mode, f=False)*

        As an example of both instantiating a geometric algebra and multivectors consider the
        following code fragment for a 3-d Euclidian geometric algebra.

        .. code-block:: python

            from sympy import symbols
            from ga import Ga
            (x, y, z) = coords = symbols('x,y,z',real=True)
            o3d = Ga('e_x e_y e_z', g=[1,1,1], coords=coords)
            (ex, ey, ez) = o3d.mv()
            V = o3d.mv('V','vector',f=True)

        First consider the multivector instantiation *V = o3d.mv('V','vector',f=True)*.  Here
        a 3-dimensional multivector field that is a function of *x*, *y*, and *z* (*f=True*) is
        being instantiated.  If latex output were used (to be discussed later) the multivector
        *V* would be displayed as

        .. math::
          :nowrap:

          \begin{equation}
            A^{x}\boldsymbol{e}_{x} + A^{y}\boldsymbol{e}_{y} + A^{z}\boldsymbol{e}_{z}
          \end{equation}

        Where the coefficients of the basis vectors are generalized *sympy* functions of the
        coordinates.  The superscripts (Denoted in text output by *A__x*, etc. so
        that for text output *A* would be printed as *A__x\*e_x+A__y\*e_y+A__z\*e_z*) are formed
        from the coordinate symbols or if there are no coordinates from the subscripts of
        the basis vectors.  The types of name and modes available for multivector instantiation are

            .. image:: tabels/instanciate_mv.png
                :width: 750px
                :align: center

        Line 5 of the previous listing illustrates the case of using the *mv* member function with
        no arguments. The code does not return a multivector, but rather a tuple or the basis vectors of the geometric algebra *o3d*.
        The elements of the tuple then can
        be used to construct multivectors, or multivector fields through the operations
        of addition, subtraction, multiplication (geometric, inner, and outer products and left and right contraction).
        As an example we could construct the vector function

        .. code-block:: python

            F = x**2*ex + z*ey + x*y*ez

        or the bivector function

        .. code-block:: python

            B = z*(ex^ey) + y*(ey^ez) + y*(ex^ez).

    If one wished to calculate the left and right geometric derivatives of *F* and *B* the required code would be

    .. code-block:: python

        (grad,rgrad) = o3d.grads()
        dF = grad*F
        dB = grad*B
        dFr = F*rgrad
        dBr = B*rgrad

    *dF*, *dB*, *dFr*, and *dBr* are all multivector functions. For the code where the order of the operations are
    reversed

    .. code-block:: python

        (grad,rgrad) = o3d.grads()
        dFop = F*grad
        dBop = B*grad
        dFrop = rgrad*F
        dBrop = rgrad*B

    *dFop*, *dBop*, *dFrop*, and *dBrop* are all multivector differential operators (again see section :ref:`DOPS`).


Basic Multivector Class Functions
---------------------------------

    *convert_to_blades(self)*

       Convert multivector from the base representation to the blade representation.
       If multivector is already in blade representation nothing is done.


    *convert_from_blades(self)*

       Convert multivector from the blade representation to the base representation.
       If multivector is already in base representation nothing is done.


    *diff(self,var)*

       Calculate derivative of each multivector coefficient with resepect to
       variable *var* and form new multivector from coefficients.


    *dual(self)*

       Return dual of multivector which is multivector left multiplied by
       pseudoscalar *Mv.i* (Hestenes,p22).

    *even(self)*

       Return the even grade components of the multivector.


    *exp(self,hint='+')*

        Return exponential of a multivector :math:`A` if :math:`A^{2}` is a scalar (if :math:`A^{2}` is not a scalar an
        error message is generated).  If :math:`A` is the multivector then :math:`\boldsymbol{e}^{A}` is returned
        where the default *hint*, *+*, assumes :math:`A^{2} > 0` so that


        .. math::
          :nowrap:

            \begin{equation}
                    \boldsymbol{e}^{A} = \cosh\sqrt{A^{2}}+\sinh\sqrt{A^{2}}\left (\frac{A}{\sqrt{A^{2}}}\right )
            \end{equation}

        If the mode is not *+* then :math:`A^{2} < 0` is assumed so that


        .. math::
          :nowrap:

            \begin{equation}
                    \boldsymbol{e}^{A} = \cos \sqrt{-A^{2}}+\sin\sqrt{-A^{2}}\left (\frac{A}{\sqrt{-A^{2}}}\right ).
            \end{equation}

        The hint is required for symbolic multivectors :math:`A` since in general *sympy* cannot determine if
        :math:`A^{2}` is positive or negative.  If :math:`A` is purely numeric the hint is ignored.


    *expand(self)*

       Return multivector in which each coefficient has been expanded using
       sympy *expand()* function.


    *factor(self)*

       Apply the sympy *factor* function to each coefficient of the multivector.


    *Fmt(self, fmt=1,title=None)*

        Function to print multivectors in different formats where

            .. image :: tabels/fmt_opts.png
                :width: 400px
                :align: center

        *title* appends a title string to the beginning of the output.  An equal sign in
        the title string is not required, but is added as a default.


    *func(self,fct)*

       Apply the *sympy* scalar function *fct* to each coefficient of the multivector.


    *grade(self,igrade=0)*

        Return a multivector that consists of the part of the multivector of
        grade equal to *igrade*.  If the multivector has no *igrade* part
        return a zero multivector.


    *inv(self)*

       Return the inverse of the multivector :math:`M` (*M.inv()*) if :math:`MM^{\dagger}` is a nonzero
       scalar.  If :math:`MM^{\dagger}`
       is not a scalar the program exits with an error message.


    *norm(self)*

       Return the norm of the multivector :math:`M` (*M.norm()*) defined by :math:`\sqrt{MM^{\dagger}}` if
       :math:`MM^{\dagger}` is a
       scalar (a sympy scalar
       is returned).  If :math:`MM^{\dagger}` is not a scalar the program exits with an error message.


    *norm2(self)*

       Return the square of the norm of the multivector :math:`M` (*M.norm2()*) defined by :math:`MM^{\dagger}`
       if :math:`MM^{\dagger}`
       is a scalar (a sympy scalar
       is returned).  If :math:`MM^{\dagger}` is not a scalar the program exits with an error message.


    *proj(self,bases_lst)*

       Return the projection of the multivector :math:`M` (*M.proj(bases_lst)*) onto the subspace defined by the list of bases
       (*bases_lst*).


    *scalar(self)*

        Return the coefficient (sympy scalar) of the scalar part of a
        multivector.


    *simplify(self,mode=simplify)*

       *mode* is a sympy simplification function of a list/tuple of sympy
       simplification functions that are applied in sequence (if more than
       one function) each coefficient of the multivector.  For example if
       we wished to applied *trigsimp* and *ratsimp* sympy functions to the
       mulitvector *F* the code would be


       .. code-block:: python

          Fsimp = F.simplify(mode=[trigsimp,ratsimp]).


       Actually *simplify* could be used to apply any scalar sympy function to
       the coefficients of the multivector.


    *subs(self,x)*

       Return multivector where sympy *subs* function has been applied to each
       coefficient of multivector for argument dictionary/list *x*.


    *rev(self)*

       Return the reverse of the multivector.


    *set_coef(self,grade,base,value)*

       Set the multivector coefficient of index *(grade,base)* to *value*.


    *trigsimp(self,\*\*kwargs)*

       Apply the sympy trignometric simplification function *trigsimp* to
       each coefficient of the multivector. *\*\*kwargs* are the arguments of
       *trigsimp*.  See sympy documentation on *trigsimp* for more information.

Basic Multivector Functions
---------------------------

    *Com(A,B)*

       Calulate commutator of multivectors *A* and *B*.  Returns *(AB-BA)/2*.

    *GAeval(s,pstr=False)*

       Returns multivector expression for string *s* with operator precedence for
       string *s* defined by inputs to function *def_prec()*.  if *pstr=True*
       *s* and *s* with parenthesis added to enforce operator precedence are printed.

    *Nga(x,prec=5)*

       If *x* is a multivector with coefficients that contain floating point numbers, *Nga()*
       rounds all these numbers to a precision of *prec* and returns the rounded multivector.

    *ReciprocalFrame(basis,mode='norm')*

       If *basis* is a list/tuple of vectors, *ReciprocalFrame()* returns a tuple of reciprocal
       vectors.  If *mode=norm* the vectors are normalized.  If *mode* is anything other than
       *norm* the vectors are unnormalized and the normalization coefficient is added to the
       end of the tuple.  One must divide by the coefficient to normalize the vectors.

    *ScalarFunction(TheFunction)*

       If *TheFuction* is a real *sympy* fuction a scalar multivector function is returned.

    *cross(v1,v2)*

       If *v1* and *v2* are 3-dimensional euclidian vectors the vector cross product is
       returned, :math:`v_{1}\times v_{2} = -I\left ( v_{1}\wedge v_{2} \right )`.

    *def_prec(gd,op_ord='<>|,^,\*')*

       This is used with the *GAeval()* function to evaluate a string representing a multivector
       expression with a revised operator precedence.  *def_prec()* redefines the operator
       precedence for multivectors. *def_prec()* must be called in the main program an the
       argument *gd* must be *globals()*.  The argument *op_ord* defines the order of operator
       precedence from high to low with groups of equal precedence separated by commas. the default
       precedence *op_ord='<>|,^,\*'* is that used by Hestenes.

    *dual(M)*

       Return the dual of the multivector *M*, math:`MI^{-1}`.

    *inv(B)*

       If for the multivector :math:`B`, :math:`BB^{\dagger}` is a nonzero scalar, return
       :math:`B^{-1} = B^{\dagger}/(BB^{\dagger})`.

    *proj(B,A)*

       Project blade *A* on blade *B* returning :math:`\left ( A\lfloor B\right ) B^{-1}`.

    *refl(B,A)*

       Reflect blade *A* in blade *B*. If *r* is grade of *A* and *s* is grade of *B*
       returns :math:`(-1)^{s(r+1)}BAB^{-1}`.

    *rot(itheta,A)*

       Rotate blade *A* by 2-blade *itheta*.  Is is assumed that *itheta\*itheta > 0* so that
       the rotation is Euclidian and not hyperbolic so that the angle of
       rotation is *theta = itheta.norm()*.  Ther in 3-dimensional Euclidian space. *theta* is the angle of rotation (scalar in radians) and
       *n* is the vector axis of rotation.  Returned is the rotor *cos(theta)+sin(theta)*N* where *N* is
       the normalized dual of *n*.



Multivector Derivatives
-----------------------

    The various derivatives of a multivector function is accomplished by
    multiplying the gradient operator vector with the function.  The gradiant
    operation vector is returned by the *Ga.mv()* function if coordinates
    are defined.  For example if we have for a 3-D vector space

    .. code-block:: python

        X = (x,y,z) = symbols('x y z')
        o3d = Ga('e*x|y|z',metric='[1,1,1]',coords=X)
        (ex,ey,ez) = o3d.mv()
        (grad,rgrad) = o3d.grads()


    Then the gradient operator vector is *grad* (actually the user can give
    it any name he wants to).  Then the derivatives of the multivector
    function *F = o3d.mv('F','mv',f=True)* are given by multiplying by the
    left geometric derivative operator and the right geometric derivative operator
    *grad* :math:`= \nabla` and *rgrad* :math:`= \bar{\nabla}`.  Another option
    is to use the gradiant operator members of the geometric algebra directly where we have
    :math:`\nabla =` *o3d.grad* and :math:`\bar{\nabla} =` *o3d.rgrad*.

    .. math::
      :nowrap:

      \begin{align*}
            \nabla F &=  \mbox{grad $*$ F} \\
            F \bar{\nabla} &=  \mbox{F $*$ rgrad} \\
            \nabla \wedge F &=  \mbox{grad ^ F} \\
            F \wedge \bar{\nabla} &=  \mbox{F ^ rgrad} \\
            \nabla \cdot F &=  \mbox{grad $|$ F} \\
            F \cdot \bar{\nabla} &=  \mbox{F $|$ rgrad} \\
            \nabla \lfloor F &=  \mbox{grad $<$ F} \\
            F \lfloor \bar{\nabla} &=  \mbox{F $<$ rgrad} \\
            \nabla \rfloor F &=  \mbox{grad $>$ F} \\
            F \rfloor \bar{\nabla} &= \mbox{F $>$ rgrad}
      \end{align*}


    The preceding code block gives examples of all possible multivector
    derivatives of the multivector function *F* where the operation returns
    a multivector function. The complementary operations

        .. math::
          :nowrap:

          \begin{align*}
                F \nabla &=  \mbox{F $*$ grad} \\
                \bar{\nabla} F &=  \mbox{rgrad $*$ F} \\
                F \wedge \nabla &=  \mbox{F ^ grad} \\
                \bar{\nabla} \wedge F &=  \mbox{rgrad ^ F} \\
                F \cdot \nabla &=  \mbox{F $|$ grad} \\
                \bar{\nabla}\cdot F &=  \mbox{rgrad $|$ F} \\
                F \lfloor \nabla &=  \mbox{F $<$ grad} \\
                \bar{\nabla} \lfloor F &=  \mbox{rgrad $<$ F} \\
                F \rfloor \nabla &=  \mbox{F $>$ grad} \\
                \bar{\nabla} \rfloor F &= \mbox{rgrad $>$ F}
          \end{align*}

    all return multivector linear differential operators.


Submanifolds
------------

    In general the geometric algebra that the user defines exists on the tangent space of
    a manifold.  The submanifold class, *Sm*, is derived from
    the *Ga* class and allows one
    to define a submanifold of a manifold by defining a coordinate mapping between the submanifold
    coordinates and the manifold coordinates.  What is returned as the submanifold is the geometric
    algebra of the tangent space of the submanifold. The submanifold for a geometric algebra is
    instantiated with

    *Ga.sm(map,coords,root='e',norm=False)*

        To define the submanifold we must define a coordinate map from the coordinates of the submanifold to
        each of the coordinates of the base manifold.  Thus the arguments *map* and *coords* are
        respectively lists of functions and symbols.  The list of symbols, *coords*, are the coordinates of the
        submanifold and are of length equal to the dimension of the submanifold.  The list of functions, *map*,
        define the mapping from the coordinate space of the submanifold to the coordinate space of the
        base manifold.  The length of *map* is equal to the dimension of the base manifold and each function in
        *map* is a function of the coordinates of the submanifold. As a concrete example consider the
        following code.

        .. literalinclude:: latex/manifold_src.py
            :language: python

        The program output is

        .. figure:: latex/manifold_src.png
           :width: 350px
           :align: center

        The base manifold, *sp3d*, is a 3-d Euclidian space using standard spherical coordinates. The submanifold
        *sph2d* of *sp3d* is a spherical surface of radius :math:`1`.  To take the sumanifold operation one step further
        the submanifold *cir1d* of *sph2d* is a circle in *sph2d* where the latitude of the circle is :math:`\pi/8`.

        In each case, for demonstration purposes, a scalar and vector function on each manifold is defined (*f* and *F*
        for the 2-d manifold and *h* and *H* for the 1-d manifold) and the geometric derivative of each function is taken.  The
        manifold mapping and the metric tensor for *cir1d* of *sph2d* are also shown. Note that if the submanifold
        basis vectors are not normalized the program output is.

        .. figure:: latex/man_unnorm.png
           :width: 500px
           :align: center

Linear Transformations
----------------------

    The mathematical background for linear transformations is in section :ref:`Ltrans`.  Linear transformations on the tangent space of
    the manifold are instantiated with the *Ga* member function *lt* (the actual class being instantiated is *Lt*) as shown in
    lines 12, 20, 26, and 44 of the code (Ltrans.py)

    .. literalinclude:: pysource/Ltrans.py
       :language: python


    In all of the examples in the above code the default instantiation is used which produces a general (all the
    coefficients of the linear transformation are symbolic constants) linear transformation. Note that to
    instantiate linear transformations
    coordinates, :math:`\left \{\boldsymbol{e}_{i} \right \}`, must be defined when the geometric algebra associated with the
    linear transformation is instantiated.
    This is due to the naming conventions of the general linear transformation (coordinate names are used)
    and for the calculation
    of the trace of the linear transformation which requires taking a divergence.
    To instantiate a specific linear transformation
    the usage of *lt()* is  *Ga.lt(M,f=False)*.

    *M* is an expression that can define the coefficients of the linear transformation in various ways
    defined as follows.

    .. image:: tabels/lt_ops.png
        :width: 800px
        :align: center

    *f* is *True* or *False*. If *True* the symbolic coefficients of the general linear transformation are
    instantiated as functions of the coordinates.

    The different methods of instantiation are demonstrated in the code (*LtransInst.py*)

    .. literalinclude:: latex/LtransInst.py
        :language: python

    with output

    .. figure:: latex/LtransInst.png
       :width: 500px
       :align: center

    The member functions of the *Lt* class are

    *Lt(A)*

        Returns the image of the multivector :math:`A` under the linear transformation :math:`L` where
        :math:`L(A)` is defined by the
        linearity of :math:`L`, the vector values :math:`L\left ( e_{i}\right )`, and the definition
        :math:`L \left ( e_{i_{1}}\wedge\dots\wedge e_{i_{r}}\right ) = L\left ( e_{i_{1}} \right )
        \wedge\dots\wedge L\left ( e_{i_{r}}\right )`.

    *Lt.det()*

        Returns the determinant (a scalar) of the linear transformation, :math:`L`, defined by
        :math:`\det (L)I = L(I)`.

    *Lt.adj()*

        Returns the adjoint (a linear transformation) of the linear transformation, :math:`L`, defined by
        :math:`a\cdot L(b) = b\cdot \bar{L}(a)` where
        :math:`a` and :math:`b` are any two vectors in the tangent space and :math:`\bar{L}` is the adjoint of :math:`L`.

    *Lt.tr()*

        Returns the trace (a scalar) of the linear transformation, :math:`L`, defined by
        :math:`\f{\operatorname{tr}}{L}=\nabla_{a}\cdot\f{L}{a}` where :math:`a` is a vector in the tangent space.

    *Lt.matrix()*

        Returns the matrix representation (sympy Matrix) of the linear transformation, :math:`L`, defined by
        :math:`L\left ( e_{i}\right ) = L_{ij} e_{j}` where :math:`L_{ij}` is the matrix representation.

    The *Ltrans.py* demonstrate the use of the various *Lt* member functions and operators. The operators that can
    be used with
    linear transformations are *+*, *-*, and *\**. If :math:`A` and :math:`B` are linear transformations,
    :math:`V` a multivector,
    and :math:`\alpha` a
    scalar then :math:`(A\pm B)(V) = A(V)\pm B(V)`, :math:`(AB)(V) = A(B(V))`, and
    :math:`(\alpha A)(V) = \alpha A(V)`.

    The ouput of *Ltrans.py* is

    .. figure:: latex/Ltrans.png
       :width: 750px
       :align: center

Differential Operators
----------------------

    For the mathematical treatment of linear mulivector differential operators see section :ref:`DOPS`.  The is a differential
    operator class *Dop*. However, one never needs to use it directly.  The operators are constructed from linear
    combinations of multivector products of the operators *Ga.grad* and *Ga.rgrad* as shown in the following code for
    both orthogonal rectangular and spherical 3-d coordinate systems.

    .. literalinclude:: latex/dop.py
        :language: python

    The output of this code is.

    .. figure:: latex/dop.png
       :width: 600px
       :align: center

Instantiating a Multi-linear Functions (Tensors)
------------------------------------------------

The mathematical background for multi-linear functions is in section :ref:`MLtrans`.  To instantiate a multi-linear function use

*Mlt(self, f, Ga, nargs=None, fct=False)* where the arguments are

        .. image:: tabels/tensor_inst.png
            :width: 600px
            :align: center

Examples of different types of instantiation are

    .. code-block:: python

        #TensorDef.py

        import sys
        from sympy import symbols,sin,cos
        from printer import Format,xpdf,Get_Program,Print_Function
        from ga import Ga
        from lt import Mlt

        coords = symbols('t x y z',real=True)
        (st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=coords)

        A = st4d.mv('T','bivector')

        def TA(a1,a2):
            global A
            return A | (a1 ^ a2)

        T = Mlt(TA,st4d) # Define multi-linear function

Basic Multilinear Function Class Functions
------------------------------------------

If we can instantiate multilinear functions we can use all the multilinear function class functions as described as follows.
See section :ref:`MLtrans` for the mathematical description of each operation.

    *self(kargs)*

       Calling function to evaluates multilinear function for \T{kargs} list of vector arguments and returns a value.  Note that a
       sympy scalar is returned, \emph{not} a multilinear function.


    *self.contract(slot1,slot2)*

        Returns contraction of tensor between *slot1* and *slot2* where *slot1* is the index of the first vector argument and
        *slot2* is the index of the second vector argument of the tensor. For example if we have a rank two tensor, *T(a1,a2)*,
        then *T.contract(1,2)* is the contraction of *T*.  For this case since there are only two slots there can only be one
        contraction.


    *self.pdiff(slot)*

        Returns gradient of tensor, *T*, with respect to slot vector.  For example if the tensor is :math:`\f{T}{a_{1},a_{2}}` then
        *T.pdiff(2)* is :math:`\nabla_{a_{2}}T`.  Since *T* is a scalar function, *T.pdiff(2)* is a vector function.


    *self.cderiv()*

        Returns covariant derivative of tensor field. If *T* is a tensor of rank :math:`k` then *T.cderiv()* is a tensor of rank :math:`k+1`.
        The operation performed is defined in section :ref:`MLtrans`.

Standard Printing
-----------------

    Printing of multivectors is handled by the module *printer* which contains
    a string printer class derived from the sympy string printer class and a latex
    printer class derived from the sympy latex printer class.  Additionally, there
    is an *Eprint* class that enhances the console output of sympy to make
    the printed output multivectors, functions, and derivatives more readable.
    *Eprint* requires an ansi console such as is supplied in linux/osx or the
    program *ansicon* (github.com/adoxa/ansicon) for windows which replaces *cmd.exe*.

    For a windows user the simplest way to implement ansicon is to use the *geany*
    editor and in the Edit :math:`\rightarrow` Preferences :math:`\rightarrow` Tools menu replace *cmd.exe* with
    *ansicon.exe* (be sure to supply the path to *ansicon*).

    If *Eprint* is called in a program (linux) when multivectors are printed
    the basis blades or bases are printed in bold text, functions are printed in red,
    and derivative operators in green.

    For formatting the multivector output there is the member function

    *Fmt(self,fmt=1,title=None)*

        *Fmt* is used to control how the multivector is printed with the argument
        *fmt*.  If *fmt=1* the entire multivector is printed on one line.  If
        *fmt=2* each grade of the multivector is printed on one line.  If *fmt=3*
        each component (base) of the multivector is printed on one line.  If a
        *title* is given then *title=multivector* is printed.  If the usual print
        command is used the entire multivector is printed on one line.

Latex Printing
--------------

    For latex printing one uses one functions from the *ga* module and one
    function from the *printer* module.  The
    functions are

    *Format(Fmode=True,Dmode=True,ipy=False)*

       This function from the *ga* module turns on latex printing with the
       following options

        .. image:: tabels/latex_prnt.png
            :width: 700px
            :align: center

    *xpdf(filename=None,debug=False,paper=(14,11),crop=False)*

       This function from the *printer* module post-processes the output captured from
       print statements, writes the resulting latex strings to the file *filename},
       processes the file with pdflatex, and displays the resulting pdf file.   All latex files except
       the pdf file are deleted. If *debug = True* the file *filename* is printed to
       standard output for debugging purposes and *filename* (the tex file) is saved.  If *filename* is not entered the default
       filename is the root name of the python program being executed with *.tex* appended.
       The *paper* option defines the size of the paper sheet for latex. The format for the *paper* is

        .. image:: tabels/latex_paper.png
            :width: 600px
            :align: center

       The default of *paper=(14,11)* was chosen so that long multivector expressions would not be truncated on the display.

       If the *crop* input is *True* the linux *pdfcrop* program is used to crop the pdf output (if output is one page).  This only works
       for linux installations (where *pdfcrop* is installed).

       The *xpdf* function requires that latex and a pdf viewer be installed on
       the computer.

       *xpdf* is not required when printing latex in IPython notebook.


    As an example of using the latex printing options when the following code is
    executed

    .. literalinclude:: latex/print_example1.py
        :language: python

    The following is displayed

    .. figure:: latex/print_example1.png
       :width: 700px
       :align: center

    For the cases of derivatives the code is

    .. literalinclude:: latex/print_example2.py
        :language: python

    and the latex displayed output is (:math:`f` is a scalar function)

    .. figure:: latex/print_example2.png
       :width: 700px
       :align: center

    This example also demonstrates several other features of the latex printer.  In the
    case that strings are input into the latex printer such as *grad\*\\bm{A}*,
    *grad^\\bm{A}*, or *grad\*\\bm{A}'!.  The text symbols *grad*, *^*, *|*, and
    *\** are mapped by the *xpdf()* post-processor as follows if the string contains
    an *=*.

    .. image:: tabels/latex_replc.png
        :width: 375px
        :align: center

    If the first character in the string to be printed is a *\%* none of the above substitutions
    are made before the latex processor is applied.  In general for the latex
    printer strings are assumed to be in a math environment (equation or
    align) unless the first character in the string is a *\#*. [#f3]_

      Except where noted the conventions for latex printing follow those of the
      latex printing module of sympy. This includes translating sympy variables
      with Greek name (such as *alpha*) to the equivalent Greek symbol
      (math:`\alpha`) for the purpose of latex printing.  Also a single
      underscore in the variable name (such as *X_j*) indicates a subscript
      (:math:`X_{j}`), and a double underscore (such as *X__k*) a
      superscript (:math:`X^{k}`).  The only other change with regard to the
      sympy latex printer is that matrices are printed full size (equation
      displaystyle).

For formatting the tensor (*Mlt*) output there is the member function

*Fmt(self,cnt=1,title=None)*

    *Fmt* is used to control how the tensor is printed with the argument
    *cnt*.  If *cnt=1* the each tensor component is printed on one line.  If
    *fmt=n* :math:`n` tensor components are printed on one line.  If a
    *title* is given then *title=tensor* is printed.  If the usual print
    command is used one tensor component is printed on one line. If *cnt* is
    greater or equal to the number of tensor components then the entire tensor
    is printer on one line.

Examples
========


Algebra
-------

BAC-CAB Formulas
^^^^^^^^^^^^^^^^

This example demonstrates the most general metric tensor

.. math::
  :nowrap:

  \begin{equation}
  g_{ij} = \left [ \begin{array}{cccc} \left ( a\cdot a\right )  & \left ( a\cdot b\right )  & \left ( a\cdot c\right )  & \left ( a\cdot d\right )  \\
  \left ( a\cdot b\right )  & \left ( b\cdot b\right )  & \left ( b\cdot c\right )  & \left ( b\cdot d\right )  \\
  \left ( a\cdot c\right )  & \left ( b\cdot c\right )  & \left ( c\cdot c\right )  & \left ( c\cdot d\right )  \\
  \left ( a\cdot d\right )  & \left ( b\cdot d\right )  & \left ( c\cdot d\right )  & \left ( d\cdot d\right )
  \end{array}\right ]
  \end{equation}

and how the *galgebra* module can be used to verify and expand geometric algebra identities consisting of relations between
the abstract vectors :math:`a`, :math:`b`, :math:`c`, and :math:`d`.

    .. literalinclude:: latex/baccab.py
        :language: python

The preceeding code block also demonstrates the mapping of *\ **, *^*, and *|* to appropriate latex
symbols.

.. note::

  The :math:`\times` symbol is the commutator product of two multivectors, :math:`A\times B = (AB-BA)/2`.

The output of the code is

.. image:: latex/baccab.png
    :width: 500px
    :align: center

Reciprocal Frame
^^^^^^^^^^^^^^^^

The reciprocal frame of vectors with respect to the basis vectors is required
for the evaluation of the geometric dervative.  The following example demonstrates
that for the case of an arbitrary 3-dimensional Euclidian basis the reciprocal
basis vectors are correctly calculated.

    .. literalinclude:: latex/recp_frame.py
        :language: python

The preceeding code also demonstrated the use of the *\%* directive in
printing a string so that *^* is treated literally and not translated
to *\\wedge*. Note that "%E^{2} =" is printed as ":math:`E^{2} =`"
and not as ":math:`E\wedge {2} =`".

.. image:: latex/recp_frame.png
    :width: 800px
    :align: center

The formulas derived for :math:`E1`, :math:`E2`, :math:`E3`, and :math:`E^{2}` could
also be applied to the numerical calculations of crystal properties.

Lorentz-Transformation
^^^^^^^^^^^^^^^^^^^^^^

A simple physics demonstation of geometric algebra is the derivation of
the Lorentz-Transformation.  In this demonstration a 2-dimensional
Minkowski space is defined and the Lorentz-Transformation is generated
from a rotation of a vector in the Minkowski space using the rotor
:math:`R`.

    .. literalinclude:: latex/lorentz_trans.py
        :language: python

The preceeding code also demonstrates how to use the sympy *subs* functions
to perform the hyperbolic half angle transformation.  The code also shows
the use of both the *#* and *\%* directives in the text string
``r"#%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} = t'\bm{\gamma'_{t}}+x'\bm{\gamma'_{x}} = R\left ( t'\bm{\gamma_{t}}+x'\bm{\gamma_{x}}\right ) R^{\dagger}"``.
Both the *#* and *\%* are needed in this text string for two reasons.  First, the text string contains an *=* sign.  The latex preprocessor
uses this a key to combine the text string with a sympy expression to be printed after the text string.  The *#* is required to inform
the preprocessor that there is no sympy expression to follow.  Second, the *\%* is requires to inform the preprocessor that the text
string is to be displayed in latex math mode and not in text mode (if *#* is present the default latex mode is text mode unless
overridden by the *\%* directive).

.. image:: latex/lorentz_trans.png
    :width: 500px
    :align: center

Linear Transformations
^^^^^^^^^^^^^^^^^^^^^^

Examples of linear transformations produced by the following code -

    .. literalinclude:: latex/Ltrans.py
        :language: python

Note the trace of the general basis 2-d case and the adjoint for the Minkowski 4-d case.  Also note
the results of the adjoint test for both cases.

.. image:: latex/Ltrans.png
    :width: 800px
    :align: center

Calculus
--------

Derivatives in Spherical Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following code shows how to use *galgebra* to use spherical coordinates.
The gradient of a scalar function, :math:`f`, the divergence and curl
of a vector function, :math:`A`, and the exterior derivative (curl) of
a bivector function, :math:`B` are calculated.  Note that to get the
standard curl of a 3-dimension function the result is multiplied by
:math:`-I` the negative of the pseudoscalar.

.. note::

    In geometric calculus the operator :math:`\nabla^{2}` is well defined
    on its own as the geometic derivative of the geometric derivative.
    However, if needed we have for the vector function :math:`A` the relations
    (since :math:`\nabla\cdot A` is a scalar it's curl is equal to it's
    geometric derivative and it's divergence is zero) -

    .. math::
        :nowrap:

        \begin{align*}
        \nabla A =& \nabla\wedge A + \nabla\cdot A \\
        \nabla^{2} A =& \nabla\left ( {{\nabla\wedge A}} \right ) + \nabla\left ( {{\nabla\cdot A}} \right ) \\
        \nabla^{2} A =& \nabla\wedge\left ( {{\nabla\wedge A}} \right ) + \nabla\cdot\left ( {{\nabla\wedge A}} \right )
        +\nabla\wedge\left ( {{\nabla\cdot A}} \right ) + \nabla\cdot\left ( {{\nabla\cdot A}} \right ) \\
        \nabla^{2} A =& \nabla\wedge\left ( {{\nabla\wedge A}} \right ) + \left ( {{\nabla\cdot\nabla}} \right ) A
        - \nabla\left ( {{\nabla\cdot A}} \right ) + \nabla\left ( {{\nabla\cdot A}} \right ) \\
        \nabla^{2} A =& \nabla\wedge\nabla\wedge A + \left ( {{\nabla\cdot\nabla}} \right )A
        \end{align*}



    In the derivation we have used that :math:`\nabla\cdot\left ( {{\nabla\wedge A}} \right ) = \left ( {{\nabla\cdot\nabla}} \right )A - \nabla\left ( {{\nabla\cdot A}} \right )`
    which is implicit in the second *BAC-CAB* formula.
    No parenthesis is needed for the geometric curl of the curl (exterior derivative of exterior derivative)
    since the :math:`\wedge` operation is associative unlike the vector curl operator and :math:`\nabla\cdot\nabla` is the usual Laplacian
    operator.

.. literalinclude:: latex/spherical.py
    :language: python

Results of code

.. image:: latex/spherical.png
    :width: 800px
    :align: center

Maxwell's Equations
^^^^^^^^^^^^^^^^^^^

The geometric algebra formulation of Maxwell's equations is deomonstrated
with the formalism developed in "Geometric Algebra for Physicists" [Doran]_.
In this formalism the signature of the metric is :math:`(1,-1,-1,-1)` and the
basis vectors are :math:`\gamma_{t}`, :math:`\gamma_{x}`, :math:`\gamma_{y}`,
and :math:`\gamma_{z}`.  The if :math:`\boldsymbol{E}` and :math:`\boldsymbol{B}` are the
normal electric and magnetic field vectors the electric and magnetic
bivectors are given by :math:`E = \boldsymbol{E}\gamma_{t}` and :math:`B = \boldsymbol{B}\gamma_{t}`.
The electromagnetic bivector is then :math:`F = E+IB` where
:math:`I = \gamma_{t}\gamma_{x}\gamma_{y}\gamma_{z}` is the pesudo-scalar
for the Minkowski space.  Note that the electromagnetic bivector is isomorphic
to the electromagnetic tensor.  Then if :math:`J` is the 4-current all of
Maxwell's equations are given by :math:`\boldsymbol{\nabla}F = J`.  For more details
see [Doran]_ chapter 7.

.. literalinclude:: latex/maxwell.py
    :language: python

Code output

.. image:: latex/maxwell.png
    :width: 950px
    :align: center

Dirac Equation
^^^^^^^^^^^^^^

In [Doran]_ equation 8.89 (page 283) is the geometric algebra formulation of the Dirac equation.  In this equation
:math:`\psi` is an 8-component real spinor which is to say that it is a multivector with sacalar, bivector, and
pseudo-vector components in the space-time geometric algebra (it consists only of even grade components).

.. literalinclude:: latex/dirac.py
    :language: python

The equations displayed are the partial differential equations for each component of the Dirac equation
in rectangular coordinates we the driver for the equations is the 4-potential :math:`A`.  One utility
of these equations is to setup a numerical solver for the Dirac equation.

.. image:: latex/dirac.png
    :width: 950px
    :align: center

Manifolds and Tensors
^^^^^^^^^^^^^^^^^^^^^

The example code for tensors on a manifold is

    .. literalinclude:: latex/maxwell.py
        :language: python

This code uses a unit sphere as the example manifold and shows the
geometric derivative and directional derivative operators.  Additonally,
rank 1 and rank 2 general tensors are defined.  The rank 2 tensor is
contracted and evaluated to show multilinear properties.  The covariant
derivatives of both tensors are caculated.  The results are

    .. image:: latex/manifold-0.png
       :width: 700px
       :align: center



    .. image:: latex/manifold-1.png
       :width: 800px
       :align: center

.. rubric:: Citations

.. [Doran]  `<http://www.mrao.cam.ac.uk/~cjld1/pages/book.htm>`_
    ``Geometric Algebra for Physicists`` by C. Doran and A. Lasenby, Cambridge
    University Press, 2003.

.. [Hestenes]  `<http://geocalc.clas.asu.edu/html/CA_to_GC.html>`_
    ``Clifford Algebra to Geometric Calculus`` by D.Hestenes and G. Sobczyk, Kluwer
    Academic Publishers, 1984.

.. [Macdonald] '<http://faculty.luther.edu/~macdonal>'_
   ``Linear and Geometric Algebra`` by Alan Macdonald, `<http://www.amazon.com/Alan-Macdonald/e/B004MB2QJQ>`_


.. rubric:: Footnotes

.. [#f0] In this case :math:`D_{B}^{j_{1}\dots j_{n}} = F` and :math:`\partial_{j_{1}\dots j_{n}} = 1`.

.. [#f1] Since :math:`\underline{T}` is linear we do not require :math:`I^{2} = \pm 1`.

.. [#f2] In this case :math:`y` is a vector in the tangent space and not a coordinate vector so that the
         basis vectors are {\em not} a function of :math:`y`.

.. [#f3] Preprocessing do not occur for the Ipython notebook and the string post processing commands
         *\%* and *\#* are not used in this case
