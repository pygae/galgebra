What is Geometric Algebra?
==========================

.. raw:: html

   <script type="text/x-mathjax-config">
   MathJax.Hub.Config({TeX: { equationNumbers: { autoNumber: "AMS" } }});
   </script>

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

Basics of Geometric Algebra
---------------------------

Geometric algebra is the Clifford algebra of a real finite dimensional vector space or the algebra that results when the vector space is extended with a product of vectors (geometric product) that is associative, left and right distributive, and yields a real number for the square (geometric product) of any vector :cite:`Hestenes`, :cite:`Doran`. The elements of the geometric algebra are called multivectors and consist of the linear combination of scalars, vectors, and the geometric product of
two or more vectors. The additional axioms for the geometric algebra are that for any vectors :math:`a`, :math:`b`, and :math:`c` in the base vector space (:cite:`Doran`,p85):

.. math::

   \begin{array}{c}
     a\lp bc \rp = \lp ab \rp c \\
     a\lp b+c \rp = ab+ac \\
     \lp a + b \rp c = ac+bc \\
     aa = a^{2} \in \Re.
     \end{array}

If the dot (inner) product of two vectors is defined by (:cite:`Doran`,p86)

.. math:: \be a\cdot b \equiv (ab+ba)/2, \ee

then we have

.. math::

   \begin{aligned}
        c &= a+b \\
        c^{2} &= (a+b)^{2} \\
        c^{2} &= a^{2}+ab+ba+b^{2} \\
        a\cdot b &= (c^{2}-a^{2}-b^{2})/2 \in \Re
     \end{aligned}

Thus :math:`a\cdot b` is real. The objects generated from linear combinations of the geometric products of vectors are called multivectors. If a basis for the underlying vector space are the vectors :math:`{\left \{{{{\eb}}_{1},\dots,{{\eb}}_{n}} \rbrc}` (we use boldface :math:`\eb`\ ’s to denote basis vectors) a complete basis for the geometric algebra is given by the scalar :math:`1`, the vectors :math:`{{\eb}}_{1},\dots,{{\eb}}_{n}` and all geometric products of vectors

.. math:: \be {{\eb}}_{i_{1}}{{\eb}}_{i_{2}}\dots {{\eb}}_{i_{r}} \mbox{ where } 0\le r \le n\mbox{, }0 \le i_{j} \le n \mbox{ and } i_{1}<i_{2}<\dots<i_{r} \ee

Each base of the complete basis is represented by a non-commutative symbol (except for the scalar 1) with name :math:`{{\eb}}_{i_{1}}\dots {{\eb}}_{i_{r}}` so that the general multivector :math:`{\boldsymbol{A}}` is represented by (:math:`A` is the scalar part of the multivector and the :math:`A^{i_{1},\dots,i_{r}}` are scalars)

.. math::

   \be {\boldsymbol{A}} = A + \sum_{r=1}^{n}\sum_{\substack{i_{1},\dots,i_{r}\\ 0\le i_{j}<i_{j+1} \le n}}
                  A^{i_{1},\dots,i_{r}}{{\eb}}_{i_{1}}{{\eb}}_{i_{2}}\dots {{\eb}}_{r} \ee

The critical operation in setting up the geometric algebra is reducing the geometric product of any two bases to a linear combination of bases so that we can calculate a multiplication table for the bases. Since the geometric product is associative we can use the operation (by definition for two vectors :math:`a\cdot b \equiv (ab+ba)/2` which is a scalar)

.. math::

   \be \label{reduce}
         {{\eb}}_{i_{j+1}}{{\eb}}_{i_{j}} = 2{{\eb}}_{i_{j+1}}\cdot {{\eb}}_{i_{j}} - {{\eb}}_{i_{j}}{{\eb}}_{i_{j+1}} \ee

These processes are repeated until every basis list in :math:`{\boldsymbol{A}}` is in normal (ascending) order with no repeated elements. As an example consider the following

.. math::

   \begin{aligned}
         {{\eb}}_{3}{{\eb}}_{2}{{\eb}}_{1} &= (2({{\eb}}_{2}\cdot {{\eb}}_{3}) - {{\eb}}_{2}{{\eb}}_{3}){{\eb}}_{1} \\
                         &= 2{\lp {{{\eb}}_{2}\cdot {{\eb}}_{3}} \rp }{{\eb}}_{1} - {{\eb}}_{2}{{\eb}}_{3}{{\eb}}_{1} \\
                         &= 2{\lp {{{\eb}}_{2}\cdot {{\eb}}_{3}} \rp }{{\eb}}_{1} - {{\eb}}_{2}{\lp {2{\lp {{{\eb}}_{1}\cdot {{\eb}}_{3}} \rp }-{{\eb}}_{1}{{\eb}}_{3}} \rp } \\
                         &= 2{\lp {{\lp {{{\eb}}_{2}\cdot {{\eb}}_{3}} \rp }{{\eb}}_{1}-{\lp {{{\eb}}_{1}\cdot {{\eb}}_{3}} \rp }{{\eb}}_{2}} \rp }+{{\eb}}_{2}{{\eb}}_{1}{{\eb}}_{3} \\
                         &= 2{\lp {{\lp {{{\eb}}_{2}\cdot {{\eb}}_{3}} \rp }{{\eb}}_{1}-{\lp {{{\eb}}_{1}\cdot {{\eb}}_{3}} \rp }{{\eb}}_{2}+
                            {\lp {{{\eb}}_{1}\cdot {{\eb}}_{2}} \rp }{{\eb}}_{3}} \rp }-{{\eb}}_{1}{{\eb}}_{2}{{\eb}}_{3}
      \end{aligned}

which results from repeated application of eq. (:math:`\ref{reduce}`). If the product of basis vectors contains repeated factors eq. (:math:`\ref{reduce}`) can be used to bring the repeated factors next to one another so that if :math:`{{\eb}}_{i_{j}} = {{\eb}}_{i_{j+1}}` then :math:`{{\eb}}_{i_{j}}{{\eb}}_{i_{j+1}} = {{\eb}}_{i_{j}}\cdot {{\eb}}_{i_{j+1}}` which is a scalar that commutes with all the terms in the product and can be brought to the front of the product. Since every repeated pair
of vectors in a geometric product of :math:`r` factors reduces the number of non-commutative factors in the product by :math:`r-2`. The number of bases in the multivector algebra is :math:`2^{n}` and the number containing :math:`r` factors is :math:`{n\choose r}` which is the number of combinations or :math:`n` things taken :math:`r` at a time (binomial coefficient).

The other construction required for formulating the geometric algebra is the outer or wedge product (symbol :math:`{\wedge}`) of :math:`r` vectors denoted by :math:`a_{1}{\wedge}\dots{\wedge}a_{r}`. The wedge product of :math:`r` vectors is called an :math:`r`-blade and is defined by (:cite:`Doran`,p86)

.. math:: \be a_{1}{\wedge}\dots{\wedge}a_{r} \equiv \sum_{i_{j_{1}}\dots i_{j_{r}}} \epsilon^{i_{j_{1}}\dots i_{j_{r}}}a_{i_{j_{1}}}\dots a_{i_{j_{1}}} \ee

where :math:`\epsilon^{i_{j_{1}}\dots i_{j_{r}}}` is the contravariant permutation symbol which is :math:`+1` for an even permutation of the superscripts, :math:`0` if any superscripts are repeated, and :math:`-1` for an odd permutation of the superscripts. From the definition :math:`a_{1}{\wedge}\dots{\wedge}a_{r}` is antisymmetric in all its arguments and the following relation for the wedge product of a vector :math:`a` and an :math:`r`-blade :math:`B_{r}` can be derived

.. math::

   \be \label{wedge}
         a{\wedge}B_{r} = (aB_{r}+(-1)^{r}B_{r}a)/2 \ee

Using eq. (:math:`\ref{wedge}`) one can represent the wedge product of all the basis vectors in terms of the geometric product of all the basis vectors so that one can solve (the system of equations is lower diagonal) for the geometric product of all the basis vectors in terms of the wedge product of all the basis vectors. Thus a general multivector :math:`{\boldsymbol{B}}` can be represented as a linear combination of a scalar and the basis blades.

.. math:: \be {\boldsymbol{B}} = B + \sum_{r=1}^{n}\sum_{i_{1},\dots,i_{r},\;\forall\; 0\le i_{j} \le n} B^{i_{1},\dots,i_{r}}{{\eb}}_{i_{1}}{\wedge}{{\eb}}_{i_{2}}{\wedge}\dots{\wedge}{{\eb}}_{r} \ee

Using the blades :math:`{{\eb}}_{i_{1}}{\wedge}{{\eb}}_{i_{2}}{\wedge}\dots{\wedge}{{\eb}}_{r}` creates a graded algebra where :math:`r` is the grade of the basis blades. The grade-:math:`r` part of :math:`{\boldsymbol{B}}` is the linear combination of all terms with grade :math:`r` basis blades.

Grade Projection
~~~~~~~~~~~~~~~~

The scalar part of :math:`{\boldsymbol{B}}` is defined to be grade-:math:`0`. Now that the blade expansion of :math:`{\boldsymbol{B}}` is defined we can also define the grade projection operator :math:`{\left <{{\boldsymbol{B}}} \right >_{r}}` by

.. math:: \be {\left <{{\boldsymbol{B}}} \right >_{r}} = \sum_{i_{1},\dots,i_{r},\;\forall\; 0\le i_{j} \le n} B^{i_{1},\dots,i_{r}}{{\eb}}_{i_{1}}{\wedge}{{\eb}}_{i_{2}}{\wedge}\dots{\wedge}{{\eb}}_{r} \ee

and

.. math:: \be {\left <{{\boldsymbol{B}}} \right >_{}} \equiv {\left <{{\boldsymbol{B}}} \right >_{0}} = B \ee

Multivector Products
~~~~~~~~~~~~~~~~~~~~

Then if :math:`{\boldsymbol{A}}_{r}` is an :math:`r`-grade multivector and :math:`{\boldsymbol{B}}_{s}` is an :math:`s`-grade multivector we have

.. math::

   \be {\boldsymbol{A}}_{r}{\boldsymbol{B}}_{s} = {\left <{{\boldsymbol{A}}_{r}{\boldsymbol{B}}_{s}} \right >_{{\left |{r-s}\right |}}}+{\left <{{\boldsymbol{A}}_{r}{\boldsymbol{B}}_{s}} \right >_{{\left |{r-s}\right |}+2}}+\cdots
                                {\left <{{\boldsymbol{A}}_{r}{\boldsymbol{B}}_{s}} \right >_{r+s}} \ee

and define (:cite:`Hestenes`,p6)

.. math::

   \begin{aligned}
         {\boldsymbol{A}}_{r}{\wedge}{\boldsymbol{B}}_{s} &\equiv {\left <{{\boldsymbol{A}}_{r}{\boldsymbol{B}}_{s}} \right >_{r+s}} \\
         {\boldsymbol{A}}_{r}\cdot{\boldsymbol{B}}_{s} &\equiv {\left \{ { \begin{array}{cc}
         r\mbox{ and }s \ne 0: & {\left <{{\boldsymbol{A}}_{r}{\boldsymbol{B}}_{s}} \right >_{{\left |{r-s}\right |}}}  \\
         r\mbox{ or }s = 0: & 0 \end{array}} \right \}}
      \end{aligned}

where :math:`{\boldsymbol{A}}_{r}\cdot{\boldsymbol{B}}_{s}` is called the dot or inner product of two pure grade multivectors. For the case of two non-pure grade multivectors

.. math::

   \begin{aligned}
         {\boldsymbol{A}}{\wedge}{\boldsymbol{B}} &= \sum_{r,s}{\left <{{\boldsymbol{A}}} \right >_{r}}{\wedge}{\left <{{\boldsymbol{B}}} \right >_{{s}}} \\
         {\boldsymbol{A}}\cdot{\boldsymbol{B}} &= \sum_{r,s\ne 0}{\left <{{\boldsymbol{A}}} \right >_{r}}\cdot{\left <{{\boldsymbol{B}}} \right >_{{s}}}
      \end{aligned}

Two other products, the left (:math:`\rfloor`) and right (:math:`\lfloor`) contractions, are defined by

.. math::

   \begin{aligned}
         {\boldsymbol{A}}\lfloor{\boldsymbol{B}} &\equiv \sum_{r,s}{\left \{ {\begin{array}{cc} {\left <{{\boldsymbol{A}}_r{\boldsymbol{B}}_{s}} \right >_{r-s}} & r \ge s \\
                                                     0                                               & r < s \end{array}} \right \}}  \\
         {\boldsymbol{A}}\rfloor{\boldsymbol{B}} &\equiv \sum_{r,s}{\left \{ {\begin{array}{cc} {\left <{{\boldsymbol{A}}_{r}{\boldsymbol{B}}_{s}} \right >_{s-r}} & s \ge r \\
                                                     0                                               & s < r\end{array}} \right \}}
      \end{aligned}

Reverse of Multivector
~~~~~~~~~~~~~~~~~~~~~~

A final operation for multivectors is the reverse. If a multivector :math:`{\boldsymbol{A}}` is the geometric product of :math:`r` vectors (versor) so that :math:`{\boldsymbol{A}} = a_{1}\dots a_{r}` the reverse is defined by

.. math::

   \begin{aligned}
         {\boldsymbol{A}}^{{\dagger}} \equiv a_{r}\dots a_{1}
      \end{aligned}

where for a general multivector we have (the the sum of the reverse of versors)

.. math:: \be {\boldsymbol{A}}^{{\dagger}} = A + \sum_{r=1}^{n}(-1)^{r(r-1)/2}\sum_{i_{1},\dots,i_{r},\;\forall\; 0\le i_{j} \le n} A^{i_{1},\dots,i_{r}}{{\eb}}_{i_{1}}{\wedge}{{\eb}}_{i_{2}}{\wedge}\dots{\wedge}{{\eb}}_{r} \ee

note that if :math:`{\boldsymbol{A}}` is a versor then :math:`{\boldsymbol{A}}{\boldsymbol{A}}^{{\dagger}}\in\Re` and (:math:`AA^{{\dagger}} \ne 0`)

.. math:: \be {\boldsymbol{A}}^{-1} = {\displaystyle\frac{{\boldsymbol{A}}^{{\dagger}}}{{\boldsymbol{AA}}^{{\dagger}}}} \ee

The reverse is important in the theory of rotations in :math:`n`-dimensions. If :math:`R` is the product of an even number of vectors and :math:`RR^{{\dagger}} = 1` then :math:`RaR^{{\dagger}}` is a composition of rotations of the vector :math:`a`. If :math:`R` is the product of two vectors then the plane that :math:`R` defines is the plane of the rotation. That is to say that :math:`RaR^{{\dagger}}` rotates the component of :math:`a` that is projected into the plane defined by :math:`a` and
:math:`b` where :math:`R=ab`. :math:`R` may be written :math:`R = e^{\frac{\theta}{2}U}`, where :math:`\theta` is the angle of rotation and :math:`U` is a unit blade :math:`\lp U^{2} = \pm 1\rp` that defines the plane of rotation.

Reciprocal Frames
~~~~~~~~~~~~~~~~~

If we have :math:`M` linearly independent vectors (a frame), :math:`a_{1},\dots,a_{M}`, then the reciprocal frame is :math:`a^{1},\dots,a^{M}` where :math:`a_{i}\cdot a^{j} = \delta_{i}^{j}`, :math:`\delta_{i}^{j}` is the Kronecker delta (zero if :math:`i \ne j` and one if :math:`i = j`). The reciprocal frame is constructed as follows:

.. math:: \be E_{M} = a_{1}{\wedge}\dots{\wedge}a_{M} \ee

.. math:: \be E_{M}^{-1} = {\displaystyle\frac{E_{M}}{E_{M}^{2}}} \ee

Then

.. math:: \be a^{i} = \lp -1\rp ^{i-1}\lp a_{1}{\wedge}\dots{\wedge}\breve{a}_{i} {\wedge}\dots{\wedge}a_{M}\rp E_{M}^{-1} \ee

where :math:`\breve{a}_{i}` indicates that :math:`a_{i}` is to be deleted from the product. In the standard notation if a vector is denoted with a subscript the reciprocal vector is denoted with a superscript. The set of reciprocal vectors will be calculated if a coordinate set is given when a geometric algebra is instantiated since they are required for geometric differentiation when the ``Ga`` member function ``Ga.mvr()`` is called to return the reciprocal basis in terms of the basis vectors.

.. _sect_manifold:

Manifolds and Submanifolds
--------------------------

A :math:`m`-dimensional vector manifold\ [4]_, :math:`\mathcal{M}`, is defined by a coordinate tuple (tuples are indicated by the vector accent “:math:`\vec{\;\;\;}`”)

.. math:: \be \vec{x} = \paren{x^{1},\dots,x^{m}}, \ee

and the differentiable mapping (:math:`U^{m}` is an :math:`m`-dimensional subset of :math:`\Re^{m}`)

.. math:: \be \f{\bm{e}^{\mathcal{M}}}{\vec{x}}\colon U^{m}\subseteq\Re^{m}\rightarrow \mathcal{V}, \ee

where :math:`\mathcal{V}` is a vector space with an inner product\ [5]_ (:math:`\cdot`) and is of :math:`{{\dim}\lp {\mathcal{V}} \rp } \ge m`.

Then a set of basis vectors for the tangent space of :math:`\mathcal{M}` at :math:`\vec{x}`, :math:`{{{\mathcal{T}_{\vec{x}}}\lp {\mathcal{M}} \rp }}`, are

.. math:: \be \bm{e}_{i}^{\mathcal{M}} = \pdiff{\bm{e}^{\mathcal{M}}}{x^{i}} \ee

and

.. math:: \be \f{g_{ij}^{\mathcal{M}}}{\vec{x}} = \bm{e}_{i}^{\mathcal{M}}\cdot\bm{e}_{j}^{\mathcal{M}}. \ee

A :math:`n`-dimensional (:math:`n\le m`) submanifold :math:`\mathcal{N}` of :math:`\mathcal{M}` is defined by a coordinate tuple

.. math:: \be \vec{u} = \paren{u^{1},\dots,u^{n}}, \ee

and a differentiable mapping

.. math::

   \be \label{eq_79}
       \f{\vec{x}}{\vec{u}}\colon U^{n}\subseteq\Re^{n}\rightarrow U^{m}\subseteq\Re^{m},
    \ee

Then the basis vectors for the tangent space :math:`{{{\mathcal{T}_{\vec{u}}}\lp {\mathcal{N}} \rp }}` are (using :math:`{{{{\eb}}^{\mathcal{N}}}\lp {\vec{u}} \rp } = {{{{\eb}}^{\mathcal{M}}}\lp {{{\vec{x}}\lp {\vec{u}} \rp }} \rp }` and the chain rule)\ [6]_

.. math::

   \be     \f{\bm{e}_{i}^{\mathcal{N}}}{\vec{u}} = \pdiff{\f{\bm{e}^{\mathcal{N}}}{\vec{u}}}{u^{i}}
                                                 = \pdiff{\f{\bm{e}^{\mathcal{M}}}{\vec{x}}}{x^{j}}\pdiff{x^{j}}{u^{i}}
                                                 = \f{\bm{e}_{j}^{\mathcal{M}}}{\f{\vec{x}}{\vec{u}}}\pdiff{x^{j}}{u^{i}}, \ee

and

.. math::

   \be \label{eq_81}
       \f{g_{ij}^{\mathcal{N}}}{\vec{u}} = \pdiff{x^{k}}{u^{i}}\pdiff{x^{l}}{u^{j}}
                                               \f{g_{kl}^{\mathcal{M}}}{\f{\vec{x}}{\vec{u}}}.
    \ee

Going back to the base manifold, :math:`\mathcal{M}`, note that the mapping :math:`{{{\eb}^{\mathcal{M}}}\lp {\vec{x}} \rp }\colon U^{n}\subseteq\Re^{n}\rightarrow \mathcal{V}` allows us to calculate an unnormalized pseudo-scalar for :math:`{{{\mathcal{T}_{\vec{x}}}\lp {\mathcal{M}} \rp }}`,

.. math::

   \be     \f{I^{\mathcal{M}}}{\vec{x}} = \f{\bm{e}_{1}^{\mathcal{M}}}{\vec{x}}
                                          \W\dots\W\f{\bm{e}_{m}^{\mathcal{M}}}{\vec{x}}. \ee

With the pseudo-scalar we can define a projection operator from :math:`\mathcal{V}` to the tangent space of :math:`\mathcal{M}` by

.. math::

   \be     \f{P_{\vec{x}}}{\bm{v}} = (\bm{v}\cdot \f{I^{\mathcal{M}}}{\vec{x}})
                                 \paren{\f{I^{\mathcal{M}}}{\vec{x}}}^{-1} \;\forall\; \bm{v}\in\mathcal{V}. \ee

In fact for each tangent space :math:`{{{\mathcal{T}_{\vec{x}}}\lp {\mathcal{M}} \rp }}` we can define a geometric algebra :math:`{{\mathcal{G}}\lp {{{{\mathcal{T}_{\vec{x}}}\lp {\mathcal{M}} \rp }}} \rp }` with pseudo-scalar :math:`I^{\mathcal{M}}` so that if :math:`A \in {{\mathcal{G}}\lp {\mathcal{V}} \rp }` then

.. math::

   \be     \f{P_{\vec{x}}}{A} = \paren{A\cdot \f{I^{\mathcal{M}}}{\vec{x}}}
                            \paren{\f{I^{\mathcal{M}}}{\vec{x}}}^{-1}
                            \in \f{\mathcal{G}}{\Tn{\mathcal{M}}{\vec{x}}}\;\forall\;
                            A \in \f{\mathcal{G}}{\mathcal{V}} \ee

and similarly for the submanifold :math:`\mathcal{N}`.

If the embedding :math:`{{{\eb}^{\mathcal{M}}}\lp {\vec{x}} \rp }\colon U^{n}\subseteq\Re^{n}\rightarrow \mathcal{V}` is not given, but the metric tensor :math:`{{g_{ij}^{\mathcal{M}}}\lp {\vec{x}} \rp }` is given the geometric algebra of the tangent space can be constructed. Also the derivatives of the basis vectors of the tangent space can be calculated from the metric tensor using the Christoffel symbols, :math:`{{\Gamma_{ij}^{k}}\lp {\vec{u}} \rp }`, where the derivatives of the basis
vectors are given by

.. math:: \be \pdiff{\bm{e}_{j}^{\mathcal{M}}}{x^{i}} =\f{\Gamma_{ij}^{k}}{\vec{u}}\bm{e}_{k}^{\mathcal{M}}. \ee

If we have a submanifold, :math:`\mathcal{N}`, defined by eq. (:math:`\ref{eq_79}`) we can calculate the metric of :math:`\mathcal{N}` from eq. (:math:`\ref{eq_81}`) and hence construct the geometric algebra and calculus of the tangent space, :math:`{{{\mathcal{T}_{\vec{u}}}\lp {\mathcal{N}} \rp }}\subseteq {{{\mathcal{T}_{{{\vec{x}}\lp {\vec{u}} \rp }}}\lp {\mathcal{M}} \rp }}`.

**Note:**

If the base manifold is normalized (use the hat symbol to denote normalized tangent vectors, :math:`\hat{{\eb}}_{i}^{\mathcal{M}}`, and the resulting metric tensor, :math:`\hat{g}_{ij}^{\mathcal{M}}`) we have :math:`\hat{{\eb}}_{i}^{\mathcal{M}}\cdot\hat{{\eb}}_{i}^{\mathcal{M}} = \pm 1` and :math:`\hat{g}_{ij}^{\mathcal{M}}` does not posses enough information to calculate :math:`g_{ij}^{\mathcal{N}}`. In that case we need to know :math:`g_{ij}^{\mathcal{M}}`, the metric tensor of the base
manifold before normalization. Likewise, for the case of a vector manifold unless the mapping, :math:`{{{\eb}^{\mathcal{M}}}\lp {\vec{x}} \rp }\colon U^{m}\subseteq\Re^{m}\rightarrow \mathcal{V}`, is constant the tangent vectors and metric tensor can only be normalized after the fact (one cannot have a mapping that automatically normalizes all the tangent vectors).

Geometric Derivative
--------------------

The directional derivative of a multivector field :math:`{{F}\lp {x} \rp }` is defined by (:math:`a` is a vector and :math:`h` is a scalar)

.. math:: \be \paren{a\cdot\nabla_{x}}F \equiv \lim_{h\rightarrow 0}\bfrac{\f{F}{x+ah}-\f{F}{x}}{h}. \label{eq_50} \ee

Note that :math:`a\cdot\nabla_{x}` is a scalar operator. It will give a result containing only those grades that are already in :math:`F`. :math:`{\lp {a\cdot\nabla_{x}} \rp }F` is the best linear approximation of :math:`{{F}\lp {x} \rp }` in the direction :math:`a`. Equation (:math:`\ref{eq_50}`) also defines the operator :math:`\nabla_{x}` which for the basis vectors, :math:`{\left \{{{\eb}_{i}} \rbrc}`, has the representation (note that the :math:`{\left \{{{\eb}^{j}} \rbrc}` are reciprocal
basis vectors)

.. math:: \be \nabla_{x} F = {\eb}^{j}{\displaystyle\frac{\partial F}{\partial x^{j}}} \ee

If :math:`F_{r}` is a :math:`r`-grade multivector (if the independent vector, :math:`x`, is obvious we suppress it in the notation and just write :math:`\nabla`) and :math:`F_{r} = F_{r}^{i_{1}\dots i_{r}}{\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}}` then

.. math:: \be \nabla F_{r} = {\displaystyle\frac{\partial F_{r}^{i_{1}\dots i_{r}}}{\partial x^{j}}}{\eb}^{j}\lp {\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}} \rp  \ee

Note that :math:`{\eb}^{j}\lp {\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}} \rp` can only contain grades :math:`r-1` and :math:`r+1` so that :math:`\nabla F_{r}` also can only contain those grades. For a grade-:math:`r` multivector :math:`F_{r}` the inner (div) and outer (curl) derivatives are

.. math:: \be \nabla\cdot F_{r} = \left < \nabla F_{r}\right >_{r-1} = {\eb}^{j}\cdot {{\displaystyle\frac{\partial {F_{r}}}{\partial {x^{j}}}}} \ee

and

.. math:: \be \nabla{\wedge}F_{r} = \left < \nabla F_{r}\right >_{r+1} = {\eb}^{j}{\wedge}{{\displaystyle\frac{\partial {F_{r}}}{\partial {x^{j}}}}} \ee

For a general multivector function :math:`F` the inner and outer derivatives are just the sum of the inner and outer derivatives of each grade of the multivector function.

Geometric Derivative on a Manifold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the case of a manifold the derivatives of the :math:`{\eb}_{i}`\ ’s are functions of the coordinates, :math:`{\left \{{x^{i}} \rbrc}`, so that the geometric derivative of a :math:`r`-grade multivector field is

.. math::

   \begin{aligned}
       \nabla F_{r} &= {\eb}^{i}{{\displaystyle\frac{\partial {F_{r}}}{\partial {x^{i}}}}} = {\eb}^{i}{{\displaystyle\frac{\partial {}}{\partial {x^{i}}}}}
                      {\lp {F_{r}^{i_{1}\dots i_{r}}{\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}}} \rp } \nonumber \\
                    &= {{\displaystyle\frac{\partial {F_{r}^{i_{1}\dots i_{r}}}}{\partial {x^{i}}}}}{\eb}^{i}{\lp {{\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}}} \rp }
                       +F_{r}^{i_{1}\dots i_{r}}{\eb}^{i}{{\displaystyle\frac{\partial {}}{\partial {x^{i}}}}}{\lp {{\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}}} \rp }\end{aligned}

where the multivector functions :math:`{\eb}^{i}{{\displaystyle\frac{\partial {}}{\partial {x^{i}}}}}{\lp {{\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}}} \rp }` are the connection for the manifold\ [7]_.

The directional (material/convective) derivative, :math:`{\lp {v\cdot\nabla} \rp }F_{r}` is given by

.. math::

   \begin{aligned}
       {\lp {v\cdot\nabla} \rp } F_{r} &= v^{i}{{\displaystyle\frac{\partial {F_{r}}}{\partial {x^{i}}}}} = v^{i}{{\displaystyle\frac{\partial {}}{\partial {x^{i}}}}}
                      {\lp {F_{r}^{i_{1}\dots i_{r}}{\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}}} \rp } \nonumber \\
                    &= v^{i}{{\displaystyle\frac{\partial {F_{r}^{i_{1}\dots i_{r}}}}{\partial {x^{i}}}}}{\lp {{\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}}} \rp }
                       +v^{i}F_{r}^{i_{1}\dots i_{r}}{{\displaystyle\frac{\partial {}}{\partial {x^{i}}}}}{\lp {{\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}}} \rp },\end{aligned}

so that the multivector connection functions for the directional derivative are :math:`{{\displaystyle\frac{\partial {}}{\partial {x^{i}}}}}{\lp {{\eb}_{i_{1}}{\wedge}\dots{\wedge}{\eb}_{i_{r}}} \rp }`. Be careful and note that :math:`{\lp {v\cdot\nabla} \rp } F_{r} \ne v\cdot {\lp {\nabla F_{r}} \rp }` since the dot and geometric products are not associative with respect to one another (:math:`v\cdot\nabla` is a scalar operator).

Normalizing Basis for Derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The basis vector set, :math:`{\left \{ {{\eb}_{i}} \rbrc}`, is not in general normalized. We define a normalized set of basis vectors, :math:`{\left \{{{\boldsymbol{\hat{e}}}_{i}} \rbrc}`, by

.. math:: \be {\boldsymbol{\hat{e}}}_{i} = {\displaystyle\frac{{\eb}_{i}}{\sqrt{{\left |{{\eb}_{i}^{2}}\right |}}}} = {\displaystyle\frac{{\eb}_{i}}{{\left |{{\eb}_{i}}\right |}}}. \ee

This works for all :math:`{\eb}_{i}^{2} \neq 0`. Note that :math:`{\boldsymbol{\hat{e}}}_{i}^{2} = \pm 1`.

Thus the geometric derivative for a set of normalized basis vectors is (where :math:`F_{r} = F_{r}^{i_{1}\dots i_{r}} \bm{\hat{e}}_{i_{1}}\W\dots\W\bm{\hat{e}}_{i_{r}}` and [no summation] :math:`\hat{F}_{r}^{i_{1}\dots i_{r}} = F_{r}^{i_{1}\dots i_{r}} \abs{\bm{\hat{e}}_{i_{1}}}\dots\abs{\bm{\hat{e}}_{i_{r}}}`).

.. math::

   \be     \nabla F_{r} = \eb^{i}\pdiff{F_{r}}{x^{i}} =
                      \pdiff{F_{r}^{i_{1}\dots i_{r}}}{x^{i}}\bm{e}^{i}
                      \paren{\bm{\hat{e}}_{i_{1}}\W\dots\W\bm{\hat{e}}_{i_{r}}}
                       +F_{r}^{i_{1}\dots i_{r}}\bm{e}^{i}\pdiff{}{x^{i}}
                       \paren{\bm{\hat{e}}_{i_{1}}\W\dots\W\bm{\hat{e}}_{i_{r}}}. \ee

To calculate :math:`{\eb}^{i}` in terms of the :math:`{\boldsymbol{\hat{e}}}_{i}`\ ’s we have

.. math::

   \begin{aligned}
       {\eb}^{i} &= g^{ij}{\eb}_{j} \nonumber \\
       {\eb}^{i} &= g^{ij}{\left |{{\eb}_{j}}\right |}{\boldsymbol{\hat{e}}}_{j}.\end{aligned}

This is the general (non-orthogonal) formula. If the basis vectors are orthogonal then (no summation over repeated indexes)

.. math::

   \begin{aligned}
       {\eb}^{i} &= g^{ii}{\left |{{\eb}_{i}}\right |}{\boldsymbol{\hat{e}}}_{i} \nonumber \\
       {\eb}^{i} &= {\displaystyle\frac{{\left |{{\eb}_{i}}\right |}}{g_{ii}}}{\boldsymbol{\hat{e}}}_{i} = {\displaystyle\frac{{\left |{{\boldsymbol{\hat{e}}}_{i}}\right |}}{{\eb}_{i}^{2}}}{\boldsymbol{\hat{e}}}_{i}.\end{aligned}

Additionally, one can calculate the connection of the normalized basis as follows

.. math::

   \begin{aligned}
       {{\displaystyle\frac{\partial {{\lp {{\left |{{\eb}_{i}}\right |}{\boldsymbol{\hat{e}}}_{i}} \rp }}}{\partial {x^{j}}}}} =& {{\displaystyle\frac{\partial {{\eb}_{i}}}{\partial {x^{j}}}}}, \nonumber \\
       {{\displaystyle\frac{\partial {{\left |{{\eb}_{i}}\right |}}}{\partial {x^{j}}}}}{\boldsymbol{\hat{e}}}_{i}
                                         +{\left |{{\eb}_{i}}\right |}{{\displaystyle\frac{\partial {{\boldsymbol{\hat{e}}}_{i}}}{\partial {x^{j}}}}} =& {{\displaystyle\frac{\partial {{\eb}_{i}}}{\partial {x^{j}}}}}, \nonumber \\
       {{\displaystyle\frac{\partial {{\boldsymbol{\hat{e}}}_{i}}}{\partial {x^{j}}}}} =& {\displaystyle\frac{1}{{\left |{{\eb}_{i}}\right |}}}{\lp {{{\displaystyle\frac{\partial {{\eb}_{i}}}{\partial {x^{j}}}}}
                                          -{{\displaystyle\frac{\partial {{\left |{{\eb}_{i}}\right |}}}{\partial {x^{j}}}}}{\boldsymbol{\hat{e}}}_{i}} \rp },\nonumber \\
                                       =& {\displaystyle\frac{1}{{\left |{{\eb}_{i}}\right |}}}{{\displaystyle\frac{\partial {{\eb}_{i}}}{\partial {x^{j}}}}}
                                          -{\displaystyle\frac{1}{{\left |{{\eb}_{i}}\right |}}}{{\displaystyle\frac{\partial {{\left |{{\eb}_{i}}\right |}}}{\partial {x^{j}}}}}{\boldsymbol{\hat{e}}}_{i},\nonumber \\
                                       =& {\displaystyle\frac{1}{{\left |{{\eb}_{i}}\right |}}}{{\displaystyle\frac{\partial {{\eb}_{i}}}{\partial {x^{j}}}}}
                                          -{\displaystyle\frac{1}{2g_{ii}}}{{\displaystyle\frac{\partial {g_{ii}}}{\partial {x^{j}}}}}{\boldsymbol{\hat{e}}}_{i},\end{aligned}

where :math:`{{\displaystyle\frac{\partial {{\eb}_{i}}}{\partial {x^{j}}}}}` is expanded in terms of the :math:`{\boldsymbol{\hat{e}}}_{i}`\ ’s.

.. _ldops:

Linear Differential Operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First a note on partial derivative notation. We shall use the following notation for a partial derivative where the manifold coordinates are :math:`x_{1},\dots,x_{n}`:

.. math::

   \be\label{eq_66a}
       \bfrac{\partial^{j_{1}+\cdots+j_{n}}}{\partial x_{1}^{j_{1}}\dots\partial x_{n}^{j_{n}}} = \partial_{j_{1}\dots j_{n}}.
   \ee

If :math:`j_{k}=0` the partial derivative with respect to the :math:`k^{th}` coordinate is not taken. If :math:`j_{k} = 0` for all :math:`1 \le k \le n` then the partial derivative operator is the scalar one. If we consider a partial derivative where the :math:`x`\ ’s are not in normal order such as

.. math:: \be {\displaystyle\frac{\partial^{j_{1}+\cdots+j_{n}}}{\partial x_{i_{1}}^{j_{1}}\dots\partial x_{i_{n}}^{j_{n}}}}, \ee

and the :math:`i_{k}`\ ’s are not in ascending order. The derivative can always be put in the form in eq (:math:`\ref{eq_66a}`) since the order of differentiation does not change the value of the partial derivative (for the smooth functions we are considering). Additionally, using our notation the product of two partial derivative operations is given by

.. math:: \be \partial_{i_{1}\dots i_{n}}\partial_{j_{1}\dots j_{n}} = \partial_{i_{1}+j_{1},\dots, i_{n}+j_{n}}. \ee

A general general multivector linear differential operator is a linear combination of multivectors and partial derivative operators denoted by

.. math::

   \be\label{eq_66b}
       D \equiv D^{i_{1}\dots i_{n}}\partial_{i_{1}\dots i_{n}}.
   \ee

Equation (:math:`\ref{eq_66b}`) is the normal form of the differential operator in that the partial derivative operators are written to the right of the multivector coefficients and do not operate upon the multivector coefficients. The operator of eq (:math:`\ref{eq_66b}`) can operate on mulitvector functions, returning a multivector function via the following definitions.

:math:`F` as

.. math:: \be D\circ F = D^{j_{1}\dots j_{n}}\circ\partial_{j_{1}\dots j_{n}}F,\label{eq_67a}  \ee

, or

.. math:: \be F\circ D = \partial_{j_{1}\dots j_{n}}F\circ D^{j_{1}\dots j_{n}},\label{eq_68a} \ee

where the :math:`D^{j_{1}\dots j_{n}}` are multivector functions and :math:`\circ` is any of the multivector multiplicative operations.

Equations (:math:`\ref{eq_67a}`) and (:math:`\ref{eq_68a}`) are not the most general multivector linear differential operators, the most general would be

.. math:: \be D \left( F \right) = {D^{j_{1}\dots j_{n}}}\left({\partial_{j_{1}\dots j_{n}}F}\right), \ee

where :math:`{{D^{j_{1}\dots j_{n}}}\lp {} \rp }` are linear multivector functionals.

The definition of the sum of two differential operators is obvious since any multivector operator, :math:`\circ`, is a bilinear operator :math:`{\lp {{\lp {D_{A}+D_{B}} \rp }\circ F = D_{A}\circ F+D_{B}\circ F} \rp }`, the product of two differential operators :math:`D_{A}` and :math:`D_{B}` operating on a multivector function :math:`F` is defined to be (:math:`\circ_{1}` and :math:`\circ_{2}` are any two multivector multiplicative operations)

.. math::

   \begin{aligned}
       {\lp {D_{A}\circ_{1}D_{B}} \rp }\circ_{2}F &\equiv {\lp {D_{A}^{i_{1}\dots i_{n}}\circ_{1}
                                                     \partial_{i_{1}\dots i_{n}}{\lp {D_{B}^{j_{1}\dots j_{n}}
                                                     \partial_{j_{1}\dots j_{n}}} \rp }} \rp }\circ_{2}F \nonumber \\
                                             &= {\lp {D_{A}^{i_{1}\dots i_{n}}\circ_{1}
                                                {\lp {{\lp {\partial_{i_{1}\dots i_{n}}D_{B}^{j_{1}\dots j_{n}}} \rp }
                                                \partial_{j_{1}\dots j_{n}}+
                                                D_{B}^{j_{1}\dots j_{n}}} \rp }
                                                \partial_{i_{1}+j_{1},\dots, i_{n}+j_{n}}} \rp }\circ_{2}F \nonumber \\
                                             &= {\lp {D_{A}^{i_{1}\dots i_{n}}\circ_{1}{\lp {\partial_{i_{1}\dots i_{n}}D_{B}^{j_{1}\dots j_{n}}} \rp }} \rp }
                                                \circ_{2}\partial_{j_{1}\dots j_{n}}F+
                                                {\lp {D_{A}^{i_{1}\dots i_{n}}\circ_{1}D_{B}^{j_{1}\dots j_{n}}} \rp }
                                                \circ_{2}\partial_{i_{1}+j_{1},\dots, i_{n}+j_{n}}F,\end{aligned}

where we have used the fact that the :math:`\partial` operator is a scalar operator and commutes with :math:`\circ_{1}` and :math:`\circ_{2}`.

Thus for a pure operator product :math:`D_{A}\circ D_{B}` we have

.. math::

   \be D_{A}\circ D_{B} = \paren{D_{A}^{i_{1}\dots i_{n}}\circ\paren{\partial_{i_{1}\dots i_{n}}D_{B}^{j_{1}\dots j_{n}}}}
                                                \partial_{j_{1}\dots j_{n}}+
                                                \paren{D_{A}^{i_{1}\dots i_{n}}\circ_{1}D_{B}^{j_{1}\dots j_{n}}}
                                                \partial_{i_{1}+j_{1},\dots, i_{n}+j_{n}} \label{eq_71a}  \ee

and the form of eq (:math:`\ref{eq_71a}`) is the same as eq (:math:`\ref{eq_67a}`). The basis of eq (:math:`\ref{eq_71a}`) is that the :math:`\partial` operator operates on all object to the right of it as products so that the product rule must be used in all differentiations. Since eq (:math:`\ref{eq_71a}`) puts the product of two differential operators in standard form we also evaluate :math:`F\circ_{2}{\lp {D_{A}\circ_{1}D_{B}} \rp }`.

We now must distinguish between the following cases. If :math:`D` is a differential operator and :math:`F` a multivector function should :math:`D\circ F` and :math:`F\circ D` return a differential operator or a multivector. In order to be consistent with the standard vector analysis we have :math:`D\circ F` return a multivector and :math:`F\circ D` return a differential operator. Then we define the complementary differential operator :math:`\bar{D}` which is identical to :math:`D` except that
:math:`\bar{D}\circ F` returns a differential operator according to eq (:math:`\ref{eq_71a}`)\ [8]_ and :math:`F\circ\bar{D}` returns a multivector according to eq (:math:`\ref{eq_68a}`).

A general differential operator is built from repeated applications of the basic operator building blocks :math:`{\lp {\bar{\nabla}\circ A} \rp }`, :math:`{\lp {A\circ\bar{\nabla}} \rp }`, :math:`{\lp {\bar{\nabla}\circ\bar{\nabla}} \rp }`, and :math:`{\lp {A\pm \bar{\nabla}} \rp }`. Both :math:`\nabla` and :math:`\bar{\nabla}` are represented by the operator

.. math::

   \be 
       \nabla = \bar{\nabla} = e^{i}\pdiff{}{x^{i}},
    \ee

but are flagged to produce the appropriate result.

In the our notation the directional derivative operator is :math:`a\cdot\nabla`, the Laplacian :math:`\nabla\cdot\nabla` and the expression for the Riemann tensor, :math:`R^{i}_{jkl}`, is

.. math:: \be \paren{\nabla\W\nabla}\eb^{i} = \half R^{i}_{jkl}\paren{\eb^{j}\W\eb^{k}}\eb^{l}. \ee

We would use the complement if we wish a quantum mechanical type commutator defining

.. math::

   \be
       \com{x,\nabla} \equiv x\nabla - \bar{\nabla}x,
   \ee

, or if we wish to simulate the dot notation (Doran and Lasenby)

.. math::

   \be
       \dot{F}\dot{\nabla} = F\bar{\nabla}.
   \ee

Split Differential Operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To implement the general “dot” notation for differential operators in python is not possible. Another type of symbolic notation is required. I propose what one could call the “split differential operator.” For :math:`\nabla` denote the corresponding split operator by two operators :math:`{{\nabla}_{\mathcal{G}}}` and :math:`{{\nabla}_{\mathcal{D}}}` where in practice :math:`{{\nabla}_{\mathcal{G}}}` is a tuple of vectors and :math:`{{\nabla}_{\mathcal{D}}}` is a tuple of corresponding partial
derivatives. Then the equivalent of the “dot” notation would be

.. math:: \be \dot{\nabla}{\lp {A\dot{B}C} \rp } = {{\nabla}_{\mathcal{G}}}{\lp {A{\lp {{{\nabla}_{\mathcal{D}}}B} \rp }C} \rp }.\label{splitopV} \ee

We are using the :math:`\mathcal{G}` subscript to indicate the geometric algebra parts of the multivector differential operator and the :math:`\mathcal{D}` subscript to indicate the scalar differential operator parts of the multivector differential operator. An example of this notation in 3D Euclidean space is

.. math::

   \begin{aligned}
       {{\nabla}_{\mathcal{G}}} &= {\lp {{{\eb}}_{x},{{\eb}}_{y},{{\eb}}_{z}} \rp }, \\
       {{\nabla}_{\mathcal{D}}} &= {\lp {{{\displaystyle\frac{\partial {}}{\partial {x}}}},{{\displaystyle\frac{\partial {}}{\partial {y}}}},{{\displaystyle\frac{\partial {}}{\partial {x}}}}} \rp },\end{aligned}

To implement :math:`{{\nabla}_{\mathcal{G}}}` and :math:`{{\nabla}_{\mathcal{D}}}` we have in the example

.. math::

   \begin{aligned}
       {{\nabla}_{\mathcal{D}}}B &= {\lp {{{\displaystyle\frac{\partial {B}}{\partial {x}}}},{{\displaystyle\frac{\partial {B}}{\partial {y}}}},{{\displaystyle\frac{\partial {B}}{\partial {z}}}}} \rp } \\
       {\lp {{{\nabla}_{\mathcal{D}}}B} \rp }C &= {\lp {{{\displaystyle\frac{\partial {B}}{\partial {x}}}}C,{{\displaystyle\frac{\partial {B}}{\partial {y}}}}C,{{\displaystyle\frac{\partial {B}}{\partial {z}}}}C} \rp } \\
       A{\lp {{{\nabla}_{\mathcal{D}}}B} \rp }C &= {\lp {A{{\displaystyle\frac{\partial {B}}{\partial {x}}}}C,A{{\displaystyle\frac{\partial {B}}{\partial {y}}}}C,A{{\displaystyle\frac{\partial {B}}{\partial {z}}}}C} \rp }.\end{aligned}

Then the final evaluation is

.. math:: \be {{\nabla}_{\mathcal{G}}}{\lp {A{\lp {{{\nabla}_{\mathcal{D}}}B} \rp }C} \rp } = {{\eb}}_{x}A{{\displaystyle\frac{\partial {B}}{\partial {x}}}}C+{{\eb}}_{y}A{{\displaystyle\frac{\partial {B}}{\partial {y}}}}C+{{\eb}}_{z}A{{\displaystyle\frac{\partial {B}}{\partial {z}}}}C, \ee

which could be called the “dot” product of two tuples. Note that :math:`\nabla = {{\nabla}_{\mathcal{G}}}{{\nabla}_{\mathcal{D}}}` and :math:`\dot{F}\dot{\nabla} = F\bar{\nabla} = {\lp {{{\nabla}_{\mathcal{D}}}F} \rp }{{\nabla}_{\mathcal{G}}}`.

For the general multivector differential operator, :math:`D`, the split operator parts are :math:`{{D}_{\mathcal{G}}}`, a tuple of basis blade multivectors and :math:`{{D}_{\mathcal{D}}}`, a tuple of scalar differential operators that correspond to the coefficients of the basis-blades in the total operator :math:`D` so that

.. math:: \be \dot{D}{\lp {A\dot{B}C} \rp } = {{D}_{\mathcal{G}}}{\lp {A{\lp {{{D}_{\mathcal{D}}}B} \rp }C} \rp }. \label{splitopM} \ee

If the index set for the basis blades of a geometric algebra is denoted by :math:`{\left \{{n} \rbrc}` where :math:`{\left \{{n} \rbrc}` contains :math:`2^{n}` indices for an :math:`n` dimensional geometric algebra then the most general multivector differential operator can be written\ [9]_

.. math::

   \be D = {{\displaystyle}\sum_{l\in{\left \{
   {n} \rbrc}}{{\eb}}^{l}D_{{\left \{
   {l} \rbrc}}} \ee

.. math::

   \be \dot{D}{\lp {A\dot{B}C} \rp } = {{D}_{\mathcal{G}}}{\lp {A{\lp {{{D}_{\mathcal{D}}}B} \rp }C} \rp } = {{\displaystyle}\sum_{l\in{\left \{
   {n} \rbrc}}{{\eb}}^{l}{\lp {A{\lp {D_{l}B} \rp }C} \rp }} \ee

or

.. math::

   \be {\lp {A\dot{B}C} \rp }\dot{D} = {\lp {A{\lp {{{D}_{\mathcal{D}}}B} \rp }C} \rp }{{D}_{\mathcal{G}}} = {{\displaystyle}\sum_{l\in{\left \{
   {n} \rbrc}}{\lp {A{\lp {D_{l}B} \rp }C} \rp }{{\eb}}^{l}}. \ee

The implementation of equations :math:`\ref{splitopV}` and :math:`\ref{splitopM}` is described in sections :ref:`makeMV` and :ref:`makeMVD`.

.. _Ltrans:

Linear Transformations/Outermorphisms
-------------------------------------

In the tangent space of a manifold, :math:`\mathcal{M}`, (which is a vector space) a linear transformation is the mapping :math:`\underline{T}\colon{{{\mathcal{T}_{\vec{x}}}\lp {\mathcal{M}} \rp }}\rightarrow{{{\mathcal{T}_{\vec{x}}}\lp {\mathcal{M}} \rp }}` (we use an underline to indicate a linear transformation) where for all :math:`x,y\in {{{\mathcal{T}_{\vec{x}}}\lp {\mathcal{M}} \rp }}` and :math:`\alpha\in\Re` we have

.. math::

   \begin{aligned}
       {{\underline{T}}\lp {x+y} \rp } =& {{\underline{T}}\lp {x} \rp } + {{\underline{T}}\lp {y} \rp } \\
       {{\underline{T}}\lp {\alpha x} \rp } =& \alpha{{\underline{T}}\lp {x} \rp }\end{aligned}

The outermorphism induced by :math:`\underline{T}` is defined for :math:`x_{1},\dots,x_{r}\in{{{\mathcal{T}_{\vec{x}}}\lp {\mathcal{M}} \rp }}` where :math:`\newcommand{\f}[2]{{#1}\lp {#2} \rp } \newcommand{\Tn}[2]{\f{\mathcal{T}_{#2}}{#1}} r\le\f{\dim}{\Tn{\mathcal{M}}{\vec{x}}}`

.. math::

   \be \newcommand{\f}[2]{{#1}\lp {#2} \rp }
   \newcommand{\W}{\wedge}
   \f{\underline{T}}{x_{1}\W\dots\W x_{r}} \equiv \f{\underline{T}}{x_{1}}\W\dots\W\f{\underline{T}}{x_{r}} \ee

If :math:`I` is the pseudo scalar for :math:`{{{\mathcal{T}_{\vec{x}}}\lp {\mathcal{M}} \rp }}` we also have the following definitions for determinate, trace, and adjoint (:math:`\overline{T}`) of :math:`\underline{T}`

.. math::

   \begin{align}
       \f{\underline{T}}{I} \equiv&\; \f{\det}{\underline{T}}I\text{,} \label{eq_82}\\
       \f{\tr}{\underline{T}} \equiv&\; \nabla_{y}\cdot\f{\underline{T}}{y}\text{,} \label{eq_83}\\ 
       x\cdot \f{\overline{T}}{y} \equiv&\; y\cdot \f{\underline{T}}{x}.\ \label{eq_84}\\
   \end{align}

If :math:`{\left \{{{{\eb}}_{i}} \rbrc}` is a basis for :math:`{{{\mathcal{T}_{\vec{x}}}\lp {\mathcal{M}} \rp }}` then we can represent :math:`\underline{T}` with the matrix :math:`\underline{T}_{i}^{j}` used as follows (Einstein summation convention as usual) -

.. math:: \be     \f{\underline{T}}{\eb_{i}} = \underline{T}_{i}^{j}\eb_{j}, \label{eq_85} \ee

The let :math:`{\lp {\underline{T}^{-1}} \rp }_{m}^{n}` be the inverse matrix of :math:`\underline{T}_{i}^{j}` so that :math:`{\lp {\underline{T}^{-1}} \rp }_{m}^{k}\underline{T}_{k}^{j} = \delta^{j}_{m}` and

.. math:: \be \underline{T}^{-1}{\lp {a^{i}{{\eb}}_{i}} \rp } = a^{i}{\lp {\underline{T}^{-1}} \rp }_{i}^{j}{{\eb}}_{j} \label{eq_85a} \ee

and calculate

.. math::

   \begin{aligned}
       \underline{T}^{-1}{\lp {\underline{T}{\lp {a} \rp }} \rp } &= \underline{T}^{-1}{\lp {\underline{T}{\lp {a^{i}{{\eb}}_{i}} \rp }} \rp } \nonumber \\
           &= \underline{T}^{-1}{\lp {a^{i}\underline{T}_{i}^{j}{{\eb}}_{j}} \rp } \nonumber \\
           &= a^{i}{\lp {\underline{T}^{-1}} \rp }_{i}^{j} \underline{T}_{j}^{k}{{\eb}}_{k} \nonumber \\
           &= a^{i}\delta_{i}^{j}{{\eb}}_{j} = a^{i}{{\eb}}_{i} = a.\end{aligned}

Thus if eq :math:`\ref{eq_85a}` is used to define the :math:`\underline{T}_{i}^{j}` then the linear transformation defined by the matrix :math:`{\lp {\underline{T}^{-1}} \rp }_{m}^{n}` is the inverse of :math:`\underline{T}`.

In eq. (:math:`\ref{eq_85}`) the matrix, :math:`\underline{T}_{i}^{j}`, only has it’s usual meaning if the :math:`{\left \{{{{\eb}}_{i}} \rbrc}` form an orthonormal Euclidean basis (Minkowski spaces not allowed). Equations (:math:`\ref{eq_82}`) through (:math:`\ref{eq_84}`) become

.. math::

   \begin{aligned}
       {{\det}\lp {\underline{T}} \rp } =&\; {{\underline{T}}\lp {{{\eb}}_{1}{\wedge}\dots{\wedge}{{\eb}}_{n}} \rp }{\lp {{{\eb}}_{1}{\wedge}\dots{\wedge}{{\eb}}_{n}} \rp }^{-1},\\
       {{{\mbox{tr}}}\lp {\underline{T}} \rp } =&\; \underline{T}_{i}^{i},\\
       \overline{T}_{j}^{i} =&\;  g^{il}g_{jp}\underline{T}_{l}^{p}.\end{aligned}

A important form of linear transformation with a simple representation is the spinor transformation. If :math:`S` is an even multivector we have :math:`SS^{{\dagger}} = \rho^{2}`, where :math:`\rho^{2}` is a scalar. Then :math:`S` is a spinor transformation is given by (:math:`v` is a vector)

.. math:: \be {{S}\lp {v} \rp } = SvS^{{\dagger}} \ee

if :math:`{{S}\lp {v} \rp }` is a vector and

.. math:: \be {{S^{-1}}\lp {v} \rp } = \frac{S^{{\dagger}}vS}{\rho^{4}}. \ee

Thus

.. math::

   \begin{aligned}
       {{S^{-1}}\lp {{{S}\lp {v} \rp }} \rp } &= \frac{S^{{\dagger}}SvS^{{\dagger}}S}{\rho^{4}} \nonumber \\
                            &= \frac{\rho^{2}v\rho^{2}}{\rho^{4}} \nonumber \\
                            &= v. \end{aligned}

One more topic to consider is whether or not :math:`T^{i}_{j}` should be called the matrix representation of :math:`T` ? The reason that this is a question is that for a general metric :math:`g_{ij}` is that because of the dependence of the dot product on the metric :math:`T^{i}_{j}` does not necessarily show the symmetries of the underlying transformation :math:`T`. Consider the expression

.. math::

   \begin{aligned}
       a\cdot{{T}\lp {b} \rp } &= a^{i}{{\eb}}_{i}\cdot{{T}\lp {b^{j}{{\eb}}_{j}} \rp } \nonumber \\
                      &= a^{i}{{\eb}}_{i}\cdot {{T}\lp {{{\eb}}_{j}} \rp }b^{j} \nonumber \\
                      &= a^{i}{{\eb}}_{i}\cdot{{\eb}}_{k} T_{j}^{k}b^{j} \nonumber \\
                      &= a^{i}g_{ik}T_{j}^{k}b^{j}.\end{aligned}

It is

.. math:: \be T_{ij} = g_{ik}T_{j}^{k} \ee

that has the proper symmetry for self adjoint transformations :math:`(a\cdot{{T}\lp {b} \rp } = b\cdot{{T}\lp {a} \rp })` in the sense that if :math:`T = \overline{T}` then :math:`T_{ij} = T_{ji}`. Of course if we are dealing with a manifold where the :math:`g_{ij}`\ ’s are functions of the coordinates then the matrix representation of a linear transformation will also be a function of the coordinates. Assuming we use :math:`T_{ij}` for the matrix representation of the linear transformation,
:math:`T`, then if we given the matrix representation, :math:`T_{ij}`, we can construct the linear transformation given by :math:`T^{i}_{j}` as follows

.. math::

   \begin{aligned}
       T_{ij} &= g_{ik}T_{j}^{k} \nonumber \\
       g^{li}T_{ij} &= g^{li}g_{ik}T_{j}^{k} \nonumber \\
       g^{li}T_{ij} &= \delta_{k}^{l}T_{j}^{k} \nonumber \\
       g^{li}T_{ij} &= T_{j}^{l}.\end{aligned}

Any program/code that represents :math:`T` should allow one to define :math:`T` in terms of :math:`T_{ij}` or :math:`T_{j}^{l}` and likewise given a linear transformation :math:`T` obtain both :math:`T_{ij}` and :math:`T_{j}^{l}` from it. Please note that these considerations come into play for any non-Euclidean metric with respect to the trace and adjoint of a linear transformation since calculating either requires a dot product.

.. _MLtrans:

Multilinear Functions
---------------------

A multivector multilinear function\ [10]_ is a multivector function :math:`{{T}\lp {A_{1},\dots,A_{r}} \rp }` that is linear in each of it arguments\ [11]_ (it could be implicitly non-linearly dependent on a set of additional arguments such as the position coordinates, but we only consider the linear arguments). :math:`T` is a *tensor* of degree :math:`r` if each variable :math:`A_{j}` is restricted to the vector space :math:`\mathcal{V}_{n}`. More generally if each
:math:`A_{j}\in{{\mathcal{G}}\lp {\mathcal{V}_{n}} \rp }` (the geometric algebra of :math:`\mathcal{V}_{n}`), we call :math:`T` an *extensor* of degree-:math:`r` on :math:`{{\mathcal{G}}\lp {\mathcal{V}_{n}} \rp }`.

If the values of :math:`{{T} \lp {a_{1},\dots,a_{r}} \rp }` :math:`\lp a_{j}\in\mathcal{V}_{n}\;\forall\; 1\le j \le r \rp` are :math:`s`-vectors (pure grade :math:`s` multivectors in :math:`{{\mathcal{G}}\lp {\mathcal{V}_{n}} \rp }`) we say that :math:`T` has grade :math:`s` and rank :math:`r+s`. A tensor of grade zero is called a *multilinear form*.

In the normal definition of tensors as multilinear functions the tensor is defined as a mapping

.. math:: T:{\huge \times}_{i=1}^{r}\mathcal{V}_{i}\rightarrow\Re,

\ so that the standard tensor definition is an example of a grade zero degree/rank$ r $ tensor in our definition.

Algebraic Operations
~~~~~~~~~~~~~~~~~~~~

The properties of tensors are (:math:`\alpha\in\Re`, :math:`a_{j},b\in\mathcal{V}_{n}`, :math:`T` and :math:`S` are tensors of rank :math:`r`, and :math:`\circ` is any multivector multiplicative operation)

.. math::

   \begin{aligned}
       {{T}\lp {a_{1},\dots,\alpha a_{j},\dots,a_{r}} \rp } =& \alpha{{T}\lp {a_{1},\dots,a_{j},\dots,a_{r}} \rp }, \\
       {{T}\lp {a_{1},\dots,a_{j}+b,\dots,a_{r}} \rp } =& {{T}\lp {a_{1},\dots,a_{j},\dots,a_{r}} \rp }+ {{T}\lp {a_{1},\dots,a_{j-1},b,a_{j+1},\dots,a_{r}} \rp }, \\
       {{\lp T\pm S\rp }\lp {a_{1},\dots,a_{r}} \rp } \equiv& {{T}\lp {a_{1},\dots,a_{r}} \rp }\pm{{S}\lp {a_{1},\dots,a_{r}} \rp }.\end{aligned}

Now let :math:`T` be of rank :math:`r` and :math:`S` of rank :math:`s` then the product of the two tensors is

.. math:: \be \f{\lp T\circ S\rp}{a_{1},\dots,a_{r+s}} \equiv \f{T}{a_{1},\dots,a_{r}}\circ\f{S}{a_{r+1},\dots,a_{r+s}}, \ee

where “:math:`\circ`” is any multivector multiplicative operation.

Covariant, Contravariant, and Mixed Representations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The arguments (vectors) of the multilinear function can be represented in terms of the basis vectors or the reciprocal basis vectors

.. math::

   \begin{aligned}
       a_{j} =& a^{i_{j}}{{\eb}}_{i_{j}}, \label{vrep}\\
             =& a_{i_{j}}{{\eb}}^{i_{j}}. \label{rvrep}\end{aligned}

Equation (:math:`\ref{vrep}`) gives :math:`a_{j}` in terms of the basis vectors and eq (:math:`\ref{rvrep}`) in terms of the reciprocal basis vectors. The index :math:`j` refers to the argument slot and the indices :math:`i_{j}` the components of the vector in terms of the basis. The covariant representation of the tensor is defined by

:math:`\newcommand{\indices}[1]{#1}\begin{aligned}  T\indices{_{i_{1}\dots i_{r}}} \equiv& {{T}\lp {{{\eb}}_{i_{1}},\dots,{{\eb}}_{i_{r}}} \rp } \\  {{T}\lp {a_{1},\dots,a_{r}} \rp } =& {{T}\lp {a^{i_{1}}{{\eb}}_{i_{1}},\dots,a^{i_{r}}{{\eb}}_{i_{r}}} \rp } \nonumber \\  =& {{T}\lp {{{\eb}}_{i_{1}},\dots,{{\eb}}_{i_{r}}} \rp }a^{i_{1}}\dots a^{i_{r}} \nonumber \\  =& T\indices{_{i_{1}\dots i_{r}}}a^{i_{1}}\dots a^{i_{r}}.\end{aligned}`\ $

Likewise for the contravariant representation

.. math::

   \begin{aligned}
   T\indices{^{i_{1}\dots i_{r}}} \equiv& {{T}\lp {{{\eb}}^{i_{1}},\dots,{{\eb}}^{i_{r}}} \rp } \\
       {{T}\lp {a_{1},\dots,a_{r}} \rp } =& {{T}\lp {a_{i_{1}}{{\eb}}^{i_{1}},\dots,a_{i_{r}}{{\eb}}^{i_{r}}} \rp } \nonumber \\
                                =& {{T}\lp {{{\eb}}^{i_{1}},\dots,{{\eb}}^{i_{r}}} \rp }a_{i_{1}}\dots a_{i_{r}} \nonumber \\
                                =& T\indices{^{i_{1}\dots i_{r}}}a_{i_{1}}\dots a_{i_{r}}.\end{aligned}

One could also have a mixed representation

.. math::

   \begin{aligned}
   T\indices{_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{r}}} \equiv& {{T}\lp {{{\eb}}_{i_{1}},\dots,{{\eb}}_{i_{s}},{{\eb}}^{i_{s+1}}\dots{{\eb}}^{i_{r}}} \rp } \\
       {{T}\lp {a_{1},\dots,a_{r}} \rp } =& {{T}\lp {a^{i_{1}}{{\eb}}_{i_{1}},\dots,a^{i_{s}}{{\eb}}_{i_{s}},
                                   a_{i_{s+1}}{{\eb}}^{i_{s}}\dots,a_{i_{r}}{{\eb}}^{i_{r}}} \rp } \nonumber \\
                                =& {{T}\lp {{{\eb}}_{i_{1}},\dots,{{\eb}}_{i_{s}},{{\eb}}^{i_{s+1}},\dots,{{\eb}}^{i_{r}}} \rp }
                                   a^{i_{1}}\dots a^{i_{s}}a_{i_{s+1}},\dots a^{i_{r}} \nonumber \\
                                =& T\indices{_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{r}}}a^{i_{1}}\dots a^{i_{s}}a_{i_{s+1}}\dots a^{i_{r}}.\end{aligned}

In the representation of :math:`T` one could have any combination of covariant (lower) and contravariant (upper) indexes.

To convert a covariant index to a contravariant index simply consider

.. math::

   \begin{aligned}
       \f{T}{\eb_{i_{1}},\dots,\eb^{i_{j}},\dots,\eb_{i_{r}}} =& \f{T}{\eb_{i_{1}},\dots,g^{i_{j}k_{j}}\eb_{k_{j}},\dots,\eb_{i_{r}}} \nonumber \\
                                                              =& g^{i_{j}k_{j}}\f{T}{\eb_{i_{1}},\dots,\eb_{k_{j}},\dots,\eb_{i_{r}}} \nonumber \\
       T_{i_{1}\dots}{}^{i_{j}}{}_{\dots i_{r}} =& g^{i_{j}k_{j}}T\indices{_{i_{1}\dots i_{j}\dots i_{r}}}.
   \end{aligned}

Similarly one could lower an upper index with :math:`g_{i_{j}k_{j}}`.

Contraction and Differentiation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The contraction of a tensor between the :math:`j^{th}` and :math:`k^{th}` variables (slots) is

.. math:: \be \f{T}{a_{i},\dots,a_{j-1},\nabla_{a_{k}},a_{j+1},\dots,a_{r}} = \nabla_{a_{j}}\cdot\lp \nabla_{a_{k}}\f{T}{a_{1},\dots,a_{r}}\rp . \ee

This operation reduces the rank of the tensor by two. This definition gives the standard results for *metric contraction* which is proved as follows for a rank :math:`r` grade zero tensor (the circumflex “:math:`\breve{\:\:}`” indicates that a term is to be deleted from the product).

.. math::

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
                                                 \breve{a}^{i_{j}}\dots\breve{a}^{i_{m}}\dots a^{i_{r}} \label{eq108}
   \end{align}

Equation (:math:`\ref{eq108}`) is the correct formula for the metric contraction of a tensor.

If we have a mixed representation of a tensor, :math:`T\indices{_{i_{1}\dots}{}^{i_{j}}{}_{\dots i_{k}\dots i_{r}}}`, and wish to contract between an upper and lower index (:math:`i_{j}` and :math:`i_{k}`) first lower the upper index and then use eq (:math:`\ref{eq108}`) to contract the result. Remember lowering the index does *not* change the tensor, only the *representation* of the tensor, while contraction results in a *new* tensor. First lower index

.. math:: \be T\indices{_{i_{1}\dots}{}^{i_{j}}{}_{\dots i_{k}\dots i_{r}}} \xRightarrow{\small Lower Index} g_{i_{j}k_{j}}T\indices{_{i_{1}\dots}{}^{k_{j}}{}_{\dots i_{k}\dots i_{r}}} \ee

Now contract between :math:`i_{j}` and :math:`i_{k}` and use the properties of the metric tensor.

.. math::

   \begin{aligned}
       g_{i_{j}k_{j}}T\indices{_{i_{1}\dots}{}^{k_{j}}{}_{\dots i_{k}\dots i_{r}}} \xRightarrow{\small Contract}&
                   g^{i_{j}i_{k}}g_{i_{j}k_{j}}T\indices{_{i_{1}\dots}{}^{k_{j}}{}_{\dots i_{k}\dots i_{r}}} \nonumber \\
                   =& \delta_{k_{j}}^{i_{k}}T\indices{_{i_{1}\dots}{}^{k_{j}}{}_{\dots i_{k}\dots i_{r}}}. \label{114a}\end{aligned}

Equation (:math:`\ref{114a}`) is the standard formula for contraction between upper and lower indexes of a mixed tensor.

Finally if :math:`{{T}\lp {a_{1},\dots,a_{r}} \rp }` is a tensor field (implicitly a function of position) the tensor derivative is defined as

.. math::

   \begin{aligned}
       {{T}\lp {a_{1},\dots,a_{r};a_{r+1}} \rp } \equiv \lp a_{r+1}\cdot\nabla\rp {{T}\lp {a_{1},\dots,a_{r}} \rp },\end{aligned}

assuming the :math:`a^{i_{j}}` coefficients are not a function of the coordinates.

This gives for a grade zero rank :math:`r` tensor

.. math::

   \begin{aligned}
       \lp a_{r+1}\cdot\nabla\rp {{T}\lp {a_{1},\dots,a_{r}} \rp } =& a^{i_{r+1}}\partial_{x^{i_{r+1}}}a^{i_{1}}\dots a^{i_{r}}
                                                           T_{i_{1}\dots i_{r}}, \nonumber \\
                                                        =& a^{i_{1}}\dots a^{i_{r}}a^{i_{r+1}}
                                                           \partial_{x^{i_{r+1}}}T_{i_{1}\dots i_{r}}.\end{aligned}

From Vector to Tensor
~~~~~~~~~~~~~~~~~~~~~

A rank one tensor is a vector since it satisfies all the axioms for a vector space, but a vector in not necessarily a tensor since not all vectors are multilinear (actually in the case of vectors a linear function) functions. However, there is a simple isomorphism between vectors and rank one tensors defined by the mapping :math:`{{v}\lp {a} \rp }:\mathcal{V}\rightarrow\Re` such that if :math:`v,a \in\mathcal{V}`

.. math:: \be \f{v}{a} \equiv v\cdot a. \ee

So that if :math:`v = v^{i}{{\eb}}_{i} = v_{i}{{\eb}}^{i}` the covariant and contravariant representations of :math:`v` are (using :math:`{{\eb}}^{i}\cdot{{\eb}}_{j} = \delta^{i}_{j}`)

.. math:: \be \f{v}{a} = v_{i}a^{i} = v^{i}a_{i}. \ee

Parallel Transport and Covariant Derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The covariant derivative of a tensor field :math:`{{T}\lp {a_{1},\dots,a_{r};x} \rp }` (:math:`x` is the coordinate vector of which :math:`T` can be a non-linear function) in the direction :math:`a_{r+1}` is (remember :math:`a_{j} = a_{j}^{k}{{\eb}}_{k}` and the :math:`{{\eb}}_{k}` can be functions of :math:`x`) the directional derivative of :math:`{{T}\lp {a_{1},\dots,a_{r};x} \rp }` where all the arguments of :math:`T` are parallel transported. The definition of parallel transport is if
:math:`a` and :math:`b` are tangent vectors in the tangent spaced of the manifold then

.. math:: \be     \paren{a\cdot\nabla_{x}}b = 0 \label{eq108a} \ee

if :math:`b` is parallel transported. Since :math:`b = b^{i}{{\eb}}_{i}` and the derivatives of :math:`{{\eb}}_{i}` are functions of the :math:`x^{i}`\ ’s then the :math:`b^{i}`\ ’s are also functions of the :math:`x^{i}`\ ’s so that in order for eq (:math:`\ref{eq108a}`) to be satisfied we have

.. math::

   \begin{aligned}
       {\lp {a\cdot\nabla_{x}} \rp }b =& a^{i}\partial_{x^{i}}{\lp {b^{j}{{\eb}}_{j}} \rp } \nonumber \\
                                 =& a^{i}{\lp {{\lp {\partial_{x^{i}}b^{j}} \rp }{{\eb}}_{j} + b^{j}\partial_{x^{i}}{{\eb}}_{j}} \rp } \nonumber \\
                                 =& a^{i}{\lp {{\lp {\partial_{x^{i}}b^{j}} \rp }{{\eb}}_{j} + b^{j}\Gamma_{ij}^{k}{{\eb}}_{k}} \rp } \nonumber \\
                                 =& a^{i}{\lp {{\lp {\partial_{x^{i}}b^{j}} \rp }{{\eb}}_{j} + b^{k}\Gamma_{ik}^{j}{{\eb}}_{j}} \rp }\nonumber \\
                                 =& a^{i}{\lp {{\lp {\partial_{x^{i}}b^{j}} \rp } + b^{k}\Gamma_{ik}^{j}} \rp }{{\eb}}_{j} = 0.\end{aligned}

Thus for :math:`b` to be parallel transported we must have

.. math:: \be     \partial_{x^{i}}b^{j} = -b^{k}\Gamma_{ik}^{j}. \label{eq121a} \ee

The geometric meaning of parallel transport is that for an infinitesimal rotation and dilation of the basis vectors (cause by infinitesimal changes in the :math:`x^{i}`\ ’s) the direction and magnitude of the vector :math:`b` does not change.

If we apply eq (:math:`\ref{eq121a}`) along a parametric curve defined by :math:`{{x^{j}}\lp {s} \rp }` we have

.. math::

   \begin{align}
       \deriv{b^{j}}{s}{} =& \deriv{x^{i}}{s}{}\pdiff{b^{j}}{x^{i}} \nonumber \\
                          =& -b^{k}\deriv{x^{i}}{s}{}\Gamma_{ik}^{j}, \label{eq122a}
   \end{align}

and if we define the initial conditions :math:`{{b^{j}}\lp {0} \rp }{{\eb}}_{j}`. Then eq (:math:`\ref{eq122a}`) is a system of first order linear differential equations with initial conditions and the solution, :math:`{{b^{j}}\lp {s} \rp }{{\eb}}_{j}`, is the parallel transport of the vector :math:`{{b^{j}}\lp {0} \rp }{{\eb}}_{j}`.

An equivalent formulation for the parallel transport equation is to let :math:`{{\gamma}\lp {s} \rp }` be a parametric curve in the manifold defined by the tuple :math:`{{\gamma}\lp {s} \rp } = {\lp {{{x^{1}}\lp {s} \rp },\dots,{{x^{n}}\lp {s} \rp }} \rp }`. Then the tangent to :math:`{{\gamma}\lp {s} \rp }` is given by

.. math:: \be \deriv{\gamma}{s}{} \equiv \deriv{x^{i}}{s}{}\eb_{i} \ee

and if :math:`{{v}\lp {x} \rp }` is a vector field on the manifold then

.. math::

   \begin{align}
       \paren{\deriv{\gamma}{s}{}\cdot\nabla_{x}}v =& \deriv{x^{i}}{s}{}\pdiff{}{x^{i}}\paren{v^{j}\eb_{j}} \nonumber \\
            =&\deriv{x^{i}}{s}{}\paren{\pdiff{v^{j}}{x^{i}}\eb_{j}+v^{j}\pdiff{\eb_{j}}{x^{i}}} \nonumber \\
            =&\deriv{x^{i}}{s}{}\paren{\pdiff{v^{j}}{x^{i}}\eb_{j}+v^{j}\Gamma^{k}_{ij}\eb_{k}} \nonumber \\
            =&\deriv{x^{i}}{s}{}\pdiff{v^{j}}{x^{i}}\eb_{j}+\deriv{x^{i}}{s}{}v^{k}\Gamma^{j}_{ik}\eb_{j} \nonumber \\
            =&\paren{\deriv{v^{j}}{s}{}+\deriv{x^{i}}{s}{}v^{k}\Gamma^{j}_{ik}}\eb_{j} \nonumber \\
            =& 0. \label{eq124a}
   \end{align}

Thus eq (:math:`\ref{eq124a}`) is equivalent to eq (:math:`\ref{eq122a}`) and parallel transport of a vector field along a curve is equivalent to the directional derivative of the vector field in the direction of the tangent to the curve being zero.

If the tensor component representation is contra-variant (superscripts instead of subscripts) we must use the covariant component representation of the vector arguments of the tensor, :math:`a = a_{i}{{\eb}}^{i}`. Then the definition of parallel transport gives

.. math::

   \begin{aligned}
       {\lp {a\cdot\nabla_{x}} \rp }b =& a^{i}\partial_{x^{i}}{\lp {b_{j}{{\eb}}^{j}} \rp } \nonumber \\
                                 =& a^{i}{\lp {{\lp {\partial_{x^{i}}b_{j}} \rp }{{\eb}}^{j} + b_{j}\partial_{x^{i}}{{\eb}}^{j}} \rp },\end{aligned}

and we need

.. math:: \be     \paren{\partial_{x^{i}}b_{j}}\eb^{j} + b_{j}\partial_{x^{i}}\eb^{j} = 0. \label{eq111a} \ee

To satisfy equation (:math:`\ref{eq111a}`) consider the following

.. math::

   \begin{aligned}
       \partial_{x^{i}}{\lp {{{\eb}}^{j}\cdot{{\eb}}_{k}} \rp } =& 0 \nonumber \\
       {\lp {\partial_{x^{i}}{{\eb}}^{j}} \rp }\cdot{{\eb}}_{k} + {{\eb}}^{j}\cdot{\lp {\partial_{x^{i}}{{\eb}}_{k}} \rp } =& 0  \nonumber \\
       {\lp {\partial_{x^{i}}{{\eb}}^{j}} \rp }\cdot{{\eb}}_{k} + {{\eb}}^{j}\cdot{{\eb}}_{l}\Gamma_{ik}^{l} =& 0 \nonumber \\
       {\lp {\partial_{x^{i}}{{\eb}}^{j}} \rp }\cdot{{\eb}}_{k} + \delta_{l}^{j}\Gamma_{ik}^{l} =& 0 \nonumber \\
       {\lp {\partial_{x^{i}}{{\eb}}^{j}} \rp }\cdot{{\eb}}_{k} + \Gamma_{ik}^{j} =& 0 \nonumber \\
       {\lp {\partial_{x^{i}}{{\eb}}^{j}} \rp }\cdot{{\eb}}_{k} =& -\Gamma_{ik}^{j}\end{aligned}

Now dot eq (:math:`\ref{eq111a}`) into :math:`{{\eb}}_{k}` giving

.. math::

   \begin{aligned}
       {\lp {\partial_{x^{i}}b_{j}} \rp }{{\eb}}^{j}\cdot{{\eb}}_{k} + b_{j}{\lp {\partial_{x^{i}}{{\eb}}^{j}} \rp }\cdot{{\eb}}_{k} =& 0  \nonumber \\
       {\lp {\partial_{x^{i}}b_{j}} \rp }\delta_{j}^{k} - b_{j}\Gamma_{ik}^{j} =& 0 \nonumber \\
       {\lp {\partial_{x^{i}}b_{k}} \rp } = b_{j}\Gamma_{ik}^{j}.\end{aligned}

Thus if we have a mixed representation of a tensor

.. math::

   \be \f{T}{a_{1},\dots,a_{r};x} =
       \f{T\indices{_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{r}}}}{x}a^{i_{1}}\dots a^{i_{s}}a_{i_{s+1}}\dots a_{i_{r}}, \ee

the covariant derivative of the tensor is

.. math::

   \begin{align}
       {\lp {a_{r+1}\cdot D} \rp } {{T}\lp {a_{1},\dots,a_{r};x} \rp } =&
           {{\displaystyle\frac{\partial {T\indices{_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{r}}}}}{\partial {x^{r+1}}}}}a^{i_{1}}\dots a^{i_{s}}a_{i_{s+1}}\dots a^{r}_{i_{r}}
           a^{i_{r+1}} \nonumber \\
           &\hspace{-0.5in}+ \sum_{p=1}^{s}{{\displaystyle\frac{\partial {a^{i_{p}}}}{\partial {x^{i_{r+1}}}}}}T\indices{_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{r}}}a^{i_{1}}\dots
           \breve{a}^{i_{p}}\dots a^{i_{s}}a_{i_{s+1}}\dots a_{i_{r}}a^{i_{r+1}} \nonumber \\
           &\hspace{-0.5in}+ \sum_{q=s+1}^{r}{{\displaystyle\frac{\partial {a_{i_{p}}}}{\partial {x^{i_{r+1}}}}}}T\indices{_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{r}}}a^{i_{1}}\dots
           a^{i_{s}}a_{i_{s+1}}\dots\breve{a}_{i_{q}}\dots a_{i_{r}}a^{i_{r+1}} \nonumber \\
           =& {{\displaystyle\frac{\partial {T\indices{_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{r}}}}}{\partial {x^{r+1}}}}}a^{i_{1}}\dots a^{i_{s}}a_{i_{s+1}}\dots a^{r}_{i_{r}}
           a^{i_{r+1}} \nonumber \\
           &\hspace{-0.5in}- \sum_{p=1}^{s}\Gamma_{i_{r+1}l_{p}}^{i_{p}}T\indices{_{i_{1}\dots i_{p}\dots i_{s}}^{i_{s+1}
           \dots i_{r}}}a^{i_{1}}\dots
           a^{l_{p}}\dots a^{i_{s}}a_{i_{s+1}}\dots a_{i_{r}}a^{i_{r+1}} \nonumber \\
           &\hspace{-0.5in}+ \sum_{q=s+1}^{r}\Gamma_{i_{r+1}i_{q}}^{l_{q}}T\indices{_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{q}
           \dots i_{r}}}a^{i_{1}}\dots
           a^{i_{s}}a_{i_{s+1}}\dots a_{l_{q}}\dots a_{i_{r}}a^{i_{r+1}}   .   \label{eq126a} \\
   \end{align}

From eq (:math:`\ref{eq126a}`) we obtain the components of the covariant derivative to be

.. math::

   \begin{aligned}
       {{\displaystyle\frac{\partial {T\indices{_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{r}}}}}{\partial {x^{r+1}}}}}
       - \sum_{p=1}^{s}\Gamma_{i_{r+1}l_{p}}^{i_{p}}T\indices{_{i_{1}\dots i_{p}\dots i_{s}}^{i_{s+1}\dots i_{r}}}
       + \sum_{q=s+1}^{r}\Gamma_{i_{r+1}i_{q}}^{l_{q}}T\indices{_{i_{1}\dots i_{s}}^{i_{s+1}\dots i_{q}\dots i_{r}}}.\end{aligned}

The component free form of the covariant derivative (the one used to calculate it in the code) is

.. math::

   \be \mathcal{D}_{a_{r+1}} {{T}\lp {a_{1},\dots,a_{r};x} \rp } \equiv \nabla T
           - \sum_{k=1}^{r}{{T}\lp {a_{1},\dots,{\lp {a_{r+1}\cdot\nabla} \rp } a_{k},\dots,a_{r};x} \rp }. \ee

Representation of Multivectors in *sympy*
-----------------------------------------

The *sympy* python module offers a simple way of representing multivectors using linear combinations of commutative expressions (expressions consisting only of commuting *sympy* objects) and non-commutative symbols. We start by defining :math:`n` non-commutative *sympy* symbols as a basis for the vector space

``(e_1,...,e_n) = symbols('e_1,...,e_n',commutative=False,real=True)``

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

``def_prec(gd,op_ord='<>|,^,*')``

   Define the precedence of the multivector operations. The function ``def_prec()`` must be called from the main program and the first argument ``gd`` must be set to ``globals()``. The second argument ``op_ord`` determines the operator precedence for expressions input to the function ``GAeval()``. The default value of ``op_ord`` is ``'<>|,^,*'``. For the default value the ``<``, ``>``, and ``|`` operations have equal precedence followed by ``^``, and ``^`` is followed by ``*``.

``GAeval(s,pstr=False)``

   The function ``GAeval()`` returns a multivector expression defined by the string ``s`` where the operations in the string are parsed according to the precedences defined by ``def_prec()``. ``pstr`` is a flag to print the input and output of ``GAeval()`` for debugging purposes. ``GAeval()`` works by adding parenthesis to the input string ``s`` with the precedence defined by ``op_ord='<>|,,*'``. Then the parsed string is converted to a *sympy* expression using the python ``eval()`` function.
   For example consider where ``X``, ``Y``, ``Z``, and ``W`` are multivectors

   .. code:: python

      def_prec(globals())
      V = GAeval('X|Y^Z*W')

   The *sympy* variable ``V`` would evaluate to ``((X|Y)^Z)*W``.

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

Note that the strings shown in the above equation are only used when the values of :math:`g_{ij}` are output (printed). In the ga module (library) the :math:`g_{ij}` symbols are stored in a member of the geometric algebra instance so that if ``o3d`` is a geometric algebra then ``o3d.g`` is the metric tensor ( :math:`g_{ij} =` ``o3d.g[i,j]``) for that algebra.

The default definition of :math:`g` can be overwritten by specifying a string that will define :math:`g`. As an example consider a symbolic representation for conformal geometry. Define for a basis

.. code:: python

   basis = 'a0 a1 a2 n nbar'

and for a metric

.. code:: python

   g = '# # # 0 0, # # # 0 0, # # # 0 0, 0 0 0 0 2, 0 0 0 2 0'

then calling ``cf3d = Ga(basis,g=g)`` would initialize the metric tensor

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

   o3d = Ga('a_1 a_2 a_3',g=[1,1,1])
   (a_1,a_2,a_3) = o3d.mv()

or use the ``Ga.build()`` function -

.. code:: python

   (o3d,a_1,a_2,a_3) = Ga.build('a_1 a_2 a_3',g=[1,1,1])

Note that the python variable name for a basis vector does not have to correspond to the name give in ``Ga()`` or ``Ga.build()``, one may wish to use a shortened python variable name to reduce programming (typing) errors, for example one could use -

.. code:: python

   (o3d,a1,a2,a3) = Ga.build('a_1 a_2 a_3',g=[1,1,1])

or

.. code:: python

   (st4d,g0,g1,g2,g3) = Ga.build('gamma_0 gamma_1 gamma_2 gamma_3',\
                                   g=[1,-1,-1,-1])

for Minkowski space time.

If the latex printer is used ``e1`` would print as :math:`{\boldsymbol{e_{1}}}` and ``g1`` as :math:`{\boldsymbol{\gamma_{1}}}`.

Representation and Reduction of Multivector Bases
-------------------------------------------------

In our symbolic geometric algebra all multivectors can be obtained from the symbolic basis vectors we have input, via the different operations available to geometric algebra. The first problem we have is representing the general multivector in terms terms of the basis vectors. To do this we form the ordered geometric products of the basis vectors and develop an internal representation of these products in terms of python classes. The ordered geometric products are all multivectors of the form
:math:`a_{i_{1}}a_{i_{2}}\dots a_{i_{r}}` where :math:`i_{1}<i_{2}<\dots <i_{r}` and :math:`r \le n`. We call these multivectors bases and represent them internally with non-commutative symbols so for example :math:`a_{1}a_{2}a_{3}` is represented by

.. code:: python

   Symbol('a_1*a_2*a_3',commutative=False)

In the simplest case of two basis vectors ``a_1`` and ``a_2`` we have a list of bases

.. code:: python

   self.bases = [[Symbol('a_1',commutative=False,real=True),\
                Symbol('a_2',commutative=False,real=True)],\
                [Symbol('a_1*a_2',commutative=False,real=True)]]

For the case of the basis blades we have

.. code:: python

   self.blades = [[Symbol('a_1',commutative=False,real=True),\
                 Symbol('a_2',commutative=False,real=True)],\
                 [Symbol('a_1^a_2',commutative=False,real=True)]]

The index tuples for the bases of each pseudo grade and each grade for the case of dimension 3 is

.. code:: python

   self.indexes = (((0,),(1,),(2,)),((0,1),(0,2),(1,2)),((0,1,2)))

Then the non-commutative symbol representing each base is constructed from each index tuple. For example for ``self.indexes[1][1]`` the symbol is ``Symbol('a_1*a_3',commutative=False)``.

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

Module Components
=================

The geometric algebra module consists of the following files and classes

+----------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------+
| File           | Classes            | Usage                                                                                                                                                |
+================+====================+======================================================================================================================================================+
| ``metric.py``  | ``Metric``         | Instantiates metric tensor and derivatives of basis vectors. Normalized basis if required.                                                           |
+----------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``ga.py``      | ``Ga``             | Instantiates geometric algebra (inherits ``Metric``), generates bases, blades, multiplication tables, reciprocal basis, and left and right geometric |
|                |                    | derivative operators.                                                                                                                                |
+----------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------+
|                | ``Sm``             | Instantiates geometric algebra for submainfold (inherits ``Ga``).                                                                                    |
+----------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``mv.py``      | ``Mv``             | Instantiates multivector.                                                                                                                            |
+----------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------+
|                | ``Dop``            | Instantiates linear multivector differential operator.                                                                                               |
+----------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``lt.py``      | ``Lt``             | Instantiates multivector linear transformation.                                                                                                      |
+----------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``printer.py`` | ``Eprint``         | Starts enhanced text printing on ANSI terminal (requires ``ConEmu`` on Windows).                                                                     |
+----------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------+
|                | ``GaPrinter``      | Text printer for all geometric algebra classes (inherits from ``sympy`` ``StringPrinter``).                                                          |
+----------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------+
|                | ``GaLatexPrinter`` | :math:`\LaTeX`\ printer for all geometric algebra classes (inherits from ``sympy`` ``LatexPrinter``).                                                |
+----------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------+

Instantiating a Geometric Algebra
---------------------------------

The geometric algebra class is instantiated with

``Ga(basis,g=None,coords=None,X=None,norm=False,sig='e',Isq='-',wedge=True,debug=False)``

   The ``basis`` and ``g`` parameters were described in section :ref:`BasisMetric`. If the metric is a function of position, if we have multivector fields, or we wish to calculate geometric derivatives a coordinate set, ``coords``, is required. ``coords`` is a list of *sympy* symbols. For the case of instantiating a 3-d geometric algebra in spherical coordinates we have

   .. code:: python

      (r, th, phi) = coords = symbols('r,theta,phi', real=True)
      basis = 'e_r e_theta e_phi'
      g = [1, r**2, r**2*sin(th)**2]
      sp3d = Ga(basis,g=g,coords=coords,norm=True)

   The input ``X`` allows the metric to be input as a vector manifold. ``X`` is a list of functions of ``coords`` of dimension, :math:`m`, equal to or greater than the number of coordinates. If ``g=None`` it is assumed that ``X`` is a vector in an :math:`m`-dimensional orthonormal Euclidean vector space. If it is wished the embedding vector space to be non-Euclidean that condition is specified with ``g``. For example if we wish the embedding space to be a 5-dimensional Minkowski space then
   ``g=[-1,1,1,1,1]``. Then the Ga class uses ``X`` to calculate the manifold basis vectors as a function of the coordinates and from them the metric tensor\ [12]_.

   If ``norm=True`` the basis vectors of the manifold are normalized so that the absolute values of the squares of the basis vectors are one. *Currently you should only use this option for diagonal metric tensors, and even there due so with caution, due to the possible problems with taking the square root of a general*\ sympy\* expression (one that has an unknown sign).\*

   **When a geometric algebra is created the unnormalized metric tensor is always saved so that submanifolds created from the normalized manifold can be calculated correctly.**

   ``sig`` indicates the signature of the vector space in the following ways\ [13]_.

   1. If the metric tensor is purely numerical (the components are not symbolic or functions of the coordinates) and is diagonal (orthogonal basis vectors) the signature is computed from the metric tensor.

   2. If the metric tensor is not purely numerical and orthogonal the following hints are used (dimension of vector space is :math:`n`)

      1. ``sig='e'`` the default hint assumes the signature is for a Euclidean space with signature :math:`(n,0)`.

      2. ``sig='m+'`` assumes the signature if for the Minkowski space :math:`(n-1,1)`.

      3. ``sig='m-'`` assumes the signature if for the Minkowski space :math:`(1,n-1)`.

      4. ``sig=p`` where ``p`` is an integer :math:`p\le n` and the signature it :math:`(p,n-p)`.

   If the metric tensor contains no symbolic constants, but is a function of the coordinates, it is possible to determine the signature of the metric numerically by specifying a allowed numerical coordinate tuple due to the invariance of the signature. This will be implemented in the future.

   Currently one need not be concerned about inputting ``sig`` unless one in using the *Ga* member function ``Ga.I()`` or the functions ``Mv.dual()`` or ``cross()`` which also use ``Ga.I()``.

   If :math:`I^{2}` is numeric it is calculated if it is not numeric then ``Isq='-'`` is the sign of the square of the pseudo-scalar. This is needed for some operations. The default is chosen for the case of a general 3D Euclidean metric.

   If ``wedge=True`` the basis blades of a multivector are printed using the ``^`` symbol between basis vectors. If ``wedge=False`` the subscripts of each individual basis vector (assuming that the basis vector symbols are of the form root symbol with a subscript\ [14]_). For example in three dimensions if the basis vectors are :math:`{{\eb}}_{x}`, :math:`{{\eb}}_{y}`, and :math:`{{\eb}}_{z}` the grade 3 basis blade would be printed as :math:`{{\eb}}_{xyz}`.

   If ``debug=True`` the data structures required to initialize the Ga class are printed out.

   To get the basis vectors for ``sp3d`` we would have to use the member function ``Ga.mv()`` in the form

   .. code:: python

      (er,eth,ephi) = sp3d.mv() 

To access the reciprocal basis vectors of the geometric algebra use the member function ``mvr()``

``Ga.mvr(norm='True')``

   ``Ga.mvr(norm)`` returns the reciprocal basis vectors as a tuple. This allows the programmer to attach any python variable names to the reciprocal basis vectors that is convenient. For example (demonstrating the use of both ``mv()`` and ``mvr()``)

   .. code:: python

      (e_x,e_y,e_z) = o3d.mv()
      (e__x,e__y,e__z) = o3d.mvr()

   If ``norm='True'`` or the basis vectors are orthogonal the dot product of the basis vector and the corresponding reciprocal basis vector is one :math:`{\lp {e_{i}\cdot e^{j}=\delta_{i}^{j}} \rp }`. If ``norm='False'`` and the basis is non-orthogonal The dot product of the basis vector and the corresponding reciprocal basis vector is the square of the pseudo scalar, :math:`I^{2}`, of the geometric algebra :math:`{\lp {e_{i}\cdot e^{j}=E^{2}\delta_{i}^{j}} \rp }`.

In addition to the basis vectors, if coordinates are defined for the geometric algebra, the left and right geometric derivative operators are calculated and accessed with the ``Ga`` member function ``grads()``.

``Ga.grads()``

   ``Ga.grads()`` returns a tuple with the left and right geometric derivative operators. A typical usage would be

   .. code:: python

      (grad,rgrad) = sp3d.grads()

   for the spherical 3-d geometric algebra. The left derivative :math:`{\lp {{\texttt{grad}} ={\boldsymbol{\nabla}}} \rp }` and the right derivative :math:`{\lp {{\texttt{rgrad}} = {\boldsymbol{\bar{\nabla}}}} \rp }` have been explained in section :ref:`ldops`. Again the names ``grad`` and ``rgrad`` used in a program are whatever the user chooses them to be. In the previous example ``grad`` and ``rgrad`` are used.

an alternative instantiation method is

``Ga.build(basis, g=None, coords=None, X=None, norm=False, debug=False)``

   The input parameters for ``Ga.build()`` are the same as for ``Ga()``. The difference is that in addition to returning the geometric algebra ``Ga.build()`` returns the basis vectors at the same time. Using ``Ga.build()`` in the previous example gives

..

   .. code:: python

      (r, th, phi) = coords = symbols('r,theta,phi', real=True)
      basis = 'e_r e_theta e_phi'
      g = [1, r**2, r**2*sin(th)**2]
      (sp3d,er,eth,ephi) = Ga.build(basis,g=g,coord=coords,norm=True)

To access the pseudo scalar of the geometric algebra use the member function ``I()``.

``Ga.I()``

   ``Ga.I()`` returns the normalized pseudo scalar :math:`{\lp {{\left |{I^{2}}\right |}=1} \rp }` for the geometric algebra. For example :math:`I = \mbox{{\texttt{o3d.I()}}}` for the ``o3d`` geometric algebra. This function requires the signature of the vector space (see instantiating a geometric algebra).

``Ga.E()``

   ``Ga.E()`` returns the unnormalized pseudo scalar :math:`E_{n} = {\eb}_{1}{\wedge}\dots{\wedge}{\eb}_{n}` for the geometric algebra.

In general we have defined member functions of the ``Ga`` class that will instantiate objects of other classes since the objects of the other classes are all associated with a particular geometric algebra object. Thus we have

===================== ======= =============
Object                Class   ``Ga`` method
===================== ======= =============
multivector           ``Mv``  ``mv``
submanifold           ``Sm``  ``sm``
linear transformation ``Lt``  ``lt``
differential operator ``Dop`` ``dop``
===================== ======= =============

for the instantiation of various objects from the ``Ga`` class. This means that in order to instantiate any of these objects we need only to import ``Ga`` into our program.

.. _makeMV:

Instantiating a Multivector
---------------------------

Since we need to associate each multivector with the geometric algebra that contains it we use a member function of Ga to instantiate every multivector\ [15]_ The multivector is instantiated with:

``Ga.mv(name, mode, f=False)``

   As an example of both instantiating a geometric algebra and multivectors consider the following code fragment for a 3-d Euclidean geometric algebra.

   .. code:: python

      from sympy import symbols
      from ga import Ga
      (x, y, z) = coords = symbols('x,y,z',real=True)
      o3d = Ga('e_x e_y e_z', g=[1,1,1], coords=coords)
      (ex, ey, ez) = o3d.mv()
      V = o3d.mv('V','vector',f=True)
      f = o3d.mv(x*y*z)
      B = o3d.mv('B',2)

   First consider the multivector instantiation in line 6,

   ``V = o3d.mv('V','vector',f=True)``

   .Here a 3-dimensional multivector field that is a function of ``x``, ``y``, and ``z`` (``f=True``) is being instantiated. If latex output were used (to be discussed later) the multivector ``V`` would be displayed as

   .. math:: \be V^{x}\eb_{x} + V^{y}\eb_{y} + V^{z}\eb_{z} \ee

   Where the coefficients of the basis vectors are generalized *sympy* functions of the coordinates. If ``f=(x,y)`` then the coefficients would be functions of ``x`` and ``y``. In general is ``f`` is a tuple of symbols then the coefficients of the basis would be functions of those symbols. The superscripts\ [16]_ are formed from the coordinate symbols or if there are no coordinates from the subscripts of the basis vectors. The types of name and modes available for multivector instantiation are

   ======== ========================== ===============================================
   ``name`` ``mode``                   result
   ======== ========================== ===============================================
   string s ``scalar``                 symbolic scalar of value Symbol(s)
   string s ``vector``                 symbolic vector
   string s ``grade2`` or ``bivector`` symbolic bivector
   string s ``r`` (integer)            symbolic r-grade multivector
   string s ``pseudo``                 symbolic pseudoscalar
   string s ``spinor``                 symbolic even multivector
   string s ``mv``                     symbolic general multivector
   scalar c None                       zero grade multivector with coefficient value c
   ======== ========================== ===============================================

   Line 5 of the previous listing illustrates the case of using the ``mv`` member function with no arguments. The code does not return a multivector, but rather a tuple or the basis vectors of the geometric algebra ``o3d``. The elements of the tuple then can be used to construct multivectors, or multivector fields through the operations of addition, subtraction, multiplication (geometric, inner, and outer products and left and right contraction). As an example we could construct the vector
   function

   .. code:: python

      F = x**2*ex + z*ey + x*y*ez

   or the bivector function

   .. code:: python

      B = z*(ex^ey) + y*(ey^ez) + y*(ex^ez).

   Line 7 is an example of instantiating a multivector scalar function (a multivector with only a scalar part). If we print ``f`` the result is ``x*y*z``. Line 8 is an example of instantiating a grade :math:`r` (in the example a grade 2) multivector where

   .. math:: \be B = B^{xy}{\eb}_{x}{\wedge}{\eb}_{y}+B^{yz}{\eb}_{y}{\wedge}{\eb}_{z}+B^{xz}{\eb}_{x}{\wedge}{\eb}_{z}. \ee

If one wished to calculate the left and right geometric derivatives of ``F`` and ``B`` the required code would be

.. code:: python

   (grad,rgrad) = o3d.grads()
   dF = grad*F
   dB = grad*B
   dFr = F*rgrad
   dBr = B*rgrad

``dF``, ``dB``, ``dFr``, and ``dBr`` are all multivector functions. For the code where the order of the operations are reversed

.. code:: python

   (grad,rgrad) = o3d.grads()
   dFop = F*grad
   dBop = B*grad
   dFrop = rgrad*F
   dBrop = rgrad*B

``dFop``, ``dBop``, ``dFrop``, and ``dBrop`` are all multivector differential operators (again see section :ref:`ldops`).

Backward Compatibility Class MV
-------------------------------

In order to be backward compatible with older versions of *galgebra* we introduce the class MV which is inherits it’s functions from then class Mv. To instantiate a geometric algebra using MV use the static function

``MV.setup(basis, metric=None, coords=None, rframe=False, debug=False, curv=(None,None))``

   This function allows a single geometric algebra to be created. If the function is called more than once the old geometric algebra is overwritten by the new geometric algebra. The named input ``metric`` is the same as the named input ``g`` in the current version of *galgebra*. Likewise, ``basis``, ``coords``, and ``debug`` are the same in the old and current versions of *galgebra*\ \ [17]_. Due to improvements in *sympy* the inputs ``rframe`` and ``curv[1]`` are no longer required. ``curv[0]`` is
   the vector function (list or tuple of scalar functions) of the coordinates required to define a vector manifold. For compatibility with the old version of *galgebra* if ``curv`` is used ``metric`` should be a orthonormal Euclidean metric of the same dimension as ``curv[0]``. It is strongly suggested that one use the new methods of defining a geometric algebra on a manifold.

``MV(base, mvtype, fct=False, blade_rep=True)``

   For the instantiation of multivector using ``MV`` the ``base`` and ``mvtype`` arguments are the same as for new methods of multivector instantiation. The ``fct`` input is the same and the ``g`` input in the new methods. ``blade_rep`` is not used in the new methods so setting ``blade_rep=False`` will do nothing. Effectively ``blade_rep=False`` was not used in the old examples.

``Fmt(self, fmt=1, title=None)``

   ``Fmt`` in ``MV`` has inputs identical to ``Fmt`` in ``Mv`` except that if ``A`` is a multivector then ``A.Fmt(2,'A')`` executes a print statement from ``MV`` and returns ``None``, while from ``Mv``, ``A.Fmt(2,'A')`` returns a string so that the function is compatible with use in *ipython notebook*.

Basic Multivector Class Functions
---------------------------------

If we can instantiate multivectors we can use all the multivector class functions as described as follows.

``blade_coefs(self,basis_lst)``

   Find coefficients (sympy expressions) of multivector basis blade expansion corresponding to basis blades in ``basis_lst``. For example if :math:`V = V^{x}{{\eb}}_{x}+V^{y}{{\eb}}_{x}+V^{z}{{\eb}}_{x}` Then :math:`V\text{.blade_coefs}([{{\eb}}_{z},{{\eb}}_{x}]) = [V^{z},V^{x}]` or if :math:`B = B^{xy}{{\eb}}_{x}{\wedge}{{\eb}}_{y}+V^{yz}{{\eb}}_{y}{\wedge}{{\eb}}_{z}` then :math:`B\text{.blade_coefs}([{{\eb}}_{x}{\wedge}{{\eb}}_{y}]) = [B^{xy}]`.

``convert_to_blades(self)``

   Convert multivector from the base representation to the blade representation. If multivector is already in blade representation nothing is done.

``convert_from_blades(self)``

   Convert multivector from the blade representation to the base representation. If multivector is already in base representation nothing is done.

``diff(self,var)``

   Calculate derivative of each multivector coefficient with respect to variable ``var`` and form new multivector from coefficients.

``dual(self)``

   The mode of the ``dual()`` function is set by the ``Ga`` class static member function, ``GA.dual_mode(mode='I+')`` of the ``GA`` geometric galgebra which sets the following return values (:math:`I` is the pseudo-scalar for the geometric algebra ``GA``)

   =========== ================
   ``mode``    Return Value
   =========== ================
   ``'+I'``    :math:`IA`
   ``'I+'``    :math:`AI`
   ``'-I'``    :math:`-IA`
   ``'I-'``    :math:`-AI`
   ``'+Iinv'`` :math:`I^{-1}A`
   ``'Iinv+'`` :math:`AI^{-1}`
   ``'-Iinv'`` :math:`-I^{-1}A`
   ``'Iinv-'`` :math:`-AI^{-1}`
   =========== ================

   For example if the geometric algebra is ``o3d``, ``A`` is a multivector in ``o3d``, and we wish to use ``mode='I-'``. We set the mode with the function ``o3d.dual('I-')`` and get the dual of ``A`` with the function ``A.dual()`` which returns :math:`-AI`.

   If ``o3d.dual(mode)`` is not called the default for the dual mode is ``mode='I+'`` and ``A*I`` is returned.

   Note that ``Ga.dual(mode)`` used the function ``Ga.I()`` to calculate the normalized pseudoscalar. Thus if the metric tensor is not numerical and orthogonal the correct hint for then ``sig`` input of the *Ga* constructor is required.

``even(self)``

   Return the even grade components of the multivector.

``exp(self,hint='-')``

   If :math:`A` is a multivector then :math:`e^{A}` is defined for any :math:`A` via the series expansion for :math:`e`. However as a practical matter we only have a simple closed form formula for :math:`e^{A}` if :math:`A^{2}` is a scalar\ [18]_. If :math:`A^{2}` is a scalar and we know the sign of :math:`A^{2}` we have the following formulas for :math:`e^{A}`.

   .. math::

      $\begin{aligned}
      A^{2} > 0 : & & &\\
      A &= \sqrt{A^{2}} {\displaystyle\frac{A}{\sqrt{A^{2}}}} ,& e^{A} &= {{\cosh}\lp {\sqrt{A^{2}}} \rp }+{{\sinh}\lp {\sqrt{A^{2}}} \rp }{\displaystyle\frac{A}{\sqrt{A^{2}}}} \\
      A^{2} < 0 : & & &\\
      A &= \sqrt{-A^{2}} {\displaystyle\frac{A}{\sqrt{-A^{2}}}} ,& e^{A} &= {{\cos}\lp {\sqrt{-A^{2}}} \rp }+{{\sin}\lp {\sqrt{-A^{2}}} \rp }{\displaystyle\frac{A}{\sqrt{-A^{2}}}} \\
      A^{2} = 0 : & & &\\
      A &=0 ,& e^{A} &= 1 + A 
      \end{aligned}

   The hint is required for symbolic multivectors :math:`A` since in general *sympy* cannot determine if :math:`A^{2}` is positive or negative. If :math:`A` is purely numeric the hint is ignored since the sign can be calculated.

``expand(self)``

   Return multivector in which each coefficient has been expanded using *sympy* ``expand()`` function.

``factor(self)``

   Apply the ``sympy`` ``factor`` function to each coefficient of the multivector.

``Fmt(self, fmt=1,title=None)``

   Fuction to print multivectors in different formats where

   ======= ============================================
   ``fmt`` 
   ======= ============================================
   1       Print entire multivector on one line.
   2       Print each grade of multivector on one line.
   3       Print each base of multivector on one line.
   ======= ============================================

   ``title`` appends a title string to the beginning of the output. An equal sign in the title string is not required, but is added as a default. Note that ``Fmt`` only overrides the the global multivector printing format for the particular instance being printed. To reset the global multivector printing format use the function ``Fmt()`` in the printer module.

``func(self,fct)``

   Apply the ``sympy`` scalar function ``fct`` to each coefficient of the multivector.

``grade(self,igrade=0)``

   Return a multivector that consists of the part of the multivector of grade equal to ``igrade``. If the multivector has no ``igrade`` part return a zero multivector.

``inv(self)``

   Return the inverse of the multivector :math:`M` (``M.inv()``). If :math:`M` is a non-zero scalar return :math:`1/M`. If :math:`M^{2}` is a non-zero scalar return :math:`M/{\lp {M^{2}} \rp }`, If :math:`MM^{{\dagger}}` is a non-zero scalar return :math:`M^{{\dagger}}/{\lp {MM^{{\dagger}}} \rp }`. Otherwise exit the program with an error message.

   All division operators (``/``, ``/=``) use right multiplication by the inverse.

``norm(self,hint='+')``

   Return the norm of the multivector :math:`M` (``M.norm()``) defined by :math:`\sqrt{{\left |{MM^{{\dagger}}}\right |}}`. If :math:`MM^{{\dagger}}` is a scalar (a *sympy* scalar is returned). If :math:`MM^{{\dagger}}` is not a scalar the program exits with an error message. If :math:`MM^{{\dagger}}` is a number *sympy* can determine if it is positive or negative and calculate the absolute value. If :math:`MM^{{\dagger}}` is a *sympy* expression (function) *sympy* cannot determine the sign of
   the expression so that ``hint='+'`` or ``hint='-'`` is needed to determine if the program should calculate :math:`\sqrt{MM^{{\dagger}}}` or :math:`\sqrt{-MM^{{\dagger}}}`. For example if we are in a Euclidean space and ``M`` is a vector then ``hint='+'``, if ``M`` is a bivector then let ``hint='-'``. If ``hint='0'`` and :math:`MM^{{\dagger}}` is a symbolic scalar ``sqrt(Abs(M*M.rev()))`` is returned where ``Abs()`` is the *sympy* symbolic absolute value function.

``norm2(self)``

   Return the the scalar defined by :math:`MM^{{\dagger}}` if :math:`MM^{{\dagger}}` is a scalar. If :math:`MM^{{\dagger}}` is not a scalar the program exits with an error message.

``proj(self,bases_lst)``

   Return the projection of the multivector :math:`M` (``M.proj(bases_lst)``) onto the subspace defined by the list of bases (``bases_lst``).

``proj(self,lst)``

   Return the projection of the mutivector :math:`A` onto the list, :math:`lst`, of basis blades. For example if :math:`A = A^{x}{{\eb}}_{x}+A^{y}{{\eb}}_{y}+A^{z}{{\eb}}_{z}` then :math:`A.proj{\lp {[{{\eb}}_{x},{{\eb}}_{y}]} \rp } = A^{x}{{\eb}}_{x}+A^{y}{{\eb}}_{y}`. Similarly if :math:`A = A^{xy}{{\eb}}_{x}{\wedge}{{\eb}}_{y}+A^{yz}{{\eb}}_{y}{\wedge}{{\eb}}_{z}` then :math:`A.proj{\lp {[{{\eb}}_{x}{\wedge}{{\eb}}_{y}]} \rp } = A^{xy}{{\eb}}_{x}{\wedge}{{\eb}}_{y}`.

``project_in_blade(self,blade)``

   Return the projection of the mutivector :math:`A` in subspace defined by the blade, :math:`B`, using the formula :math:`{\lp {A\rfloor B} \rp }B^{-1}` in :cite:`Macdonald1`, page 121.

``pure_grade(self)``

   If the multivector :math:`A` is pure (only contains one grade) return, :math:`A.pure\_grade()`, the index ('0' for a scalar, '1' for vector, '2' for a bi-vector, etc.) of the non-zero grade. If :math:`A` is not pure return the negative of the highest non-zero grade index.

``odd(self)``

   Return odd part of multivector.

``reflect_in_blade(self,blade)``

   Return the reflection of the mutivector :math:`A` in the subspace defined by the :math:`r`-grade blade, :math:`B_{r}`, using the formula (extended to multivectors) :math:`\sum_{i} {\lp {-1} \rp }^{r{\lp {i+1} \rp }}{B}_{r}{\left < {A} \right >}_{i}B_{r}^{-1}` in :cite:`Macdonald1`, page 129.

``rev(self)``

   Return the reverse of the multivector.

``rotate_multivector(self,itheta,hint='-')``

   Rotate the multivector :math:`A` via the operation :math:`e^{-\theta i/2}Ae^{\theta i/2}` where itheta = :math:`\theta i`, :math:`\theta` is a scalar, and :math:`i` is a unit, :math:`i^{2} = \pm 1`, 2-blade. If :math:`{\lp {\theta i} \rp }^{2}` is not a number ``hint`` is required to determine the sign of the square of ``itheta``. The default chosen, ``hint='-'``, is correct for any Euclidean space.

``scalar(self)``

   Return the coefficient (*sympy* scalar) of the scalar part of a multivector.

``simplify(self,mode=simplify)``

   ``mode`` is a *sympy* simplification function of a list/tuple of *sympy* simplification functions that are applied in sequence (if more than one function) each coefficient of the multivector. For example if we wished to applied ``trigsimp`` and ``ratsimp`` *sympy* functions to the multivector ``F`` the code would be

   .. code:: python

      Fsimp = F.simplify(mode=[trigsimp,ratsimp]).

   Actually ``simplify`` could be used to apply any scalar *sympy* function to the coefficients of the multivector.

``set_coef(self,grade,base,value)``

   Set the multivector coefficient of index ``(grade,base)`` to ``value``.

``subs(self,x)``

   Return multivector where *sympy* subs function has been applied to each coefficient of multivector for argument dictionary/list ``x``.

``trigsimp(self,**kwargs)``

   Apply the ``sympy`` trigonometric simplification function ``trigsimp`` to each coefficient of the multivector. ``**kwargs`` are the arguments of trigsimp. See ``sympy`` documentation on ``trigsimp`` for more information.

Basic Multivector Functions
---------------------------

``com(A,B)``

   Calculate commutator of multivectors :math:`A` and :math:`B`. Returns :math:`(AB-BA)/2`.

   Additionally, commutator and anti-commutator operators are defined by

   .. math::

      $\begin{aligned}
                  \texttt{A >> B} \equiv & {\displaystyle\frac{AB - BA}{2}} \\
                  \texttt{A << B} \equiv & {\displaystyle\frac{AB + BA}{2}}.
          \end{aligned}

``cross(v1,v2)``

   If ``v1`` and ``v2`` are 3-dimensional Euclidean vectors the vector cross product is returned, :math:`v_{1}\times v_{2} = -I{\lp {v_{1}{\wedge}v_{2}} \rp }`.

``def_prec(gd,op_ord='<>|,^,*')``\  [19]_

   This is used with the ``GAeval()`` function to evaluate a string representing a multivector expression with a revised operator precedence. ``def_prec()`` redefines the operator precedence for multivectors. ``def_prec()`` must be called in the main program an the argument ``gd`` must be ``globals()``. The argument ``op_ord`` defines the order of operator precedence from high to low with groups of equal precedence separated by commas. the default precedence ``op_ord='<>|,^,\*'`` is that used by
   Hestenes (:cite:`Hestenes`,p7,:cite:`Doran`,p38).

``dual(A,mode='I+')``

   Return the dual of the multivector ``A``. The default operation is :math:`AI`. For other modes see member function ``Mv.dual(mode)``

``even(A)``

   Return even part of :math:`A`.

``exp(A,hint='-')``

   If :math:`A` is a multivector then ``A.exp(hint)`` is returned. If :math:`A` is a *sympy* expression the *sympy* expression :math:`e^{A}` is returned (see ``sympy.exp(A)`` member function).

``GAeval(s,pstr=False)``\  [20]_

   Returns multivector expression for string ``s`` with operator precedence for string ``s`` defined by inputs to function ``def_prec()``. if ``pstr=True`` ``s`` and ``s`` with parenthesis added to enforce operator precedence are printed.

``grade(A,r=0)``

   If :math:`A` is a multivector :math:`{\left < {A} \right >}_{r}` is returned.

``inv(A)``

   If :math:`A` is a multivector and :math:`AA^{{\dagger}}` is a non-zero scalar then :math:`A^{-1} = A^{{\dagger}}/(AA^{{\dagger}})` is returned otherwise an exception is returned.

``Nga(x,prec=5)``

   If ``x`` is a multivector with coefficients that contain floating point numbers, ``Nga()`` rounds all these numbers to a precision of ``prec`` and returns the rounded multivector.

``norm(A,hint='-')``

   If :math:`A` is a multivector and :math:`AA^{{\dagger}}` is a number (not a scalar function) then :math:`\sqrt{{\left |{AA^{{\dagger}}}\right |}}` is returned. If :math:`AA^{{\dagger}}` is a scalar *sympy* expression, but not a number, and ``hint='-'`` then return :math:`\sqrt{-AA^{{\dagger}}}` otherwise return :math:`\sqrt{AA^{{\dagger}}}`.

``norm2(A)``

   If :math:`A` is a multivector and :math:`AA^{{\dagger}}` is a scalar return :math:`{\left |{AA^{{\dagger}}}\right |}`.

``odd(A)``

   Return odd part of :math:`A`.

``proj(B,A)``

   Project blade ``A`` on blade ``B`` returning :math:`{\lp {A\rfloor B} \rp }B^{-1}`.

``ReciprocalFrame(basis,mode='norm')``

   If ``basis`` is a list/tuple of vectors, ``ReciprocalFrame()`` returns a tuple of reciprocal vectors. If ``mode=norm`` the vectors are normalized. If ``mode`` is anything other than ``norm`` the vectors are unnormalized and the normalization coefficient is added to the end of the tuple. One must divide by this coefficient to normalize the vectors.

``refl(B,A)``

   Reflect multivector :math:`A` in blade :math:`B`. If :math:`s` is grade of :math:`B` returns :math:`\sum_{r}(-1)^{s(r+1)}B{\left < {A} \right >}_{r}B^{-1}`.

``rev(A)``

   If :math:`A` is a multivector return :math:`A^{{\dagger}}`.

``rot(itheta,A,hint='-')``

   If ``A`` is a multivector return ``A.rotate_multivector(itheta,hint)`` where ``itheta`` is the bi-vector blade defining the rotation. For the use of ``hint`` see the member function ``Mv.rotate_multivector(self,itheta,hint)``.

.. _makeMVD:

Multivector Derivatives
-----------------------

The various derivatives of a multivector function is accomplished by multiplying the gradient operator vector with the function. The gradient operation vector is returned by the ``Ga.grads()`` function if coordinates are defined. For example if we have for a 3-D vector space

.. code:: python

   X = (x,y,z) = symbols('x y z')
   o3d = Ga('e*x|y|z',metric='[1,1,1]',coords=X)
   (ex,ey,ez) = o3d.mv()
   (grad,rgrad) = o3d.grads()

Then the gradient operator vector is ``grad`` (actually the user can give it any name he wants to). The derivatives of the multivector function ``F = o3d.mv('F','mv',f=True)`` are given by multiplying by the left geometric derivative operator and the right geometric derivative operator (:math:`\T{grad} = \nabla` and :math:`\T{rgrad} = \bar{\nabla}`). Another option is to use the radiant operator members of the geometric algebra directly where we have :math:`\nabla = {\texttt{o3d.grad}}` and
:math:`\bar{\nabla} = {\texttt{o3d.rgrad}}`.

.. math::

   \begin{aligned}
               \nabla F &=  \texttt{grad*F} \\
               F \bar{\nabla} &=  \texttt{F*rgrad} \\
               \nabla {\wedge}F &=  \texttt{grad^F} \\
               F {\wedge}\bar{\nabla} &=  \texttt{F^rgrad} \\
               \nabla \cdot F &=  \texttt{grad|F} \\
               F \cdot \bar{\nabla} &=  \texttt{F|rgrad} \\
               \nabla \rfloor F &=  \texttt{grad<F} \\
               F \rfloor \bar{\nabla} &=  \texttt{F<rgrad} \\
               \nabla \lfloor F &=  \texttt{grad>F} \\
               F \lfloor \bar{\nabla} &= \texttt{F>rgrad}
         \end{aligned}

The preceding list gives examples of all possible multivector derivatives of the multivector function ``F`` where the operation returns a multivector function. The complementary operations

.. math::

   \begin{aligned}
               F \nabla &=  \texttt{F*grad} \\
               \bar{\nabla} F &=  \texttt{rgrad*F} \\
               F {\wedge}\nabla &=  \texttt{F^grad} \\
               \bar{\nabla} {\wedge}F &=  \texttt{rgrad^F} \\
               F \cdot \nabla &=  \texttt{F|grad} \\
               \bar{\nabla}\cdot F &=  \texttt{rgrad|F} \\
               F \rfloor \nabla &=  \texttt{F<grad} \\
               \bar{\nabla} \rfloor F &=  \texttt{rgrad<F} \\
               F \lfloor \nabla &=  \texttt{F>grad} \\
               \bar{\nabla} \lfloor F &= \texttt{rgrad>F}
         \end{aligned}

all return multivector linear differential operators.

Submanifolds
------------

In general the geometric algebra that the user defines exists on the tangent space of a manifold (see section :ref:`sect_manifold`). The submanifold class, ``Sm``, is derived from the ``Ga`` class and allows one to define a submanifold of a manifold by defining a coordinate mapping between the submanifold coordinates and the manifold coordinates. What is returned as the submanifold is the geometric algebra of the tangent space of the submanifold. The submanifold for a geometric algebra is
instantiated with

``Ga.sm(map,coords,root='e',norm=False)``

   To define the submanifold we must def a coordinate map from the coordinates of the submanifold to each of the coordinates of the base manifold. Thus the arguments ``map`` and ``coords`` are respectively lists of functions and symbols. The list of symbols, ``coords``, are the coordinates of the submanifold and are of length equal to the dimension of the submanifold. The list of functions, ``map``, define the mapping from the coordinate space of the submanifold to the coordinate space of the
   base manifold. The length of ``map`` is equal to the dimension of the base manifold and each function in ``map`` is a function of the coordinates of the submanifold. ``root`` is the root of the string that is used to name the basis vectors of the submanifold. The default value of ``root`` is ``e``. The result of this is that if the *sympy* symbols for the coordinates are ``u`` and ``v`` (two dimensional manifold) the text symbols for the basis vectors are ``e_u`` and ``e_v`` or in LaTeX
   :math:`e_{u}` and :math:`e_{v}`. As a concrete example consider the following code.

   .. literalinclude:: python/submanifold.py

   The output of this program (using LaTeX) is

   |image0|

   The base manifold, ``sp3d``, is a 3-d Euclidean space using standard spherical coordinates. The submanifold ``sph2d`` of ``sp3d`` is a spherical surface of radius :math:`1`. To take the sumanifold operation one step further the submanifold ``cir1d`` of ``sph2d`` is a circle in ``sph2d`` where the latitude of the circle is :math:`\pi/8`.

   In each case, for demonstration purposes, a scalar and vector function on each manifold is defined (``f`` and ``F`` for the 2-d manifold and ``h`` and ``H`` for the 1-d manifold) and the geometric derivative of each function is taken. The manifold mapping and the metric tensor for ``cir1d`` of ``sph2d`` are also shown. Note that if the submanifold basis vectors are not normalized\ [21]_ the program output is

   |image1|

Linear Transformations
----------------------

The mathematical background for linear transformations is in section :ref:`Ltrans`. Linear transformations on the tangent space of the manifold are instantiated with the ``Ga`` member function ``lt`` (the actual class being instantiated is ``Lt``) as shown in lines 12, 20, 26, and 44 of the code listing ``Ltrans.py``. In all of the examples in ``Ltrans.py`` the default instantiation is used which produces a general (all the coefficients of the linear transformation are symbolic constants) linear
transformation. *Note that to instantiate linear transformations coordinates, :math:`{\left \{ {{\eb}_{i}} \rbrc}`, must be defined when the geometric algebra associated with the linear transformation is instantiated. This is due to the naming conventions of the general linear transformation (coordinate names are used) and for the calculation of the trace of the linear transformation which requires taking a divergence.* To instantiate a specific linear transformation the usage of ``lt()`` is
``Ga.lt(M,f=False,mode='g')``

   ``M`` is an expression that can define the coefficients of the linear transformation in various ways defined as follows.

   +----------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``M``                      | Result                                                                                                                                                                                                                                                                                        |
   +============================+===============================================================================================================================================================================================================================================================================================+
   | string ``M``               | Coefficients are symbolic constants with names :math:`\T{M}^{x_{i}x_{j}}` where :math:`x_{i}` and :math:`x_{j}` are the names of the :math:`i^{th}` and :math:`j^{th}` coordinates (see output of ``Ltrans.py``).                                                                             |
   +----------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | char ``mode``              | If ``M`` is a string then ``mode`` determines whether the linear transformation is general, ``mode='g'``, symmetric, ``mode='s'``, or antisymmetric, ``mode='a'``. The default is ``mode='g'``.                                                                                               |
   +----------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | list ``M``                 | If ``M`` is a list of vectors equal in length to the dimension of the vector space then the linear transformation is :math:`\f{L}{\ebf_{i}} = \T{M}\mat{i}`. If ``M``\ is a list of lists of scalars where all lists are equal in length to the dimension of the vector space then the linear |
   |                            | transformation is\ :math:`\f{L}{\ebf_{i}} = \T{M}\mat{i}\mat{j}\ebf_{j}`.                                                                                                                                                                                                                     |
   +----------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | dict ``M``                 | If ``M`` is a dictionary the linear transformation is defined by :math:`\f{L}{\ebf_{i}} = \T{M}\mat{\ebf_{i}}`. If :math:`\ebf_{i}` is not in the dictionary then :math:`\f{L}{\ebf_{i}} =0`.                                                                                                 |
   +----------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | rotor ``M``                | If ``M`` is a rotor, :math:`\T{M}\T{M}^{\R}=1`, the linear transformation is defined by :math:`\f{L}{{\ebf}_{i}} = \T{M}{\ebf}_{i}\T{M}^{\R}` .                                                                                                                                               |
   +----------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | multivector function ``M`` | If ``M`` is a general multivector function, the function is tested for linearity, and if linear the coefficients of the linear transformation are calculated from :math:`\f{L}{\ebf_{i}} = \f{\T{M}}{\ebf_{i}}`.                                                                              |
   +----------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

   ``f`` is ``True`` or ``False``. If ``True`` the symbolic coefficients of the general linear transformation are instantiated as functions of the coordinates.

The different methods of instantiation are demonstrated in the code ``LtransInst.py``

.. literalinclude:: python/LtransInst.py

with output

|image2|

The member function of the ``Lt`` class are

``Lt(A)``

   Returns the image of the multivector :math:`A` under the linear transformation :math:`L` where :math:`{{L}\lp {A} \rp }` is defined by the linearity of :math:`L`, the vector values :math:`{{L}\lp {{{\eb}}_{i}} \rp }`, and the definition :math:`{{L}\lp {{{\eb}}_{i_{1}}{\wedge}\dots{\wedge}{{\eb}}_{i_{r}}} \rp } = {{L}\lp {{{\eb}}_{i_{1}}} \rp }{\wedge}\dots{\wedge}{{L}\lp {{{\eb}}_{i_{r}}} \rp }`.

``Lt.det()``

   Returns the determinant (a scalar) of the linear transformation, :math:`L`, defined by :math:`{{\det}\lp {L} \rp }I = {{L}\lp {I} \rp }`.

``Lt.adj()``

   Returns the adjoint (a linear transformation) of the linear transformation, :math:`L`, defined by :math:`a\cdot{{L}\lp {b} \rp } = b\cdot{{\bar{L}}\lp {a} \rp }` where :math:`a` and :math:`b` are any two vectors in the tangent space and :math:`\bar{L}` is the adjoint of :math:`L`.

``Lt.tr()``

   Returns the trace (a scalar) of the linear transformation, :math:`L`, defined by :math:`{{\operatorname{tr}}\lp {L} \rp }=\nabla_{a}\cdot{{L}\lp {a} \rp }` where :math:`a` is a vector in the tangent space.

``Lt.matrix()``

   Returns the matrix representation (*sympy* ``Matrix``) of the linear transformation, :math:`L`, defined by :math:`{{L}\lp {{{\eb}}_{i}} \rp } = L_{ij}{{\eb}}_{j}` where :math:`L_{ij}` is the matrix representation.

The ``Ltrans.py`` demonstrate the use of the various ``Lt`` member functions and operators. The operators that can be used with linear transformations are ``+``, ``-``, and ``*``. If :math:`A` and :math:`B` are linear transformations, :math:`V` a multivector, and :math:`\alpha` a scalar then :math:`{{{\lp {A\pm B} \rp }}\lp {V} \rp } = {{A}\lp {V} \rp }\pm{{B}\lp {V} \rp }`, :math:`{{{\lp {AB} \rp }}\lp {V} \rp } = {{A}\lp {{{B}\lp {V} \rp }} \rp }`, and
:math:`{{{\lp {\alpha A} \rp }}\lp {V} \rp } = \alpha{{A}\lp {V} \rp }`.

The ``matrix()`` member function returns a *sympy* ``Matrix`` object which can be printed in IPython notebook. To directly print an linear transformation in *ipython notebook* one must implement (yet to be done) a printing method similar to ``mv.Fmt()``.

Note that in ``Ltrans.py`` lines 30 and 49 are commented out since the latex output of those statements would run off the page. The use can uncomment those statements and run the code in the “LaTeX docs” directory to see the output.

.. literalinclude:: python/Ltrans.py

The output of this code is.

|image3|

Differential Operators
----------------------

For the mathematical treatment of linear multivector differential operators see section :ref:`ldops`. The is a differential operator class ``Dop``. However, one never needs to use it directly. The operators are constructed from linear combinations of multivector products of the operators ``Ga.grad`` and ``Ga.rgrad`` as shown in the following code for both orthogonal rectangular and spherical 3-d coordinate systems.

.. literalinclude:: python/Dop.py

The output of this code is.

|image4|

Note that for print an operator in the IPython notebook one must implement (yet to be done) a printing method similar to ``mv.Fmt()``.

Instantiating a Multi-linear Functions (Tensors)
------------------------------------------------

The mathematical background for multi-linear functions is in section :ref:`MLtrans`. To instantiate a multi-linear function use

``Mlt(self, f, Ga, nargs=None, fct=False)``

   Where the arguments are

   +-----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``f``     | Either a string for a general tensor (this option is included mainly for debugging of the ``Mlt`` class) or a multi-linear function of manifold tangent vectors (multi-vectors of grade one) to scalar. For example one could generate a custom |
   |           | python function such as shown in ``TensorDef.py`` .                                                                                                                                                                                             |
   +-----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``Ga``    | Geometric algebra that tensor is associated with.                                                                                                                                                                                               |
   +-----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``nargs`` | If ``f`` is a string then ``nargs`` is the number of vector arguments of the tensor. If ``f`` is anything other than a string ``nargs`` is not required since ``Mlt`` determines the number of vector arguments                                 |
   |           | from ``f``.                                                                                                                                                                                                                                     |
   +-----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``fct``   | If ``f`` is a string then ``fct=True`` forces the tensor to be a tensor field (function of the coordinates. If ``f`` anything other than a string ``fct`` is not required since ``Mlt`` determines whether the                                  |
   |           | tensor is a tensor field from ``f`` .                                                                                                                                                                                                           |
   +-----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. literalinclude:: python/TensorDef.py

Basic Multilinear Function Class Functions
------------------------------------------

If we can instantiate multilinear functions we can use all the multilinear function class functions as described as follows. See section :ref:`MLtrans` for the mathematical description of each operation.

``self(kargs)``

   Calling function to evaluates multilinear function for ``kargs`` list of vector arguments and returns a value. Note that a sympy scalar is returned, *not* a multilinear function.

``self.contract(slot1,slot2)``

   Returns contraction of tensor between ``slot1`` and ``slot2`` where ``slot1`` is the index of the first vector argument and ``slot2`` is the index of the second vector argument of the tensor. For example if we have a rank two tensor, ``T(a1,a2)``, then ``T.contract(1,2)`` is the contraction of ``T``. For this case since there are only two slots there can only be one contraction.

``self.pdiff(slot)``

   Returns gradient of tensor, ``T``, with respect to slot vector. For example if the tensor is :math:`{{T}\lp {a_{1},a_{2}} \rp }` then ``T.pdiff(2)`` is :math:`\nabla_{a_{2}}T`. Since ``T`` is a scalar function, ``T.pdiff(2)`` is a vector function.

``self.cderiv()``

   Returns covariant derivative of tensor field. If ``T`` is a tensor of rank :math:`k` then ``T.cderiv()`` is a tensor of rank :math:`k+1`. The operation performed is defined in section :ref:`MLtrans`.

Standard Printing
-----------------

Printing of multivectors is handled by the module ``printer`` which contains a string printer class derived from the *sympy* string printer class and a latex printer class derived from the *sympy* latex printer class. Additionally, there is an ``Eprint`` class that enhances the console output of *sympy* to make the printed output multivectors, functions, and derivatives more readable. ``Eprint`` requires an ansi console such as is supplied in linux or the program *ConEmu* replaces ``cmd.exe``.

For a windows user the simplest way to implement *ConEmu* is to use the *geany* editor and in the Edit\ :math:`\rightarrow`\ Preferences\ :math:`\rightarrow`\ Tools menu replace ``cmd.exe`` with\ [22]_

``"C:\Program Files\ConEmu\ConEmu64.exe" /WndW 180 /cmd %c``

and then run an example *galgeba* program that used ``Eprint``. The default background and foreground colors make the output unreadable. To change these parameters to reasonable values:\ [23]_

1. Right click on title bar of console.

2. Open *setting* window.

3. Open *colors* window.

4. Set the following parameters to the indicated values:

-  Text: #0
-  Back: #7
-  Popup: #0
-  Back: #7
-  :math:`\rlap{ \checkmark }\square` Extend foreground colors with background #13

If ``Eprint`` is called in a program (linux) when multivectors are printed the basis blades or bases are printed in bold text, functions are printed in red, and derivative operators in green.

For formatting the multivector output there is the member function ``Fmt(self,fmt=1,title=None)`` which is documented in the multivector member functions. This member function works in the same way for LaTeX printing.

There are two functions for returning string representations of multivectors. If ``A`` is a multivector then ``str(A)`` returns a string in which the scalar coefficients of the multivector bases have been simplified (grouped, factored, etc.). The member function ``A.raw_str()`` returns a string in which the scalar coefficients of the multivector bases have not been simplified.

Latex Printing
--------------

For latex printing one uses one functions from the ``ga`` module and one function from the ``printer`` module. The functions are

``Format(Fmode=True,Dmode=True)``

   This function from the ``ga`` module turns on latex printing with the following options

   +-----------+-----------+-----------------------------------------------------------------------------------------+
   | Argument  | Value     | Result                                                                                  |
   +===========+===========+=========================================================================================+
   | ``Fmode`` | ``True``  | Print functions without argument list, :math:`f`                                        |
   +-----------+-----------+-----------------------------------------------------------------------------------------+
   |           | ``False`` | Print functions with standard *sympy* latex formatting, :math:`{{f}\lp {x,y,z} \rp }`   |
   +-----------+-----------+-----------------------------------------------------------------------------------------+
   | ``Dmode`` | ``True``  | Print partial derivatives with condensed notation, :math:`\partial_{x}f`                |
   +-----------+-----------+-----------------------------------------------------------------------------------------+
   |           | ``False`` | Print partial derivatives with standard *sympy* latex formatting, :math:`\pdiff{f}{x}`  |
   +-----------+-----------+-----------------------------------------------------------------------------------------+

   ``Format()`` is also required for printing from *ipython notebook* (note that ``xpdf()`` is not needed to print from *ipython notebook*).

``Fmt(obj,fmt=1)``

   ``Fmt()`` can be used to set the global multivector printing format or to print a tuple, list, of dictionary\ [24]_. The modes and operation of ``Fmt()`` are as follows:

   +---------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
   | ``obj``             | Effect                                                                                                                                      |
   +=====================+=============================================================================================================================================+
   | ``obj=1,2,3``       | Global multivector format is set to 1, 2, or 3 depending on ``obj``. See multivector member function ``Fmt()`` for effect of ``obj`` value. |
   +---------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
   | obj=tuple/list/dict | The printing format of an object that is a tuple, list, or dict is controlled by the ``fmt`` argument in ``Fmt`` :                          |
   +---------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
   |                     | ``fmt=1``: Print complete ``obj`` on one line.                                                                                              |
   +---------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
   |                     | ``fmt=2``: Print one element of ``obj`` on each line.                                                                                       |
   +---------------------+---------------------------------------------------------------------------------------------------------------------------------------------+

``xpdf(filename=None,debug=False,paper=(14,11),crop=False)``

   This function from the ``printer`` module post-processes the output captured from print statements, writes the resulting latex strings to the file ``filename``, processes the file with pdflatex, and displays the resulting pdf file. All latex files except the pdf file are deleted. If ``debug = True`` the file ``filename`` is printed to standard output for debugging purposes and ``filename`` (the tex file) is saved. If ``filename`` is not entered the default filename is the root name of the
   python program being executed with ``.tex`` appended. The ``paper`` option defines the size of the paper sheet for latex. The format for the ``paper`` is

   ===================== =============================================================
   ``paper=(w,h)``       ``w`` is paper width in inches and
   \                     ``h`` is paper height in inches
   ``paper='letter'``    paper is standard letter size 8.5 in :math:`\times` 11 in
   ``paper='landscape'`` paper is standard letter size but 11 in :math:`\times` 8.5 in
   ===================== =============================================================

   The default of ``paper=(14,11)`` was chosen so that long multivector expressions would not be truncated on the display.

   If the ``crop`` input is ``True`` the linux ``pdfcrop`` program is used to crop the pdf output (if output is one page). This only works for linux installations (where ``pdfcrop`` is installed).

   The ``xpdf`` function requires that latex and a pdf viewer be installed on the computer.

   ``xpdf`` *is not required when printing latex in IPython notebook.*

As an example of using the latex printing options when the following code is executed

.. code:: python

   from printer import Format, xpdf
   from ga import Ga
   Format()
   g3d = Ga('e*x|y|z')
   A = g3d.mv('A','mv')
   print r'\bm{A} =',A
   print A.Fmt(2,r'\bm{A}')
   print A.Fmt(3,r'\bm{A}')
   xpdf()

The following is displayed

.. math::

   \begin{aligned}
         {\boldsymbol{A}} = & A+A^{x}{\boldsymbol{e_{x}}}+A^{y}{\boldsymbol{e_{y}}}+A^{z}{\boldsymbol{e_{z}}}+A^{xy}{\boldsymbol{e_{x}{\wedge}e_{y}}}+A^{xz}{\boldsymbol{e_{x}{\wedge}e_{z}}}+A^{yz}{\boldsymbol{e_{y}{\wedge}e_{z}}}+A^{xyz}{\boldsymbol{e_{x}{\wedge}e_{y}{\wedge}e_{z}}} \\
         {\boldsymbol{A}} =  & A \\  & +A^{x}{\boldsymbol{e_{x}}}+A^{y}{\boldsymbol{e_{y}}}+A^{z}{\boldsymbol{e_{z}}} \\  & +A^{xy}{\boldsymbol{e_{x}{\wedge}e_{y}}}+A^{xz}{\boldsymbol{e_{x}{\wedge}e_{z}}}+A^{yz}{\boldsymbol{e_{y}{\wedge}e_{z}}} \\  & +A^{xyz}{\boldsymbol{e_{x}{\wedge}e_{y}{\wedge}e_{z}}} \\
         {\boldsymbol{A}} =  & A \\  & +A^{x}{\boldsymbol{e_{x}}} \\  & +A^{y}{\boldsymbol{e_{y}}} \\  & +A^{z}{\boldsymbol{e_{z}}} \\  & +A^{xy}{\boldsymbol{e_{x}{\wedge}e_{y}}} \\  & +A^{xz}{\boldsymbol{e_{x}{\wedge}e_{z}}} \\  & +A^{yz}{\boldsymbol{e_{y}{\wedge}e_{z}}} \\  & +A^{xyz}{\boldsymbol{e_{x}{\wedge}e_{y}{\wedge}e_{z}}}\end{aligned}

For the cases of derivatives the code is

.. code:: python

   from printer import Format, xpdf
   from ga import Ga

   Format()
   X = (x,y,z) = symbols('x y z')
   o3d = Ga('e_x e_y e_z',g=[1,1,1],coords=X)

   f = o3d.mv('f','scalar',f=True)
   A = o3d.mv('A','vector',f=True)
   B = o3d.mv('B','grade2',f=True)

   print r'\bm{A} =',A
   print r'\bm{B} =',B

   print 'grad*f =',o3d.grad*f
   print r'grad|\bm{A} =',o3d.grad|A
   (o3d.grad*A).Fmt(2,r'grad*\bm{A}')

   print r'-I*(grad^\bm{A}) =',-o3g.mv_I*(o3d.grad^A)
   print (o3d.grad*B).Fmt(2,r'grad*\bm{B}')
   print r'grad^\bm{B} =',o3d.grad^B
   print r'grad|\bm{B} =',o3d.grad|B

   xpdf()

and the latex displayed output is (:math:`f` is a scalar function)

.. math:: \be {\boldsymbol{A}} = A^{x}{\boldsymbol{e_{x}}}+A^{y}{\boldsymbol{e_{y}}}+A^{z}{\boldsymbol{e_{z}}} \ee

.. math:: \be {\boldsymbol{B}} = B^{xy}{\boldsymbol{e_{x}{\wedge}e_{y}}}+B^{xz}{\boldsymbol{e_{x}{\wedge}e_{z}}}+B^{yz}{\boldsymbol{e_{y}{\wedge}e_{z}}} \ee

.. math:: \be {\boldsymbol{\nabla}}  f = \partial_{x} f{\boldsymbol{e_{x}}}+\partial_{y} f{\boldsymbol{e_{y}}}+\partial_{z} f{\boldsymbol{e_{z}}} \ee

.. math:: \be {\boldsymbol{\nabla}} \cdot {\boldsymbol{A}} = \partial_{x} A^{x} + \partial_{y} A^{y} + \partial_{z} A^{z} \ee

.. math::

   \begin{aligned}
    {\boldsymbol{\nabla}}  {\boldsymbol{A}} =  & \partial_{x} A^{x} + \partial_{y} A^{y} + \partial_{z} A^{z} \\  & +\lp - \partial_{y} A^{x} + \partial_{x} A^{y}\rp {\boldsymbol{e_{x}{\wedge}e_{y}}}+\lp - \partial_{z} A^{x} + \partial_{x} A^{z}\rp {\boldsymbol{e_{x}{\wedge}e_{z}}}+\lp - \partial_{z} A^{y} + \partial_{y} A^{z}\rp {\boldsymbol{e_{y}{\wedge}e_{z}}} \\ \end{aligned}

.. math:: \be -I ({\boldsymbol{\nabla}} {\wedge}{\boldsymbol{A}}) = \lp - \partial_{z} A^{y} + \partial_{y} A^{z}\rp {\boldsymbol{e_{x}}}+\lp \partial_{z} A^{x} - \partial_{x} A^{z}\rp {\boldsymbol{e_{y}}}+\lp - \partial_{y} A^{x} + \partial_{x} A^{y}\rp {\boldsymbol{e_{z}}} \ee

.. math::

   \begin{aligned}
    {\boldsymbol{\nabla}}  {\boldsymbol{B}} =  & \lp - \partial_{y} B^{xy} - \partial_{z} B^{xz}\rp {\boldsymbol{e_{x}}}+\lp \partial_{x} B^{xy} - \partial_{z} B^{yz}\rp {\boldsymbol{e_{y}}}+\lp \partial_{x} B^{xz} + \partial_{y} B^{yz}\rp {\boldsymbol{e_{z}}} \\  & +\lp \partial_{z} B^{xy} - \partial_{y} B^{xz} + \partial_{x} B^{yz}\rp {\boldsymbol{e_{x}{\wedge}e_{y}{\wedge}e_{z}}} \\ \end{aligned}

.. math:: \be {\boldsymbol{\nabla}} {\wedge}{\boldsymbol{B}} = \lp \partial_{z} B^{xy} - \partial_{y} B^{xz} + \partial_{x} B^{yz}\rp {\boldsymbol{e_{x}{\wedge}e_{y}{\wedge}e_{z}}} \ee

.. math:: \be {\boldsymbol{\nabla}} \cdot {\boldsymbol{B}} = \lp - \partial_{y} B^{xy} - \partial_{z} B^{xz}\rp {\boldsymbol{e_{x}}}+\lp \partial_{x} B^{xy} - \partial_{z} B^{yz}\rp {\boldsymbol{e_{y}}}+\lp \partial_{x} B^{xz} + \partial_{y} B^{yz}\rp {\boldsymbol{e_{z}}} \ee

This example also demonstrates several other features of the latex printer. In the case that strings are input into the latex printer such as ``r'grad*\bm{A}'``, ``r'grad^\bm{A}'``, or ``r'grad*\bm{A}'``. The text symbols ``grad``, ``^``, ``|``, and ``*`` are mapped by the ``xpdf()`` post-processor as follows if the string contains an ``=``.

========== =================== ==============================
original   replacement         displayed latex
``grad*A`` ``\bm{\nabla}A``    :math:`{\boldsymbol{\nabla}}A`
``A^B``    ``A\wedge B``       :math:`A\wedge B`
``A|B``    ``A\cdot B``        :math:`A\cdot B`
``A*B``    ``AB``              :math:`AB`
``A<B``    ``A\rfloor B``      :math:`A\rfloor B`
``A>B``    ``A\lfloor B``      :math:`A\lfloor B`
``A>>B``   ``A\times B``       :math:`A\times B`
``A<<B``   ``A\bar{\times} B`` :math:`A\bar{\times} B`
========== =================== ==============================

If the first character in the string to be printed is a ``%`` none of the above substitutions are made before the latex processor is applied. In general for the latex printer strings are assumed to be in a math environment (equation or align) unless the first character in the string is a ``#``\ [25]_.

There are two member functions for returning LaTeX string representations of multivectors. If ``A`` is a multivector then ``A.Mv_latex_str()`` returns a LaTeX string in which the scalar coefficients of the multivector bases have been simplified (grouped, factored, etc.). This function is used when using ``print`` in the LaTeX mode. The member function ``A.raw_latex_str()`` returns a LaTeX string in which the scalar coefficients of the multivector bases have not been simplified.

Printing Lists/Tuples of Multivectors/Differential Operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since the expressions for multivectors or differential operators can be very long printing lists or tuples of such items can easily exceed the page with when printing in LaTeX or in “ipython notebook.” I order to alleviate this problem the function ``Fmt`` can be used.

``Fmt(obj,fmt=0)``

   This function from the ``printer`` module allows the formatted printing of lists/tuples or multivectors/differential operators.

   +-----------+---------------------------------------------------------------------------------+
   | ``obj``   | ``obj`` is a list or tuple of multivectors and/or differential operators.       |
   +-----------+---------------------------------------------------------------------------------+
   | ``fmt=0`` | ``fmt=0`` prints each element of the list/tuple on an individual lines\ [26]_.  |
   +-----------+---------------------------------------------------------------------------------+
   |           | ``fmt=1`` prints all elements of the list/tuple on a single line\ [26]_.        |
   +-----------+---------------------------------------------------------------------------------+

   If l is a list or tuple to print in the LaTeX environment use the command

   .. code:: python

      print Fmt(l) # One element of l per line

   or

   .. code:: python

      print Fmt(l,1) # All elements of l on one line

   If you are printing in “ipython notebook” then enter

   .. code:: python

      Fmt(l) # One element of l per line

   or

..

   .. code:: python

      Fmt(l,1) # All elements of l on one line

--------------

.. [4]
   By the manifold embedding theorem any :math:`m`-dimensional manifold is isomorphic to a :math:`m`-dimensional vector manifold

.. [5]
   This product in not necessarily positive definite.

.. [6]
   In this section and all following sections we are using the Einstein summation convention unless otherwise stated.

.. [7]
   We use the Christoffel symbols of the first kind to calculate the derivatives of the basis vectors and the product rule to calculate the derivatives of the basis blades where (http://en.wikipedia.org/wiki/Christoffel_symbols)

   .. math:: \be \Gamma_{ijk} = {\frac{1}{2}}{\lp {{{\displaystyle\frac{\partial {g_{jk}}}{\partial {x^{i}}}}}+{{\displaystyle\frac{\partial {g_{ik}}}{\partial {x^{j}}}}}-{{\displaystyle\frac{\partial {g_{ij}}}{\partial {x^{k}}}}}} \rp }, \ee

   and

   .. math:: \be {{\displaystyle\frac{\partial {{{\eb}}_{j}}}{\partial {x^{i}}}}} = \Gamma_{ijk}{{\eb}}^{k}. \ee

   The Christoffel symbols of the second kind,

   .. math:: \be \Gamma_{ij}^{k} = {\frac{1}{2}}g^{kl}{\lp {{{\displaystyle\frac{\partial {g_{li}}}{\partial {x^{j}}}}}+{{\displaystyle\frac{\partial {g_{lj}}}{\partial {x^{i}}}}}-{{\displaystyle\frac{\partial {g_{ij}}}{\partial {x^{l}}}}}} \rp }, \ee

   could also be used to calculate the derivatives in term of the original basis vectors, but since we need to calculate the reciprocal basis vectors for the geometric derivative it is more efficient to use the symbols of the first kind.

.. [8]
   In this case :math:`D_{B}^{j_{1}\dots j_{n}} = F` and :math:`\partial_{j_{1}\dots j_{n}} = 1`.

.. [9]
   For example in three dimensions :math:`{\left \{{3} \rbrc} = (0,1,2,3,(1,2),(2,3),(1,3),(1,2,3))` and as an example of how the superscript would work with each grade :math:`{{\eb}}^{0}=1`, :math:`{{\eb}}^{1}={{\eb}}^{1}`, :math:`{{\eb}}^{{\lp {1,2} \rp }}={{\eb}}^{1}{\wedge}{{\eb}}^{2}`, and :math:`{{\eb}}^{{\lp {1,2,3} \rp }}={{\eb}}^{1}{\wedge}{{\eb}}^{2}{\wedge}{{\eb}}^{3}`.

.. [10]
   We are following the treatment of Tensors in section 3–10 of :cite:`Hestenes`.

.. [11]
   We assume that the arguments are elements of a vector space or more generally a geometric algebra so that the concept of linearity is meaningful.

.. [12]
   Since ``X`` or the metric tensor can be functions of coordinates the vector space that the geometric algebra is constructed from is not necessarily flat so that the geometric algebra is actually constructed on the tangent space of the manifold which is a vector space.

.. [13]
   The signature of the vector space, :math:`(p,q)`, is required to determine whether the square of the normalized pseudoscalar, :math:`I`, is :math:`+1` or :math:`-1`. In the future the metric tensor would be required to create a generalized spinor (:cite:`Hestenes`, pg106).

.. [14]
   Using LaTeX output if a basis vector is denoted by :math:`{{\eb}}_{x}` then :math:`{{\eb}}` is the root symbol and :math:`x` is the subscript

.. [15]
   There is a multivector class, ``Mv``, but in order the insure that every multivector is associated with the correct geometric algebra we always use the member function ``Ga.mv`` to instantiate the multivector.

.. [16]
   Denoted in text output by ``A__x``, etc. so that for text output ``A`` would be printed as ``A__x*e_x+A__y*e_y+A__z*e_z``.

.. [17]
   If the metric is input as a list or list or lists the object is no longer quoted (input as a string). For example the old ``metric='[1,1,1]'`` becomes ``metric=[1,1,1]``.

.. [18]
   In the future it should be possible to generate closed form expressions for :math:`e^{A}` if :math:`A^{r}` is a scalar for some interger :math:`r`.

.. [19]
   See footnote [6]_.

.. [20]
   ``GAeval`` is in the ``printer`` module\ [6]_.

.. [21]
   Remember that normalization is currently supported only for orthogonal systems (diagonal metric tensors).

.. [22]
   The 180 in the *ConEmu* command line is the width of the console you wish to display in characters. Change the number to suit you.

.. [23]
   I am not exactly sure what the different parameter setting do. I achieved the result I wished for by trial and error. I encourage the users to experiment and share their results.

.. [24]
   In *Ipython notebook* tuples, or lists, or dictionarys of multivectors do print correctly. One mode of ``Fmt()`` corrects this deficiency.

.. [25]
   Preprocessing do not occur for the Ipython notebook and the string post processing commands ``%`` and ``#`` are not used in this case.

.. [26]
   The formatting of each element is respected as applied by ``A.Fmt(fmt=1,2, or 3)`` where ``A`` is an element of ``obj``\ so that if multivector/differential operation have been formatted to print on multiple lines it will printed on multiple lines.

.. |image0| image:: images/submanifold.svg
.. |image1| image:: images/submanifold1.svg
.. |image2| image:: images/LtransInst.svg
.. |image3| image:: images/Ltrans.svg
.. |image4| image:: images/Dop.svg