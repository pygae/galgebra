Module Components
=================

.. warning::

   This page is converted from the original :math:`\LaTeX` documentation, but
   may no longer reflect the current state of the library. See the API docs at
   :doc:`api` for more up-to-date but less structured and in-depth descriptions.

   If you would like to help with merging the descriptions on the API with a
   description on this page, please head over to :issue:`300` on GitHub, where
   there's an explanation of how to do so. Even merging just one function
   explanation helps!

   Any function or class below with a ``[source]`` link to its right is
   guaranteed to be up-to-date already, as its documentation is identical on
   both pages.

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

.. class:: Ga(basis,g=None,coords=None,X=None,norm=False,sig='e',Isq='-',wedge=True,debug=False)
   :noindex:

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

.. method:: Ga.mvr(norm='True')
   :noindex:

   ``Ga.mvr(norm)`` returns the reciprocal basis vectors as a tuple. This allows the programmer to attach any python variable names to the reciprocal basis vectors that is convenient. For example (demonstrating the use of both ``mv()`` and ``mvr()``)

   .. code:: python

      (e_x,e_y,e_z) = o3d.mv()
      (e__x,e__y,e__z) = o3d.mvr()

   If ``norm='True'`` or the basis vectors are orthogonal the dot product of the basis vector and the corresponding reciprocal basis vector is one :math:`{\lp {e_{i}\cdot e^{j}=\delta_{i}^{j}} \rp }`. If ``norm='False'`` and the basis is non-orthogonal The dot product of the basis vector and the corresponding reciprocal basis vector is the square of the pseudo scalar, :math:`I^{2}`, of the geometric algebra :math:`{\lp {e_{i}\cdot e^{j}=E^{2}\delta_{i}^{j}} \rp }`.

In addition to the basis vectors, if coordinates are defined for the geometric algebra, the left and right geometric derivative operators are calculated and accessed with the ``Ga`` member function ``grads()``.

.. method:: Ga.grads()
   :noindex:

   ``Ga.grads()`` returns a tuple with the left and right geometric derivative operators. A typical usage would be

   .. code:: python

      (grad,rgrad) = sp3d.grads()

   for the spherical 3-d geometric algebra. The left derivative :math:`{\lp {{\texttt{grad}} ={\boldsymbol{\nabla}}} \rp }` and the right derivative :math:`{\lp {{\texttt{rgrad}} = {\boldsymbol{\bar{\nabla}}}} \rp }` have been explained in section :ref:`ldops`. Again the names ``grad`` and ``rgrad`` used in a program are whatever the user chooses them to be. In the previous example ``grad`` and ``rgrad`` are used.

an alternative instantiation method is

.. method:: Ga.build(basis, g=None, coords=None, X=None, norm=False, debug=False)
   :noindex:

   The input parameters for ``Ga.build()`` are the same as for ``Ga()``. The difference is that in addition to returning the geometric algebra ``Ga.build()`` returns the basis vectors at the same time. Using ``Ga.build()`` in the previous example gives

   .. code:: python

      (r, th, phi) = coords = symbols('r,theta,phi', real=True)
      basis = 'e_r e_theta e_phi'
      g = [1, r**2, r**2*sin(th)**2]
      (sp3d,er,eth,ephi) = Ga.build(basis,g=g,coord=coords,norm=True)

To access the pseudo scalar of the geometric algebra use the member function ``I()``.

.. method:: Ga.I()
   :noindex:

   ``Ga.I()`` returns the normalized pseudo scalar :math:`{\lp {{\left |{I^{2}}\right |}=1} \rp }` for the geometric algebra. For example :math:`I = \mbox{{\texttt{o3d.I()}}}` for the ``o3d`` geometric algebra. This function requires the signature of the vector space (see instantiating a geometric algebra).

.. method:: Ga.E()
   :noindex:

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

.. method:: Ga.mv(name, mode, f=False)
   :noindex:

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

.. method:: MV.setup(basis, metric=None, coords=None, rframe=False, debug=False, curv=(None,None))
   :noindex:

   This function allows a single geometric algebra to be created. If the function is called more than once the old geometric algebra is overwritten by the new geometric algebra. The named input ``metric`` is the same as the named input ``g`` in the current version of *galgebra*. Likewise, ``basis``, ``coords``, and ``debug`` are the same in the old and current versions of *galgebra*\ \ [17]_. Due to improvements in *sympy* the inputs ``rframe`` and ``curv[1]`` are no longer required. ``curv[0]`` is
   the vector function (list or tuple of scalar functions) of the coordinates required to define a vector manifold. For compatibility with the old version of *galgebra* if ``curv`` is used ``metric`` should be a orthonormal Euclidean metric of the same dimension as ``curv[0]``. It is strongly suggested that one use the new methods of defining a geometric algebra on a manifold.

.. class:: MV(base, mvtype, fct=False, blade_rep=True)
   :noindex:

   For the instantiation of multivector using ``MV`` the ``base`` and ``mvtype`` arguments are the same as for new methods of multivector instantiation. The ``fct`` input is the same and the ``g`` input in the new methods. ``blade_rep`` is not used in the new methods so setting ``blade_rep=False`` will do nothing. Effectively ``blade_rep=False`` was not used in the old examples.

.. method:: MV.Fmt(self, fmt=1, title=None)
   :noindex:

   ``Fmt`` in ``MV`` has inputs identical to ``Fmt`` in ``Mv`` except that if ``A`` is a multivector then ``A.Fmt(2,'A')`` executes a print statement from ``MV`` and returns ``None``, while from ``Mv``, ``A.Fmt(2,'A')`` returns a string so that the function is compatible with use in *ipython notebook*.

Basic Multivector Class Functions
---------------------------------

If we can instantiate multivectors we can use all the multivector class functions as described as follows.

.. method:: Mv.blade_coefs(self,basis_lst)
   :noindex:

   Find coefficients (sympy expressions) of multivector basis blade expansion corresponding to basis blades in ``basis_lst``. For example if :math:`V = V^{x}{{\eb}}_{x}+V^{y}{{\eb}}_{x}+V^{z}{{\eb}}_{x}` Then :math:`V\text{.blade_coefs}([{{\eb}}_{z},{{\eb}}_{x}]) = [V^{z},V^{x}]` or if :math:`B = B^{xy}{{\eb}}_{x}{\wedge}{{\eb}}_{y}+V^{yz}{{\eb}}_{y}{\wedge}{{\eb}}_{z}` then :math:`B\text{.blade_coefs}([{{\eb}}_{x}{\wedge}{{\eb}}_{y}]) = [B^{xy}]`.

.. method:: Mv.convert_to_blades(self)
   :noindex:

   Convert multivector from the base representation to the blade representation. If multivector is already in blade representation nothing is done.

.. method:: Mv.convert_from_blades(self)
   :noindex:

   Convert multivector from the blade representation to the base representation. If multivector is already in base representation nothing is done.

.. method:: Mv.diff(self,var)
   :noindex:

   Calculate derivative of each multivector coefficient with respect to variable ``var`` and form new multivector from coefficients.

.. method:: Mv.dual(self)
   :noindex:

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

.. method:: Mv.even(self)
   :noindex:

   Return the even grade components of the multivector.

.. method:: Mv.exp(self,hint='-')
   :noindex:

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

.. method:: Mv.expand(self)
   :noindex:

   Return multivector in which each coefficient has been expanded using *sympy* ``expand()`` function.

.. method:: Mv.factor(self)
   :noindex:

   Apply the ``sympy`` ``factor`` function to each coefficient of the multivector.

.. method:: Mv.Fmt(self, fmt=1,title=None)
   :noindex:

   Fuction to print multivectors in different formats where

   ======= ============================================
   ``fmt`` 
   ======= ============================================
   1       Print entire multivector on one line.
   2       Print each grade of multivector on one line.
   3       Print each base of multivector on one line.
   ======= ============================================

   ``title`` appends a title string to the beginning of the output. An equal sign in the title string is not required, but is added as a default. Note that ``Fmt`` only overrides the the global multivector printing format for the particular instance being printed. To reset the global multivector printing format use the function ``Fmt()`` in the printer module.

.. method:: Mv.func(self,fct)
   :noindex:

   Apply the ``sympy`` scalar function ``fct`` to each coefficient of the multivector.

.. method:: Mv.grade(self,igrade=0)
   :noindex:

   Return a multivector that consists of the part of the multivector of grade equal to ``igrade``. If the multivector has no ``igrade`` part return a zero multivector.

.. method:: Mv.inv(self)
   :noindex:

   Return the inverse of the multivector :math:`M` (``M.inv()``). If :math:`M` is a non-zero scalar return :math:`1/M`. If :math:`M^{2}` is a non-zero scalar return :math:`M/{\lp {M^{2}} \rp }`, If :math:`MM^{{\dagger}}` is a non-zero scalar return :math:`M^{{\dagger}}/{\lp {MM^{{\dagger}}} \rp }`. Otherwise exit the program with an error message.

   All division operators (``/``, ``/=``) use right multiplication by the inverse.

.. method:: Mv.norm(self,hint='+')
   :noindex:

   Return the norm of the multivector :math:`M` (``M.norm()``) defined by :math:`\sqrt{{\left |{MM^{{\dagger}}}\right |}}`. If :math:`MM^{{\dagger}}` is a scalar (a *sympy* scalar is returned). If :math:`MM^{{\dagger}}` is not a scalar the program exits with an error message. If :math:`MM^{{\dagger}}` is a number *sympy* can determine if it is positive or negative and calculate the absolute value. If :math:`MM^{{\dagger}}` is a *sympy* expression (function) *sympy* cannot determine the sign of
   the expression so that ``hint='+'`` or ``hint='-'`` is needed to determine if the program should calculate :math:`\sqrt{MM^{{\dagger}}}` or :math:`\sqrt{-MM^{{\dagger}}}`. For example if we are in a Euclidean space and ``M`` is a vector then ``hint='+'``, if ``M`` is a bivector then let ``hint='-'``. If ``hint='0'`` and :math:`MM^{{\dagger}}` is a symbolic scalar ``sqrt(Abs(M*M.rev()))`` is returned where ``Abs()`` is the *sympy* symbolic absolute value function.

.. method:: Mv.norm2(self)
   :noindex:

   Return the the scalar defined by :math:`MM^{{\dagger}}` if :math:`MM^{{\dagger}}` is a scalar. If :math:`MM^{{\dagger}}` is not a scalar the program exits with an error message.

.. method:: Mv.proj(self,bases_lst)
   :noindex:

   Return the projection of the multivector :math:`M` (``M.proj(bases_lst)``) onto the subspace defined by the list of bases (``bases_lst``).

.. method:: Mv.proj(self,lst)
   :noindex:

   Return the projection of the mutivector :math:`A` onto the list, :math:`lst`, of basis blades. For example if :math:`A = A^{x}{{\eb}}_{x}+A^{y}{{\eb}}_{y}+A^{z}{{\eb}}_{z}` then :math:`A.proj{\lp {[{{\eb}}_{x},{{\eb}}_{y}]} \rp } = A^{x}{{\eb}}_{x}+A^{y}{{\eb}}_{y}`. Similarly if :math:`A = A^{xy}{{\eb}}_{x}{\wedge}{{\eb}}_{y}+A^{yz}{{\eb}}_{y}{\wedge}{{\eb}}_{z}` then :math:`A.proj{\lp {[{{\eb}}_{x}{\wedge}{{\eb}}_{y}]} \rp } = A^{xy}{{\eb}}_{x}{\wedge}{{\eb}}_{y}`.

.. method:: Mv.project_in_blade(self,blade)
   :noindex:

   Return the projection of the mutivector :math:`A` in subspace defined by the blade, :math:`B`, using the formula :math:`{\lp {A\rfloor B} \rp }B^{-1}` in :cite:`Macdonald1`, page 121.

.. method:: Mv.pure_grade(self)
   :noindex:

   If the multivector :math:`A` is pure (only contains one grade) return, :math:`A.pure\_grade()`, the index ('0' for a scalar, '1' for vector, '2' for a bi-vector, etc.) of the non-zero grade. If :math:`A` is not pure return the negative of the highest non-zero grade index.

.. method:: Mv.odd(self)
   :noindex:

   Return odd part of multivector.

.. method:: Mv.reflect_in_blade(self,blade)
   :noindex:

   Return the reflection of the mutivector :math:`A` in the subspace defined by the :math:`r`-grade blade, :math:`B_{r}`, using the formula (extended to multivectors) :math:`\sum_{i} {\lp {-1} \rp }^{r{\lp {i+1} \rp }}{B}_{r}{\left < {A} \right >}_{i}B_{r}^{-1}` in :cite:`Macdonald1`, page 129.

.. method:: Mv.rev(self)
   :noindex:

   Return the reverse of the multivector.

.. method:: Mv.rotate_multivector(self,itheta,hint='-')
   :noindex:

   Rotate the multivector :math:`A` via the operation :math:`e^{-\theta i/2}Ae^{\theta i/2}` where itheta = :math:`\theta i`, :math:`\theta` is a scalar, and :math:`i` is a unit, :math:`i^{2} = \pm 1`, 2-blade. If :math:`{\lp {\theta i} \rp }^{2}` is not a number ``hint`` is required to determine the sign of the square of ``itheta``. The default chosen, ``hint='-'``, is correct for any Euclidean space.

.. method:: Mv.scalar(self)
   :noindex:

   Return the coefficient (*sympy* scalar) of the scalar part of a multivector.

.. method:: Mv.simplify(self,mode=simplify)
   :noindex:

   ``mode`` is a *sympy* simplification function of a list/tuple of *sympy* simplification functions that are applied in sequence (if more than one function) each coefficient of the multivector. For example if we wished to applied ``trigsimp`` and ``ratsimp`` *sympy* functions to the multivector ``F`` the code would be

   .. code:: python

      Fsimp = F.simplify(mode=[trigsimp,ratsimp]).

   Actually ``simplify`` could be used to apply any scalar *sympy* function to the coefficients of the multivector.

.. method:: Mv.set_coef(self,grade,base,value)
   :noindex:

   Set the multivector coefficient of index ``(grade,base)`` to ``value``.

.. method:: Mv.subs(self,x)
   :noindex:

   Return multivector where *sympy* subs function has been applied to each coefficient of multivector for argument dictionary/list ``x``.

.. method:: Mv.trigsimp(self,**kwargs)
   :noindex:

   Apply the ``sympy`` trigonometric simplification function ``trigsimp`` to each coefficient of the multivector. ``**kwargs`` are the arguments of trigsimp. See ``sympy`` documentation on ``trigsimp`` for more information.

Basic Multivector Functions
---------------------------

.. automethod:: galgebra.ga.Ga.com
   :noindex:

.. autofunction:: galgebra.mv.cross
   :noindex:

.. autofunction:: galgebra.printer.def_prec
   :noindex:

.. autofunction:: galgebra.mv.dual
   :noindex:

.. autofunction:: galgebra.mv.even
   :noindex:

.. autofunction:: galgebra.mv.exp
   :noindex:

.. autofunction:: galgebra.printer.GAeval
   :noindex:

.. autofunction:: galgebra.mv.grade
   :noindex:

.. autofunction:: galgebra.mv.inv
   :noindex:

.. autofunction:: galgebra.mv.Nga
   :noindex:

.. autofunction:: galgebra.mv.norm
   :noindex:

.. autofunction:: galgebra.mv.norm2
   :noindex:

.. autofunction:: galgebra.mv.odd
   :noindex:

.. autofunction:: galgebra.mv.proj
   :noindex:

.. automethod:: galgebra.ga.Ga.ReciprocalFrame(basis,mode='norm')
   :noindex:

.. autofunction:: galgebra.mv.refl(B,A)
   :noindex:

.. autofunction:: galgebra.mv.rev(A)
   :noindex:

.. autofunction:: galgebra.mv.rot(itheta,A,hint='-')
   :noindex:

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

.. method:: Ga.sm(map,coords,root='e',norm=False)
   :noindex:

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

.. method:: Ga.lt(M,f=False,mode='g')
   :noindex:

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

.. method:: Lt.__call__(A)
   :noindex:

   Returns the image of the multivector :math:`A` under the linear transformation :math:`L` where :math:`{{L}\lp {A} \rp }` is defined by the linearity of :math:`L`, the vector values :math:`{{L}\lp {{{\eb}}_{i}} \rp }`, and the definition :math:`{{L}\lp {{{\eb}}_{i_{1}}{\wedge}\dots{\wedge}{{\eb}}_{i_{r}}} \rp } = {{L}\lp {{{\eb}}_{i_{1}}} \rp }{\wedge}\dots{\wedge}{{L}\lp {{{\eb}}_{i_{r}}} \rp }`.

.. method:: Lt.det()
   :noindex:

   Returns the determinant (a scalar) of the linear transformation, :math:`L`, defined by :math:`{{\det}\lp {L} \rp }I = {{L}\lp {I} \rp }`.

.. method:: Lt.adj()
   :noindex:

   Returns the adjoint (a linear transformation) of the linear transformation, :math:`L`, defined by :math:`a\cdot{{L}\lp {b} \rp } = b\cdot{{\bar{L}}\lp {a} \rp }` where :math:`a` and :math:`b` are any two vectors in the tangent space and :math:`\bar{L}` is the adjoint of :math:`L`.

.. method:: Lt.tr()
   :noindex:

   Returns the trace (a scalar) of the linear transformation, :math:`L`, defined by :math:`{{\operatorname{tr}}\lp {L} \rp }=\nabla_{a}\cdot{{L}\lp {a} \rp }` where :math:`a` is a vector in the tangent space.

.. method:: Lt.matrix()
   :noindex:

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

.. class:: Mlt(self, f, Ga, nargs=None, fct=False)
   :noindex:

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

.. method:: Mlt.__call__(kargs)
   :noindex:

   Calling function to evaluates multilinear function for ``kargs`` list of vector arguments and returns a value. Note that a sympy scalar is returned, *not* a multilinear function.

.. method:: Mlt.contract(slot1,slot2)
   :noindex:

   Returns contraction of tensor between ``slot1`` and ``slot2`` where ``slot1`` is the index of the first vector argument and ``slot2`` is the index of the second vector argument of the tensor. For example if we have a rank two tensor, ``T(a1,a2)``, then ``T.contract(1,2)`` is the contraction of ``T``. For this case since there are only two slots there can only be one contraction.

.. method:: Mlt.pdiff(slot)
   :noindex:

   Returns gradient of tensor, ``T``, with respect to slot vector. For example if the tensor is :math:`{{T}\lp {a_{1},a_{2}} \rp }` then ``T.pdiff(2)`` is :math:`\nabla_{a_{2}}T`. Since ``T`` is a scalar function, ``T.pdiff(2)`` is a vector function.

.. method:: Mlt.cderiv()
   :noindex:

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

.. autofunction:: galgebra.printer.Format
   :noindex:

.. function:: Fmt(obj,fmt=1)
   :noindex:

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

.. function:: xpdf(filename=None,debug=False,paper=(14,11),crop=False)
   :noindex:

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

.. function:: Fmt(obj,fmt=0)
   :noindex:

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

   .. code:: python

      Fmt(l,1) # All elements of l on one line

--------------


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