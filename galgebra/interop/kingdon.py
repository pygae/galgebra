"""
Kingdon-compatible ``Cl(p, q, r)`` factory.

Uses kingdon conventions:
- Dual mode ``'Iinv+'`` (polarity dual: ``dual(A) = A * I^{-1}``)
- 0-indexed basis names for PGA (``e0, e1, ...``)
"""

from ..ga import Ga


def Cl(p: int, q: int = 0, r: int = 0, root: str = 'e', **kwargs):
    r"""
    Construct a Clifford algebra :math:`Cl(p,q,r)` using kingdon conventions.

    Differences from ``galgebra.interop.Cl``:

    - **Dual mode**: sets ``'Iinv+'`` globally so that
      ``dual(A) = A * I^{-1}``, matching kingdon's ``codegen_polarity``.
    - **Basis naming**: uses 0-indexed names (``e_0, e_1, ...``) to match
      kingdon's default for PGA algebras.

    .. note::

        This function sets the session-wide dual mode to ``'Iinv+'``
        before building the algebra.  ``galgebra.interop.Cl`` resets it
        back to ``'I+'``, so mixing the two in one session is safe as long
        as each call is followed by the code that uses that algebra before
        the next ``Cl`` call.  For full isolation, save and restore
        ``Ga.dual_mode_value`` manually::

            saved = Ga.dual_mode_value
            try:
                ga, *basis = Cl(3, 0, 1)
                # ... kingdon-convention code ...
            finally:
                Ga.dual_mode(saved)

    Parameters
    ----------
    p : int
        Number of basis vectors that square to +1.
    q : int
        Number of basis vectors that square to -1.
    r : int
        Number of basis vectors that square to 0 (degenerate).
    root : str
        Root name for basis vectors (default ``'e'``).
    **kwargs :
        Additional keyword arguments passed to :class:`Ga`.

    Returns
    -------
    tuple
        ``(ga, e_0, ..., e_{n-1})`` -- the geometric algebra and basis vectors.

    Examples
    --------
    3D Euclidean algebra::

        >>> from galgebra.interop.kingdon import Cl
        >>> ga, e1, e2, e3 = Cl(3)

    3D PGA (0-indexed, degenerate metric)::

        >>> ga, e0, e1, e2, e3 = Cl(3, 0, 1)
    """
    n = p + q + r
    if n == 0:
        raise ValueError("Total dimension p + q + r must be positive.")

    # Set kingdon dual convention
    Ga.dual_mode('Iinv+')

    # Build diagonal metric: p ones, q negative ones, r zeros
    g = [1] * p + [-1] * q + [0] * r

    # Build basis name string (0-indexed, matching kingdon)
    if n <= 10:
        basis = root + '*' + '|'.join(str(i) for i in range(n))
    else:
        basis = ' '.join(root + '_' + str(i) for i in range(n))

    return Ga.build(basis, g=g, **kwargs)
