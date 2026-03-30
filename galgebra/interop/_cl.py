"""
Core ``Cl(p, q, r)`` factory with galgebra defaults.
"""

from ..ga import Ga


def Cl(p: int, q: int = 0, r: int = 0, root: str = 'e', **kwargs):
    r"""
    Construct a Clifford algebra :math:`Cl(p,q,r)` from its signature.

    This interface was popularized by ganja.js and adopted by kingdon.
    It provides a concise way to create a geometric algebra without
    specifying a full metric tensor.

    Uses galgebra defaults (basis numbered from 1, dual mode ``'I+'``).
    For kingdon conventions, use ``galgebra.interop.kingdon.Cl`` instead.

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
        ``(ga, e_1, ..., e_n)`` -- the geometric algebra and basis vectors.

    Examples
    --------
    3D Euclidean algebra::

        >>> ga, e1, e2, e3 = Cl(3)

    Spacetime algebra (Minkowski, +---)::

        >>> ga, e1, e2, e3, e4 = Cl(1, 3)

    Projective Geometric Algebra::

        >>> ga, e1, e2, e3 = Cl(2, 0, 1)
    """
    n = p + q + r
    if n == 0:
        raise ValueError("Total dimension p + q + r must be positive.")

    # Build diagonal metric: p ones, q negative ones, r zeros
    g = [1] * p + [-1] * q + [0] * r

    # Build basis name string (1-indexed)
    if n <= 9:
        basis = root + '*' + '|'.join(str(i + 1) for i in range(n))
    else:
        basis = ' '.join(root + '_' + str(i + 1) for i in range(n))

    return Ga.build(basis, g=g, **kwargs)
