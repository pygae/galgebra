# Used by

Here is a list of projects that use GAlgebra, divided into two parts:

1. research: used to facilitate or verify research
2. functional: used as a dependency in their major functionality
3. explorational: used as an optional dependency in their examples, exercises, or tutorials

They are roughly ordered by relevance, stars and last updated time.

This data is manually mined from [forks](https://github.com/pygae/galgebra/forks) and [dependents](https://github.com/pygae/galgebra/network/dependents), as ["Used by" button is not shown](https://docs.github.com/en/code-security/supply-chain-security/understanding-your-software-supply-chain/exploring-the-dependencies-of-a-repository#changing-the-used-by-package). All usages are manually reviewed to rule out false positives (e.g. the dependency maybe added, but is not actually used), some usages are manually added to our knowledge.

## Research

(in recent years, non-exhaustive)

- [Geometric Product of Two Oriented Points in Conformal Geometric Algebra (2024)](https://link.springer.com/article/10.1007/s00006-024-01363-6)
  - two anonymous expert reviewers, who either used the python symbolic geometric algebra implementation GAlgebra or the Clifford python library to independently verify the results of this paper
- [Machine Learning Clifford invariants of ADE Coxeter elements (2024)](https://arxiv.org/abs/2310.00041)
  - uses GAlgebra for ADE invariant generation [code](https://github.com/DimaDroid/ML_Clifford_Invariants/blob/main/InvariantData/ADE8_InvariantGeneration.py)
- [Clifford spinors and root system induction: H4 and the Grand Antiprism (2021)](https://arxiv.org/abs/2103.07817)
  - perform practical computations in group theory via versors in Clifford algebra framework
  - LaTeX output of GAlgebra is used in paper [code](https://github.com/ppd22/galgebra-in-action)

## Functional

- Prof. Alan Macdonald's books: [Linear and Geometric Algebra](http://www.faculty.luther.edu/~macdonal/laga/index.html) and [Vector and Geometric Calculus](http://www.faculty.luther.edu/~macdonal/vagc/index.html)
  - GAlgebra is the books' companion library, see [GAlgebraPrimer](http://www.faculty.luther.edu/~macdonal/GAlgebraPrimer.pdf)
- [micahscopes/alglbraic](https://github.com/micahscopes/alglbraic): A python library and CLI utility for generating libraries of algebraic GLSL functions, including Clifford algebras
- [russellgoyder/sundial](https://github.com/russellgoyder/sundial): The sundial problem from a new angle
  - released as package `analemma` on PyPI
  - contains [a cheat sheet for GAlgebra](https://github.com/russellgoyder/geometric-algebra-cheat-sheet)
  - see also https://github.com/pygae/galgebra/issues/506
- [jdekozak/dirac5d](https://github.com/jdekozak/dirac5d): Five dimensional Dirac equation over the reals
  - see also https://github.com/jdekozak/dirac5d/issues/4
- [pygae/GAlgebra.jl](https://github.com/pygae/GAlgebra.jl): Julia interface to GAlgebra
  - contains test cases for [Eric Chisolm's Geometric Algebra](https://arxiv.org/abs/1205.5935)

## Explorational

- Dr. Eric Wieser's PhD thesis *Formalizing Clifford algebras and related constructions in the Lean theorem prover*
    - in section 3.3. Symbolic, the approach of GAlgebra is demonstrated and discussed
- [appliedgeometry/poissongeometry](https://github.com/appliedgeometry/poissongeometry): A Python module for (local) Poisson-Nijenhuis calculus on Poisson manifolds, with the following functions
  - used in notebook tutorials: https://github.com/appliedgeometry/poissongeometry/blob/master/docs/Tutorial_Ingles.ipynb
- [meuns/galgebra: GA4CS examples, PGA, code-generation, and more](https://github.com/meuns/galgebra)
  - see also https://github.com/pygae/galgebra/pull/68
- [hugohadfield/sym_scratch](https://github.com/hugohadfield/sym_scratch)
  - used in notebooks for CGA
