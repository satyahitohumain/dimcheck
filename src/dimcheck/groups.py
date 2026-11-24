"""Buckingham-Pi groups and Rayleigh's method of dimensional analysis."""

from __future__ import annotations

from fractions import Fraction
from math import lcm
from typing import Mapping

import sympy as sp

from .checker import DimSpec, _as_dim
from .dimension import BASES, Dimension


def _matrix(dims: Mapping[str, Dimension]) -> sp.Matrix:
    cols = [[sp.Rational(e.numerator, e.denominator) for e in d.exponents] for d in dims.values()]
    return sp.Matrix(cols).T  # rows: base dimensions, columns: variables


def _integerize(vec: sp.Matrix) -> list[int]:
    fracs = [Fraction(int(sp.fraction(v)[0]), int(sp.fraction(v)[1])) for v in vec]
    scale = lcm(*(f.denominator for f in fracs)) if fracs else 1
    ints = [int(f * scale) for f in fracs]
    # Make the first nonzero exponent positive for a canonical form.
    for x in ints:
        if x != 0:
            if x < 0:
                ints = [-y for y in ints]
            break
    return ints


def _monomial(names: list[str], exps: list[int]) -> sp.Expr:
    expr = sp.Integer(1)
    for name, e in zip(names, exps):
        expr *= sp.Symbol(name, positive=True) ** e
    return expr


def pi_groups(variables: Mapping[str, DimSpec]) -> list[sp.Expr]:
    """Return a basis of independent dimensionless products of ``variables``.

    By the Buckingham-Pi theorem there are ``n - r`` of them, where ``n`` is the number of
    variables and ``r`` the rank of their dimension matrix.

    >>> pi_groups({"rho": "kg/m^3", "v": "m/s", "D": "m", "mu": "Pa*s"})
    [D*rho*v/mu]
    """
    names = list(variables)
    dims = {n: _as_dim(v) for n, v in variables.items()}
    basis = _matrix(dims).nullspace()
    groups = [_monomial(names, _integerize(v)) for v in basis]
    return groups


def rayleigh(target: DimSpec, variables: Mapping[str, DimSpec]) -> tuple[sp.Expr, list[sp.Expr]]:
    """Find a product of powers of ``variables`` with the dimension ``target``.

    Returns ``(monomial, free_groups)``: the target quantity is ``monomial`` times an arbitrary
    function of the dimensionless ``free_groups``. Raises ``ValueError`` if no product of the
    variables can have the target dimension.

    >>> rayleigh("s", {"L": "m", "g": "m/s^2", "m": "kg"})
    (sqrt(L)/sqrt(g), [])
    """
    names = list(variables)
    dims = {n: _as_dim(v) for n, v in variables.items()}
    A = _matrix(dims)
    t = _as_dim(target)
    b = sp.Matrix([sp.Rational(e.numerator, e.denominator) for e in t.exponents])
    xs = sp.symbols(f"x0:{len(names)}")
    solutions = sp.linsolve((A, b), *xs)
    if not solutions:
        missing = [BASES[i] for i in range(len(BASES)) if t.exponents[i] != 0]
        raise ValueError(f"no product of {names} has dimension [{t}] (check {missing})")
    (sol,) = solutions
    particular = [s.subs({x: 0 for x in xs}) for s in sol]
    expr = sp.Integer(1)
    for name, e in zip(names, particular):
        expr *= sp.Symbol(name, positive=True) ** e
    return expr, pi_groups(variables)
