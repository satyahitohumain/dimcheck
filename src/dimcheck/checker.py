"""Dimensional analysis of SymPy expressions and equations."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Mapping, Union

import sympy as sp

from .dimension import Dimension
from .units import parse_dimension

DimSpec = Union[str, Dimension]
ONE = Dimension.dimensionless()

# Functions whose argument must be dimensionless and whose value is dimensionless.
_TRANSCENDENTAL = (
    sp.exp, sp.log, sp.sin, sp.cos, sp.tan, sp.cot, sp.sec, sp.csc,
    sp.asin, sp.acos, sp.atan, sp.acot, sp.sinh, sp.cosh, sp.tanh, sp.coth,
    sp.asinh, sp.acosh, sp.atanh, sp.erf, sp.erfc, sp.gamma, sp.besselj, sp.bessely,
)


class DimensionError(ValueError):
    """Raised when an expression is dimensionally inconsistent.

    ``subexpr`` is the smallest sub-expression where the problem was detected.
    """

    def __init__(self, message: str, subexpr: sp.Basic | None = None):
        super().__init__(message)
        self.subexpr = subexpr


@dataclass
class CheckResult:
    """Outcome of :func:`check`."""

    ok: bool
    dimension: Dimension | None = None
    message: str = ""
    subexpr: sp.Basic | None = None
    sides: tuple[Dimension, Dimension] | None = field(default=None)

    def __bool__(self) -> bool:
        return self.ok


def _as_dim(spec: DimSpec) -> Dimension:
    return spec if isinstance(spec, Dimension) else parse_dimension(spec)


def _rational(value: sp.Basic, where: sp.Basic) -> Fraction:
    if value.is_Rational:
        return Fraction(int(value.p), int(value.q))
    if value.is_Float:
        frac = Fraction(float(value)).limit_denominator(1000)
        if abs(float(frac) - float(value)) < 1e-12:
            return frac
    raise DimensionError(f"exponent {value} of a dimensional base must be a rational number", where)


class Analyzer:
    """Computes the dimension of SymPy expressions given dimensions for their symbols.

    Parameters
    ----------
    dims:
        Mapping from symbol (or symbol name) to a unit string (``"m/s"``), a bracketed
        dimension (``"[L T^-1]"``) or a :class:`Dimension`.
    strict:
        If true (default), symbols without a declared dimension raise an error.
        Otherwise they are treated as dimensionless.
    """

    def __init__(self, dims: Mapping[Union[str, sp.Symbol], DimSpec], strict: bool = True):
        self.dims: dict[str, Dimension] = {str(k): _as_dim(v) for k, v in dims.items()}
        self.strict = strict

    def dimension(self, expr: sp.Basic) -> Dimension:
        return self._dim(sp.sympify(expr))

    def _dim(self, e: sp.Basic) -> Dimension:
        if e.is_Number or e in (sp.pi, sp.E, sp.I, sp.EulerGamma, sp.oo, -sp.oo):
            return ONE
        if e.is_Symbol:
            name = str(e)
            if name in self.dims:
                return self.dims[name]
            if self.strict:
                raise DimensionError(f"no dimension declared for symbol {name!r}", e)
            return ONE
        if isinstance(e, sp.Function) and not isinstance(e, _TRANSCENDENTAL) and isinstance(e.func, sp.core.function.UndefinedFunction):
            name = str(e.func)
            if name in self.dims:
                return self.dims[name]
            raise DimensionError(f"no dimension declared for function {name!r}", e)
        if e.is_Add:
            return self._same(e.args, e, "terms of a sum")
        if e.is_Mul:
            out = ONE
            for a in e.args:
                out = out * self._dim(a)
            return out
        if e.is_Pow:
            base, exp = e.args
            bdim = self._dim(base)
            edim = self._dim(exp)
            if not edim.is_dimensionless:
                raise DimensionError(f"exponent {exp} has dimension [{edim}], must be dimensionless", e)
            if bdim.is_dimensionless:
                return ONE
            return bdim ** _rational(exp, e)
        if isinstance(e, _TRANSCENDENTAL):
            for a in e.args:
                adim = self._dim(a)
                if not adim.is_dimensionless:
                    raise DimensionError(
                        f"argument of {e.func.__name__}() has dimension [{adim}], must be dimensionless", e
                    )
            return ONE
        if isinstance(e, sp.Abs) or isinstance(e, (sp.re, sp.im, sp.conjugate, sp.sign)):
            inner = self._dim(e.args[0])
            return ONE if isinstance(e, sp.sign) else inner
        if isinstance(e, (sp.Min, sp.Max)):
            return self._same(e.args, e, f"arguments of {e.func.__name__}()")
        if isinstance(e, sp.Piecewise):
            return self._same([branch for branch, _ in e.args], e, "branches of a Piecewise")
        if isinstance(e, sp.Derivative):
            out = self._dim(e.expr)
            for var, count in e.variable_count:
                out = out / self._dim(var) ** int(count)
            return out
        if isinstance(e, sp.Integral):
            out = self._dim(e.function)
            for limit in e.limits:
                var = limit[0]
                vdim = self._dim(var)
                for bound in limit[1:]:
                    if bound != 0 and self._dim(bound) != vdim:
                        raise DimensionError(f"integration bound {bound} does not match variable {var}", e)
                out = out * vdim
            return out
        if isinstance(e, sp.Equality) or isinstance(e, sp.core.relational.Relational):
            return self._same(e.args, e, "sides of a relation")
        raise DimensionError(f"don't know how to find the dimension of {e.func.__name__}", e)

    def _same(self, args, where: sp.Basic, what: str) -> Dimension:
        dims = []
        for a in args:
            if a == 0:  # zero is compatible with any dimension
                continue
            dims.append((a, self._dim(a)))
        if not dims:
            return ONE
        first = dims[0][1]
        for a, d in dims[1:]:
            if d != first:
                detail = ", ".join(f"{t} ~ [{dd}]" for t, dd in dims)
                raise DimensionError(f"inconsistent {what}: {detail}", where)
        return first


def _parse(text_or_expr: Union[str, sp.Basic], names) -> sp.Basic:
    try:
        return _parse_inner(text_or_expr, names)
    except Exception as exc:  # tokenizer and sympify errors vary by Python and SymPy version
        if isinstance(exc, DimensionError):
            raise
        raise ValueError(f"could not parse {text_or_expr!r}: {exc}") from exc


def _parse_inner(text_or_expr: Union[str, sp.Basic], names) -> sp.Basic:
    if not isinstance(text_or_expr, str):
        return sp.sympify(text_or_expr)
    text = text_or_expr
    local: dict[str, sp.Basic] = {}
    for n in map(str, names):
        # A declared name used with call syntax, e.g. x(t), is an undefined function.
        called = re.search(rf"(?<![\w.]){re.escape(n)}\s*\(", text)
        local[n] = sp.Function(n) if called else sp.Symbol(n)
    if "==" in text:
        lhs, rhs = text.split("==", 1)
        return sp.Eq(sp.parse_expr(lhs, local_dict=local), sp.parse_expr(rhs, local_dict=local), evaluate=False)
    if text.count("=") == 1 and not any(op in text for op in ("<=", ">=", "!=")):
        lhs, rhs = text.split("=")
        return sp.Eq(sp.parse_expr(lhs, local_dict=local), sp.parse_expr(rhs, local_dict=local), evaluate=False)
    return sp.parse_expr(text, local_dict=local, evaluate=False)


def dimension_of(expr: Union[str, sp.Basic], dims: Mapping[str, DimSpec], strict: bool = True) -> Dimension:
    """Return the dimension of ``expr``; raise :class:`DimensionError` if it is inconsistent.

    >>> str(dimension_of("sqrt(g/L)", {"g": "m/s^2", "L": "m"}))
    'T^-1'
    """
    return Analyzer(dims, strict=strict).dimension(_parse(expr, dims.keys()))


def check(expr: Union[str, sp.Basic], dims: Mapping[str, DimSpec], strict: bool = True) -> CheckResult:
    """Check an expression or equation (``"E = m*c**2"``) for dimensional consistency.

    Never raises for dimensional problems; inspect the returned :class:`CheckResult`.
    """
    analyzer = Analyzer(dims, strict=strict)
    e = _parse(expr, dims.keys())
    try:
        if isinstance(e, sp.core.relational.Relational):
            left, right = analyzer.dimension(e.lhs), analyzer.dimension(e.rhs)
            if left != right and e.lhs != 0 and e.rhs != 0:
                return CheckResult(
                    False, None,
                    f"left side has dimension [{left}] but right side has [{right}]",
                    e, (left, right),
                )
            return CheckResult(True, left if e.lhs != 0 else right, "consistent", None, (left, right))
        return CheckResult(True, analyzer.dimension(e), "consistent")
    except DimensionError as exc:
        return CheckResult(False, None, str(exc), exc.subexpr)
