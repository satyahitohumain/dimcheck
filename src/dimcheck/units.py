"""Unit table and a small parser for unit strings such as ``"kg*m/s^2"`` or ``"J/(mol K)"``."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from fractions import Fraction

from .dimension import Dimension

M, L, T, I, K, N, J = (Dimension.base(b) for b in ("M", "L", "T", "I", "Θ", "N", "J"))
ONE = Dimension.dimensionless()


@dataclass(frozen=True)
class Unit:
    """A unit: a scale factor relative to SI and a dimension."""

    scale: float
    dim: Dimension

    def __mul__(self, other: "Unit") -> "Unit":
        return Unit(self.scale * other.scale, self.dim * other.dim)

    def __truediv__(self, other: "Unit") -> "Unit":
        return Unit(self.scale / other.scale, self.dim / other.dim)

    def __pow__(self, p: Fraction) -> "Unit":
        return Unit(self.scale ** float(p), self.dim ** p)


def _u(scale: float, dim: Dimension) -> Unit:
    return Unit(scale, dim)


_NEWTON = M * L / T**2
_JOULE = _NEWTON * L
_COULOMB = I * T
_VOLT = _JOULE / _COULOMB

# Units that accept SI prefixes.
PREFIXABLE: dict[str, Unit] = {
    "m": _u(1, L),
    "g": _u(1e-3, M),
    "s": _u(1, T),
    "A": _u(1, I),
    "K": _u(1, K),
    "mol": _u(1, N),
    "cd": _u(1, J),
    "N": _u(1, _NEWTON),
    "J": _u(1, _JOULE),
    "W": _u(1, _JOULE / T),
    "Pa": _u(1, _NEWTON / L**2),
    "C": _u(1, _COULOMB),
    "V": _u(1, _VOLT),
    "ohm": _u(1, _VOLT / I),
    "Ω": _u(1, _VOLT / I),
    "S": _u(1, I / _VOLT),
    "F": _u(1, _COULOMB / _VOLT),
    "H": _u(1, _VOLT * T / I),
    "T": _u(1, _VOLT * T / L**2),
    "Wb": _u(1, _VOLT * T),
    "Hz": _u(1, ONE / T),
    "eV": _u(1.602176634e-19, _JOULE),
    "L": _u(1e-3, L**3),
    "bar": _u(1e5, _NEWTON / L**2),
    "Wh": _u(3600, _JOULE),
}

# Units that do not take prefixes.
PLAIN: dict[str, Unit] = {
    "min": _u(60, T),
    "h": _u(3600, T),
    "hr": _u(3600, T),
    "day": _u(86400, T),
    "yr": _u(3.15576e7, T),  # Julian year
    "au": _u(1.495978707e11, L),
    "ly": _u(9.4607304725808e15, L),
    "pc": _u(3.0856775814913673e16, L),
    "atm": _u(101325, _NEWTON / L**2),
    "cal": _u(4.184, _JOULE),
    "u": _u(1.66053906660e-27, M),
    "Da": _u(1.66053906660e-27, M),
    "rad": _u(1, ONE),
    "sr": _u(1, ONE),
    "deg": _u(math.pi / 180, ONE),
    "percent": _u(0.01, ONE),
}

PREFIXES: dict[str, float] = {
    "Q": 1e30, "R": 1e27, "Y": 1e24, "Z": 1e21, "E": 1e18, "P": 1e15, "T": 1e12,
    "G": 1e9, "M": 1e6, "k": 1e3, "h": 1e2, "da": 1e1,
    "d": 1e-1, "c": 1e-2, "m": 1e-3, "u": 1e-6, "µ": 1e-6, "μ": 1e-6, "n": 1e-9,
    "p": 1e-12, "f": 1e-15, "a": 1e-18, "z": 1e-21, "y": 1e-24, "r": 1e-27, "q": 1e-30,
}


class UnitError(ValueError):
    """Raised when a unit string cannot be parsed or a unit is unknown."""


def lookup(name: str) -> Unit:
    """Resolve a single unit name, applying an SI prefix if needed (``"km"``, ``"µs"``)."""
    if name in PLAIN:
        return PLAIN[name]
    if name in PREFIXABLE:
        return PREFIXABLE[name]
    # Longest prefix first so that "da" wins over "d".
    for prefix in sorted(PREFIXES, key=len, reverse=True):
        if name.startswith(prefix) and name[len(prefix):] in PREFIXABLE:
            base = PREFIXABLE[name[len(prefix):]]
            return Unit(base.scale * PREFIXES[prefix], base.dim)
    raise UnitError(f"unknown unit {name!r}")


_TOKEN = re.compile(
    r"\s*(?:(?P<num>\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)|(?P<name>[A-Za-zΩµμ_]+)|(?P<op>\*\*|[*/^()·⋅.\-+]))"
)


def _tokenize(text: str) -> list[tuple[str, str]]:
    tokens: list[tuple[str, str]] = []
    pos = 0
    text = text.strip()
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if not m or m.end() == pos:
            raise UnitError(f"cannot parse unit string {text!r} at position {pos}")
        kind = m.lastgroup
        assert kind is not None
        value = m.group(kind)
        if kind == "op" and value in ("·", "⋅", "."):
            value = "*"
        if kind == "op" and value == "**":
            value = "^"
        tokens.append((kind, value))
        pos = m.end()
    return tokens


class _Parser:
    def __init__(self, text: str):
        self.text = text
        self.tokens = _tokenize(text)
        self.i = 0

    def peek(self) -> tuple[str, str] | None:
        return self.tokens[self.i] if self.i < len(self.tokens) else None

    def take(self) -> tuple[str, str]:
        tok = self.peek()
        if tok is None:
            raise UnitError(f"unexpected end of unit string {self.text!r}")
        self.i += 1
        return tok

    def parse(self) -> Unit:
        if not self.tokens:
            return Unit(1.0, ONE)
        unit = self.expr()
        if self.peek() is not None:
            raise UnitError(f"unexpected {self.peek()[1]!r} in unit string {self.text!r}")
        return unit

    def expr(self) -> Unit:
        unit = self.power()
        while True:
            tok = self.peek()
            if tok == ("op", "*"):
                self.take()
                unit = unit * self.power()
            elif tok == ("op", "/"):
                self.take()
                unit = unit / self.power()
            elif tok is not None and (tok[0] in ("name", "num") or tok == ("op", "(")):
                unit = unit * self.power()  # implicit multiplication: "N m"
            else:
                return unit

    def power(self) -> Unit:
        unit = self.atom()
        tok = self.peek()
        if tok == ("op", "^"):
            self.take()
            return unit ** self.exponent()
        return unit

    def exponent(self) -> Fraction:
        sign = 1
        tok = self.peek()
        if tok == ("op", "("):
            self.take()
            value = self.exponent()
            if self.take() != ("op", ")"):
                raise UnitError(f"expected ')' in exponent of {self.text!r}")
            return value
        if tok in (("op", "-"), ("op", "+")):
            self.take()
            sign = -1 if tok[1] == "-" else 1
        kind, value = self.take()
        if kind != "num":
            raise UnitError(f"expected a number as exponent in {self.text!r}, got {value!r}")
        result = Fraction(value)
        if self.peek() == ("op", "/") and self.i + 1 < len(self.tokens) and self.tokens[self.i + 1][0] == "num":
            self.take()
            result /= Fraction(self.take()[1])
        return sign * result

    def atom(self) -> Unit:
        kind, value = self.take()
        if kind == "name":
            return lookup(value)
        if kind == "num":
            return Unit(float(value), ONE)
        if value == "(":
            unit = self.expr()
            if self.take() != ("op", ")"):
                raise UnitError(f"expected ')' in {self.text!r}")
            return unit
        raise UnitError(f"unexpected {value!r} in unit string {self.text!r}")


def parse_unit(text: str) -> Unit:
    """Parse a unit string.

    Supports ``*``, ``/``, ``·``, implicit multiplication by whitespace, ``^`` or ``**``
    with integer or rational exponents (``m^(1/2)``), parentheses and SI prefixes.

    >>> str(parse_unit("kg m^2 s^-2").dim)
    'M L^2 T^-2'
    >>> parse_unit("km/h").scale
    0.2777777777777778
    """
    return _Parser(text).parse()


def parse_dimension(text: str) -> Dimension:
    """Parse either a unit string or a bracketed dimension string like ``"[M L T^-2]"``."""
    text = text.strip()
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        dim = ONE
        if inner in ("", "1"):
            return dim
        for part in inner.replace("*", " ").split():
            base, _, exp = part.partition("^")
            dim = dim * Dimension.base(base) ** Fraction(exp.strip("()") or 1)
        return dim
    return parse_unit(text).dim
