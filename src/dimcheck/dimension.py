"""Physical dimensions as exponent vectors over the seven SI base quantities."""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Union

Number = Union[int, Fraction]

BASES: tuple[str, ...] = ("M", "L", "T", "I", "Θ", "N", "J")
"""Mass, length, time, electric current, temperature, amount, luminous intensity."""

BASE_NAMES = {
    "M": "mass",
    "L": "length",
    "T": "time",
    "I": "current",
    "Θ": "temperature",
    "N": "amount",
    "J": "luminous intensity",
}


class Dimension:
    """An immutable product of base dimensions with rational exponents.

    >>> L, T = Dimension.base("L"), Dimension.base("T")
    >>> str(L / T**2)
    'L T^-2'
    """

    __slots__ = ("_exps",)

    def __init__(self, exps: Iterable[Number] = ()):
        values = tuple(Fraction(e) for e in exps)
        if len(values) > len(BASES):
            raise ValueError(f"at most {len(BASES)} exponents, got {len(values)}")
        self._exps: tuple[Fraction, ...] = values + (Fraction(0),) * (len(BASES) - len(values))

    @classmethod
    def base(cls, symbol: str) -> "Dimension":
        if symbol == "K":  # common ASCII stand-in for Θ
            symbol = "Θ"
        try:
            i = BASES.index(symbol)
        except ValueError:
            raise ValueError(f"unknown base dimension {symbol!r}; expected one of {BASES}") from None
        exps = [0] * len(BASES)
        exps[i] = 1
        return cls(exps)

    @classmethod
    def dimensionless(cls) -> "Dimension":
        return cls()

    @property
    def exponents(self) -> tuple[Fraction, ...]:
        return self._exps

    @property
    def is_dimensionless(self) -> bool:
        return not any(self._exps)

    def __mul__(self, other: "Dimension") -> "Dimension":
        return Dimension(a + b for a, b in zip(self._exps, other._exps))

    def __truediv__(self, other: "Dimension") -> "Dimension":
        return Dimension(a - b for a, b in zip(self._exps, other._exps))

    def __pow__(self, power: Number) -> "Dimension":
        p = Fraction(power)
        return Dimension(a * p for a in self._exps)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Dimension) and self._exps == other._exps

    def __hash__(self) -> int:
        return hash(self._exps)

    def __str__(self) -> str:
        if self.is_dimensionless:
            return "1"
        parts = []
        for name, e in zip(BASES, self._exps):
            if e == 0:
                continue
            parts.append(name if e == 1 else f"{name}^{e}")
        return " ".join(parts)

    def __repr__(self) -> str:
        return f"Dimension({str(self)!r})"
