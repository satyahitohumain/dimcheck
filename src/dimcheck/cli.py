"""Command-line interface: ``dimcheck check|dim|pi|solve|convert``."""

from __future__ import annotations

import argparse
import sys

from . import convert
from .checker import DimensionError, check, dimension_of
from .groups import pi_groups, rayleigh
from .units import UnitError


def _units(pairs: list[str]) -> dict[str, str]:
    out = {}
    for pair in pairs:
        name, sep, unit = pair.partition("=")
        if not sep:
            raise SystemExit(f"error: expected NAME=UNIT, got {pair!r}")
        out[name.strip()] = unit.strip()
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="dimcheck", description="Dimensional analysis from the command line.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    def with_units(p):
        p.add_argument("-u", "--unit", action="append", default=[], metavar="NAME=UNIT",
                       help="declare a symbol's unit, e.g. -u g=m/s^2 (repeatable)")
        p.add_argument("--loose", action="store_true", help="treat undeclared symbols as dimensionless")
        return p

    p = with_units(sub.add_parser("check", help="check an expression or equation for consistency"))
    p.add_argument("expr")
    p = with_units(sub.add_parser("dim", help="print the dimension of an expression"))
    p.add_argument("expr")
    p = with_units(sub.add_parser("pi", help="list independent dimensionless groups of the declared variables"))
    p = with_units(sub.add_parser("solve", help="find the power-law form of a quantity (Rayleigh's method)"))
    p.add_argument("target", help="unit of the target quantity, e.g. s")
    p = sub.add_parser("convert", help="convert a value between units")
    p.add_argument("value", type=float)
    p.add_argument("from_unit")
    p.add_argument("to_unit")

    args = parser.parse_args(argv)
    try:
        if args.cmd == "convert":
            print(f"{convert(args.value, args.from_unit, args.to_unit):.10g} {args.to_unit}")
            return 0
        units = _units(args.unit)
        if args.cmd == "check":
            result = check(args.expr, units, strict=not args.loose)
            if result.ok:
                print(f"OK  [{result.dimension}]")
                return 0
            print(f"FAIL  {result.message}")
            if result.subexpr is not None:
                print(f"      at: {result.subexpr}")
            return 1
        if args.cmd == "dim":
            print(dimension_of(args.expr, units, strict=not args.loose))
            return 0
        if args.cmd == "pi":
            groups = pi_groups(units)
            if not groups:
                print("no dimensionless groups")
            for i, g in enumerate(groups, 1):
                print(f"Pi_{i} = {g}")
            return 0
        if args.cmd == "solve":
            mono, groups = rayleigh(args.target, units)
            rest = f" * f({', '.join(map(str, groups))})" if groups else ""
            print(f"[{args.target}] ~ {mono}{rest}")
            return 0
    except (DimensionError, UnitError, ValueError, TypeError, SyntaxError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
