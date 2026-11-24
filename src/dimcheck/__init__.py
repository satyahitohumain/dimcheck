"""dimcheck: dimensional analysis for SymPy expressions, equations and LLM-generated answers."""

from .checker import Analyzer, CheckResult, DimensionError, check, dimension_of
from .dimension import Dimension
from .groups import pi_groups, rayleigh
from .units import Unit, UnitError, parse_dimension, parse_unit

__all__ = [
    "Analyzer", "CheckResult", "Dimension", "DimensionError", "Unit", "UnitError",
    "check", "dimension_of", "parse_dimension", "parse_unit", "pi_groups", "rayleigh",
    "convert",
]
__version__ = "0.1.0"


def convert(value: float, from_unit: str, to_unit: str) -> float:
    """Convert ``value`` between two units of the same dimension.

    >>> round(convert(1, "kWh", "MJ"), 6)
    3.6
    """
    src, dst = parse_unit(from_unit), parse_unit(to_unit)
    if src.dim != dst.dim:
        raise UnitError(f"cannot convert [{src.dim}] to [{dst.dim}]")
    return value * src.scale / dst.scale
