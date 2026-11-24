import sympy as sp
import pytest

from dimcheck import Dimension, DimensionError, check, dimension_of

MECH = {"m": "kg", "g": "m/s^2", "L": "m", "t": "s", "v": "m/s", "x": "m", "E": "J", "c": "m/s"}


def test_pendulum_period():
    assert check("T = 2*pi*sqrt(L/g)", {**MECH, "T": "s"})


def test_mass_energy_wrong_power():
    r = check("E = m*c**3", MECH)
    assert not r
    assert "left side" in r.message
    assert r.sides[0] != r.sides[1]


def test_sum_mismatch_reports_terms():
    r = check("x + v", MECH)
    assert not r
    assert "inconsistent terms of a sum" in r.message


def test_transcendental_argument():
    assert not check("exp(-t)", MECH)
    assert check("L*exp(-t*g/v)", MECH)
    assert not check("log(m)", MECH)


def test_symbolic_exponent_on_dimensional_base():
    r = check("L**m", {"L": "m", "m": ""})
    assert not r and "rational" in r.message


def test_dimensionless_symbolic_exponent_ok():
    assert check("(x/L)**n", {"x": "m", "L": "m", "n": ""})


def test_zero_is_compatible():
    assert check("x - x", MECH)  # simplifies to 0
    assert check("0 = v - v", MECH)


def test_strict_and_loose():
    with pytest.raises(DimensionError):
        dimension_of("q*L", MECH)
    assert dimension_of("q*L", MECH, strict=False) == Dimension.base("L")


def test_calculus():
    dims = {"x": "m", "t": "s", "omega": "1/s"}
    assert check("Derivative(x(t), t, 2) = -omega**2*x(t)", dims)
    assert not check("Derivative(x(t), t) = -omega**2*x(t)", dims)
    assert str(dimension_of("Integral(x, (t, 0, tau))", {**dims, "tau": "s"})) == "L T"


def test_sympy_objects_accepted():
    m, c = sp.symbols("m c")
    assert str(dimension_of(m * c**2, MECH)) == "M L^2 T^-2"


def test_quantum_example():
    dims = {"hbar": "J s", "m": "kg", "a": "m", "E": "J"}
    assert check("E = pi**2*hbar**2/(2*m*a**2)", dims)


def test_piecewise_and_minmax():
    assert check("Max(x, L)", MECH)
    assert not check("Min(x, t)", MECH)
    assert check("Piecewise((x, t > 0), (L, True))", MECH)


def test_parse_error_is_value_error():
    with pytest.raises(ValueError):
        dimension_of("sin(", MECH)
