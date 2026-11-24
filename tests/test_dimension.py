from fractions import Fraction

import pytest

from dimcheck import Dimension, UnitError, convert, parse_dimension, parse_unit

L, T, M = Dimension.base("L"), Dimension.base("T"), Dimension.base("M")


def test_algebra():
    assert L / T / T == L * T**-2
    assert str(M * L**2 / T**2) == "M L^2 T^-2"
    assert (L ** Fraction(1, 2)) ** 2 == L
    assert Dimension().is_dimensionless
    assert str(Dimension()) == "1"


def test_unknown_base():
    with pytest.raises(ValueError):
        Dimension.base("X")


@pytest.mark.parametrize(
    "text, expected",
    [
        ("N", M * L / T**2),
        ("kg m s^-2", M * L / T**2),
        ("kg*m/s**2", M * L / T**2),
        ("J/(mol K)", M * L**2 / T**2 / Dimension.base("N") / Dimension.base("Θ")),
        ("m^(1/2)", L ** Fraction(1, 2)),
        ("Pa·s", M / L / T),
        ("Ω", M * L**2 / T**3 / Dimension.base("I") ** 2),
        ("", Dimension()),
        ("[M L T^-2]", M * L / T**2),
        ("[1]", Dimension()),
    ],
)
def test_parse(text, expected):
    assert parse_dimension(text) == expected


def test_prefixes_and_scale():
    assert parse_unit("km").scale == pytest.approx(1e3)
    assert parse_unit("µs").scale == pytest.approx(1e-6)
    assert parse_unit("mg").scale == pytest.approx(1e-6)
    assert parse_unit("dam").scale == pytest.approx(10)
    assert parse_unit("GeV").scale == pytest.approx(1.602176634e-10)


def test_convert():
    assert convert(1, "kWh", "MJ") == pytest.approx(3.6)
    assert convert(36, "km/h", "m/s") == pytest.approx(10)
    assert convert(1, "atm", "bar") == pytest.approx(1.01325)
    with pytest.raises(UnitError):
        convert(1, "m", "s")


@pytest.mark.parametrize("bad", ["furlong", "m^", "(m", "m)", "kg**x"])
def test_bad_units(bad):
    with pytest.raises(UnitError):
        parse_unit(bad)
