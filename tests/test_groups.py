import sympy as sp

from dimcheck import pi_groups, rayleigh
from dimcheck.checker import dimension_of


def test_reynolds_number():
    groups = pi_groups({"rho": "kg/m^3", "v": "m/s", "D": "m", "mu": "Pa*s"})
    assert len(groups) == 1
    rho, v, D, mu = sp.symbols("rho v D mu", positive=True)
    g = groups[0]
    assert sp.simplify(g - rho * v * D / mu) == 0 or sp.simplify(g - mu / (rho * v * D)) == 0


def test_groups_are_dimensionless_and_counted():
    variables = {"F": "N", "rho": "kg/m^3", "v": "m/s", "D": "m", "mu": "Pa*s"}
    groups = pi_groups(variables)
    assert len(groups) == 5 - 3
    for g in groups:
        assert dimension_of(g, variables).is_dimensionless


def test_pendulum_rayleigh():
    mono, free = rayleigh("s", {"L": "m", "g": "m/s^2", "m": "kg"})
    L, g = sp.symbols("L g", positive=True)
    assert sp.simplify(mono - sp.sqrt(L / g)) == 0
    assert free == []


def test_deep_water_waves():
    # Phase speed of waves of wavelength lam on deep water: v ~ sqrt(g*lam)
    mono, free = rayleigh("m/s", {"g": "m/s^2", "lam": "m"})
    assert sp.simplify(mono - sp.sqrt(sp.Symbol("g", positive=True) * sp.Symbol("lam", positive=True))) == 0


def test_rayleigh_impossible():
    import pytest

    with pytest.raises(ValueError):
        rayleigh("kg", {"L": "m", "t": "s"})
