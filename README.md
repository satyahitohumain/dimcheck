# dimcheck

[![tests](https://github.com/satyahitohumain/dimcheck/actions/workflows/tests.yml/badge.svg)](https://github.com/satyahitohumain/dimcheck/actions/workflows/tests.yml)

Dimensional analysis for SymPy expressions. `dimcheck` tells you whether a formula or an
equation is dimensionally consistent and, if it isn't, points at the exact sub-expression that
breaks. It also finds Buckingham-Pi groups, derives power-law forms by Rayleigh's method,
and converts between units.

I built it as a fast, deterministic first filter for grading physics answers written by
people or by language models. A final answer with the wrong dimensions is wrong, and
proving that takes microseconds instead of a call to an LLM judge.

```console
$ dimcheck check "Omega = w*sqrt(1 - g/R)" -u Omega=1/s -u w=1/s -u g=m/s^2 -u R=m
FAIL  inconsistent terms of a sum: 1 ~ [1], -g/R ~ [T^-2]
      at: 1 - g/R
```

## Install

```bash
pip install git+https://github.com/satyahitohumain/dimcheck
# or, from a clone
pip install -e ".[dev]"
```

Requires Python 3.10+ and SymPy.

## Command line

```console
$ dimcheck check "T = 2*pi*sqrt(L/g)" -u T=s -u L=m -u g=m/s^2
OK  [T]

$ dimcheck check "Derivative(x(t), t, 2) = -omega**2*x(t)" -u x=m -u t=s -u omega=1/s
OK  [L T^-2]

$ dimcheck dim "hbar**2/(m*a**2)" -u hbar="J s" -u m=kg -u a=m
M L^2 T^-2

$ dimcheck pi -u F=N -u rho=kg/m^3 -u v=m/s -u D=m -u mu="Pa*s"
Pi_1 = F/(D**2*rho*v**2)
Pi_2 = F*rho/mu**2

$ dimcheck solve s -u L=m -u g=m/s^2 -u m=kg
[s] ~ sqrt(L)/sqrt(g)

$ dimcheck convert 1 kWh MJ
3.6 MJ
```

Exit codes: `0` consistent, `1` inconsistent, `2` parse or unit error, so `dimcheck check`
drops straight into shell scripts and CI.

## Python API

```python
from dimcheck import check, dimension_of, pi_groups, rayleigh, convert

check("E = m*c**2", {"E": "J", "m": "kg", "c": "m/s"})        # CheckResult(ok=True, ...)
r = check("x = A*exp(-t)", {"x": "m", "A": "m", "t": "s"})
r.ok, r.message, r.subexpr   # (False, 'argument of exp() has dimension [T], ...', exp(-t))

dimension_of("sqrt(g/L)", {"g": "m/s^2", "L": "m"})          # Dimension('T^-1')
pi_groups({"rho": "kg/m^3", "v": "m/s", "D": "m", "mu": "Pa*s"})  # [D*rho*v/mu]
convert(36, "km/h", "m/s")                                   # 10.0
```

Dimensions can be given as unit strings (`"kg m^2 s^-2"`, `"J/(mol K)"`, `"m^(1/2)"`,
`"µs"`, `"Ω"`) or as bracketed base dimensions (`"[M L T^-2]"`). `check` also accepts SymPy
objects directly.

## What it checks

| Construct | Rule |
|---|---|
| `a + b`, `Min`, `Max`, `Piecewise`, `=`, `<` | every term has the same dimension (a literal `0` matches anything) |
| `a * b`, `a / b` | dimensions multiply |
| `a ** p` | `p` is dimensionless; if `a` has dimensions, `p` must be a rational number |
| `exp`, `log`, trig, hyperbolic, `erf`, `gamma`, Bessel | argument dimensionless, result dimensionless |
| `Derivative(f, x, n)` | `[f] / [x]^n` |
| `Integral(f, (x, a, b))` | `[f] [x]`; bounds must match `[x]` |
| `x(t)` (declared name used as a function) | has the declared dimension |

Symbols without a declared dimension raise an error by default; pass `strict=False` (or
`--loose`) to treat them as dimensionless.

## Using it to grade LLM physics answers

[`examples/filter_llm_answers.py`](examples/filter_llm_answers.py) runs four candidate answers
to a rotating-hoop problem through the checker. The three dimensionally wrong ones are
rejected with a reason that can be fed back to the model:

```text
model_a  pass   Omega = sqrt(w**2 - g**2/(R**2*w**2))
model_b  REJECT Omega = sqrt(w**2 - g/(R*w**2))
         reason: inconsistent terms of a sum: w**2 ~ [T^-2], -g/(R*w**2) ~ [1]
model_c  REJECT Omega = w*sqrt(1 - g/R)
         reason: inconsistent terms of a sum: 1 ~ [1], -g/R ~ [T^-2]
model_d  REJECT Omega = sqrt(g/R)*sin(w)
         reason: argument of sin() has dimension [T^-1], must be dimensionless
```

A dimensional check is necessary, not sufficient: `model_a` passing does not prove it is right.
Pair it with a symbolic or numeric equivalence check against a reference answer.

## Design notes

- Dimensions are exponent vectors over the seven SI base quantities with `Fraction` entries,
  so half-integer powers such as `m^(1/2)` stay exact.
- Pi groups are an integer basis of the null space of the dimension matrix, scaled to the
  smallest integers with the first exponent positive.
- The unit parser is a small recursive-descent parser, not `eval`, so unit strings are safe to
  take from untrusted input. Expressions are parsed with SymPy's `parse_expr`, which does
  evaluate Python, so don't pass untrusted *expressions* to a server-side checker without a
  sandbox.

## Development

```bash
pip install -e ".[dev]"
pytest            # unit tests plus doctests
```

## References

- E. Buckingham, "On physically similar systems; illustrations of the use of dimensional
  equations", *Phys. Rev.* 4, 345 (1914).
- Bureau International des Poids et Mesures, *The International System of Units (SI
  Brochure)*, 9th ed. (2019), for base units, derived units and prefixes.

## License

MIT © 2026 Satish, Siyaramji Tech Solutions
