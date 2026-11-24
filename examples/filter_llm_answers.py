"""Reject dimensionally impossible final answers before spending money on an LLM judge.

Many wrong answers to physics problems fail a units check. Running this filter first is a cheap,
deterministic way to catch them, and it gives a precise reason that can be shown to a model
as feedback.

    python examples/filter_llm_answers.py
"""

from dimcheck import check

# Problem: small-oscillation angular frequency of a bead on a hoop of radius R spinning at
# angular speed w about a vertical diameter (for w**2 > g/R).
units = {"Omega": "1/s", "w": "1/s", "g": "m/s^2", "R": "m"}

candidates = {
    "model_a": "Omega = sqrt(w**2 - g**2/(R**2*w**2))",  # correct
    "model_b": "Omega = sqrt(w**2 - g/(R*w**2))",        # dropped a power
    "model_c": "Omega = w*sqrt(1 - g/R)",                # g/R is not dimensionless
    "model_d": "Omega = sqrt(g/R)*sin(w)",               # sin of a frequency
}

for name, answer in candidates.items():
    result = check(answer, units)
    status = "pass" if result else "REJECT"
    print(f"{name:8s} {status:6s} {answer}")
    if not result:
        print(f"         reason: {result.message}")
