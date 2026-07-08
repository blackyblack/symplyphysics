# Agent guides

## Writing proofs (derivations) for laws

Proofs are import-time symbolic derivations placed in the law module between `law` and the
`calculate_*` function. See
`symplyphysics/thermodynamics/equations_of_state/van_der_waals/critical_van_der_waals_pressure.py`
for a reference example.

Non-obvious hints:

- Derive the law only from other laws already in the repo: import their modules and use their
  `law` attributes via `subs`/`solve`. Never restate another law's formula as a local `Eq` —
  that hides the dependency and can create undetectable circular reasoning.
- Finish the proof with `assert expr_equals(derived_expr, law.rhs)` (or `expr_equals_abs` when
  only magnitudes are comparable) from `symplyphysics.core.expr_comparisons`.
- Proofs must not contain deliberate choices: every substitution must be justified by a law or a
  stated assumption in a comment. Do not hardcode intermediate formulas. Minor hardcoding is
  allowed only when the proof would otherwise require multiple additional laws that do not exist
  in the repo yet — justify each such spot with a `TODO` comment naming the missing laws.
- Prefix all proof-local names with `_` (e.g. `_isotherm_pressure`) so they are not exported.
- Comment style: `#` for section headers of the derivation, `##` for step-by-step explanations.
- Proof dependencies are Python imports, so a dependency loop between proofs fails at import
  time. If you hit such a loop, decide which law is more fundamental, keep only the derivation
  in that direction (remove or reverse the other proof), and record the chosen direction in the
  law's docstring.
- Run `pytest`, `pylint`, `mypy` and format with `yapf` before committing.
