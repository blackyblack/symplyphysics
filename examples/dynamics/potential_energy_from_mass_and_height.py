#!/usr/bin/env python3

from sympy import solve, Eq
from symplyphysics import print_expression, units, Symbol
from symplyphysics.laws.dynamics import potential_energy_from_mass_and_height as potential_energy

# Under what condition will the potential energies of two different bodies be the same
# if the mass of the first body is half the mass of the second body?

height_1 = Symbol("height_1", units.length)
height_2 = Symbol("height_2", units.length)
body_mass_1 = Symbol("body_mass_1", units.mass)

Ep1 = potential_energy.law.subs({
    potential_energy.mass: body_mass_1,
    potential_energy.height: height_1
})
Ep2 = potential_energy.law.subs({
    potential_energy.mass: 2 * body_mass_1,
    potential_energy.height: height_2
})
law = [Ep2, Ep1]
solved = solve(law, (height_1, potential_energy.potential_energy), dict=True)[0][height_1]
answer = Eq(height_1, solved)

print(f"\nFormula is:\n{print_expression(potential_energy.law)}")
print(
    f"\nSolution:\nIF potential_energy_1 = potential_energy_2 AND body_mass_2 = 2 * body_mass_1 THEN {print_expression(answer)}"
)
