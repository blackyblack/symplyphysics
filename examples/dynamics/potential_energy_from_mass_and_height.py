from symplyphysics import (
    symbols, Eq, pretty, solve
)
from symplyphysics.laws.dynamics import potential_energy_from_mass_and_height as potential_energy

# Under what condition will the potential energies of two different bodies be the same
# if the mass of the first body is half the mass of the second body?

body_mass_1, height_1, height_2 = symbols("body_mass_1 height_1 height_2")

Ep1 = potential_energy.law.subs({potential_energy.body_mass: body_mass_1, potential_energy.height: height_1})
Ep2 = potential_energy.law.subs({potential_energy.body_mass: 2 * body_mass_1, potential_energy.height: height_2})
law = [Ep2, Ep1]
solved = solve(law, (height_1, potential_energy.potential_energy_of_body), dict=True)[0][height_1]
answer = Eq(height_1, solved)

print("\nFormula is:\n{}".format(potential_energy.print(potential_energy.law)))
print("\nSolution:\nIF potential_energy_1 = potential_energy_2 AND body_mass_2 = 2 * body_mass_1 THEN {}".
      format(pretty(answer, use_unicode=False)))