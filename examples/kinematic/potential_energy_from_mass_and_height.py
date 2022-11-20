from symplyphysics import (
    units, convert_to,solve, SI
)
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity
)
from symplyphysics.laws.kinematic import potential_energy_from_mass_and_height as potential_energy

# Under what condition will the potential energies of two different bodies be the same
# if the mass of the first body is half the mass of the second body?
# h1 - height body 1
# h2 - height body 2
# m1 - mass body 1
# m2 - mass body 2
# m2 = 2 * m1
# Ep1 = Ep2
m1, m2, h1, h2 = symbols('m1 m2 h1 h2')
Ep1 = potential_energy.law.subs(potential_energy.body_mass, m1).subs(potential_energy.height,h1)
Ep2 = potential_energy.law.subs(potential_energy.body_mass, 2 * m1).subs(potential_energy.height,h2)
law = Eq(Ep2, Ep1)
solved = solve(law, h1, dict=True)[0][h1]
answer = Eq(h1,solved)
print("\nFormula is:\n{}".format(potential_energy.print()))
print("\n IF potential_energy_1 = potential_energy_2 AND body_mass_2 = 2 * body_mass_1 THEN ",
          pretty(answer, use_unicode=False))