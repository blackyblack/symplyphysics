from symplyphysics import (
    units, convert_to, SI
)
from symplyphysics.laws.kinematic import potential_energy_from_mass_and_height as potential_energy

#A bullet weighing 9 grams is fired from the gun vertically upwards and reaches
#an altitude of 1,500 meters. Under the action of gravity it falls freely to the ground.
#The potential energy of the fall is sufficient to cause serious damage to a
#person who might be hit by the bullet.
#
#https://ru.wikipedia.org/wiki/%D0%94%D1%83%D0%BB%D1%8C%D0%BD%D0%B0%D1%8F_%D1%8D%D0%BD%
#D0%B5%D1%80%D0%B3%D0%B8%D1%8F

m = units.Quantity('m')
SI.set_quantity_dimension(m, units.mass)
SI.set_quantity_scale_factor(m, 0.009 * units.kilogram)
g = units.Quantity('g')
SI.set_quantity_dimension(g, units.acceleration)
SI.set_quantity_scale_factor(g, 9.82 * units.meter / units.second**2)
h = units.Quantity('h')
SI.set_quantity_dimension(h, units.length)
SI.set_quantity_scale_factor(h, 1500 * units.meter)

print("\nFormula is:\n{}".format(potential_energy.print()))
result = potential_energy.calculate_potential_energy(m, g, h)
print("\nPotential energy = {} {}; for mass = {} {}, g = {} {} height = {} {}\n"
      .format(
          convert_to(result, units.joule).subs(units.joule, 1).evalf(7),
             units.joule,
          convert_to(m, units.kilogram).subs(units.kilogram, 1).evalf(4),
             units.kilogram,
          convert_to(g, units.meter / units.second**2).subs({units.meter: 1, units.second: 1}).evalf(3),
             units.meter / units.second**2,
          convert_to(h, units.meter).subs(units.meter,1).evalf(4),
             units.meter))