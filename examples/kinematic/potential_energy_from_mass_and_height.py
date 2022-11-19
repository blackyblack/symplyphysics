from symplyphysics import (
    units, convert_to,solve, SI
)
from symplyphysics.laws.kinematic import potential_energy_from_mass_and_height as potential_energy

#A bullet weighing 9 grams is fired from the gun vertically upwards and reaches
#an altitude of 1,500 meters. Under the action of gravity it falls freely
#to the ground and hits a shard of glass lying on the ground
#Can a bullet break a piece of glass on the ground if it takes 100 joules of work ?

m = units.Quantity('m')
SI.set_quantity_dimension(m, units.mass)
SI.set_quantity_scale_factor(m, 0.009 * units.kilogram)
h = units.Quantity('h')
SI.set_quantity_dimension(h, units.length)
SI.set_quantity_scale_factor(h, 1500 * units.meter)
print("\nFormula is:\n{}".format(potential_energy.print()))

result = potential_energy.calculate_potential_energy(m, h)
print("\nPotential energy = {} {}; for mass = {} {} and height = {} {}\n"
      .format(
          convert_to(result, units.joule).subs(units.joule, 1).evalf(7),
             units.joule,
          convert_to(m, units.kilogram).subs(units.kilogram, 1).evalf(4),
             units.kilogram,
          convert_to(h, units.meter).subs(units.meter,1).evalf(4),
             units.meter))