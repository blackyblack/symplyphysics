#!/usr/bin/env python3

from sympy import Symbol, pi
from symplyphysics import Quantity, convert_to, print_expression, units
from symplyphysics.laws.dynamics import (kinetic_energy_from_rotational_inertia_and_angular_speed as
    rotational_energy_law)
from symplyphysics.laws.kinematic.rotational_inertia.geometries import (
    solid_disk_about_central_axis as disk_formula)

# Description
## An spin testing was conducted involving a solid steel rotor of mass M = 272 kg and
## a radius R = 38.0 cm. When the sample reached an angular speed of 14000 revolutions
## per minute, the engineers heard an loud sound from the test room: the rotor had, in fact,
## exploded and badly damaged the room. How much energy was released in the explosion
## of the rotor?

mass = Symbol("mass")
radius = Symbol("radius")
angular_speed = Symbol("angular_speed")

revolution = 2 * pi
values = {
    mass: Quantity(272.0 * units.kilogram),
    radius: Quantity(38.0 * units.centimeter),
    angular_speed: Quantity(14000.0 * revolution * units.radian / units.minute)
}

rotational_inertia = disk_formula.law.rhs.subs({
    disk_formula.disk_mass: mass,
    disk_formula.radius: radius,
})

rotational_energy = rotational_energy_law.law.rhs.subs({
    rotational_energy_law.rotational_inertia: rotational_inertia,
    rotational_energy_law.angular_speed: angular_speed,
})

rotational_energy_value = convert_to(Quantity(rotational_energy.subs(values)), units.joule).evalf(3)

print(f"Formula for the rotational energy of the rotor:\n{print_expression(rotational_energy)}\n")
print(f"The energy release due to the rotor's rotation is {rotational_energy_value} J.")
