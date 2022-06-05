#!/usr/bin/env python3
from sympy import solve, pretty, symbols, pi
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.laws.dynamics import buoyant_force_from_density_and_volume as archimedes_law
from symplyphysics.laws.dynamics import acceleration_from_force as gravity_law

print("Formula is:\n{}".format(archimedes_law.print()))

height = symbols('height')

fluid_density = 0.6
cylinder_mass = 30
cylinder_height = 5
cylinder_radius = 2

cylinder_volume_function = pi * height * cylinder_radius**2
cylinder_volume = cylinder_volume_function.subs(height, cylinder_height)

solved = abs(solve(archimedes_law.law, archimedes_law.force_buoyant, dict=True)[0][archimedes_law.force_buoyant])
result_buoyant_force_above_liquid = solved.subs({
    archimedes_law.units.acceleration_due_to_gravity: 9.8,
    archimedes_law.fluid_density: fluid_density,
    archimedes_law.displaced_volume: cylinder_volume_function})

result_buoyant_force_below_liquid = solved.subs({
    archimedes_law.units.acceleration_due_to_gravity: 9.8,
    archimedes_law.fluid_density: fluid_density,
    archimedes_law.displaced_volume: cylinder_volume})

solved_gravity = solve(gravity_law.law, gravity_law.force, dict=True)[0][gravity_law.force]
result_gravity_force = solved_gravity.subs({
    gravity_law.mass: cylinder_mass,
    gravity_law.acceleration: 9.8})

print("Buoyant force above liquid function is:\n{}".format(
    pretty(result_buoyant_force_above_liquid, use_unicode=False)))
print("Buoyant force below liquid function is:\n{}".format(
    pretty(result_buoyant_force_below_liquid, use_unicode=False)))

p1 = plot(
    result_buoyant_force_above_liquid,
    (height, 0, cylinder_height),
    line_color='blue',
    title='Floating body',
    xlabel='Height below water',
    ylabel='Force',
    label='Buoyant',
    backend=MatplotlibBackend,
    legend=True,
    show=False)

p2 = plot(
    result_gravity_force,
    (height, 0, 8),
    line_color='red',
    label='Gravity',
    backend=MatplotlibBackend,
    show=False)

p3 = plot(
    result_buoyant_force_below_liquid,
    (height, cylinder_height, 8),
    line_color='blue',
    label='',
    backend=MatplotlibBackend,
    show=False)

p1.append(p2[0])
p1.append(p3[0])
p1[2].label = ''
p1.show()