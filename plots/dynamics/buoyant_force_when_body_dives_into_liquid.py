#!/usr/bin/env python3
from sympy import solve, symbols, pi
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, quantities
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.laws.dynamics import buoyant_force_from_density_and_volume as archimedes_law
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as gravity_law

print(f"Formula is:\n{print_expression(archimedes_law.law)}")

height = symbols("height")

FLUID_DENSITY = 0.6
CYLINDER_MASS = 30
CYLINDER_HEIGHT = 5
CYLINDER_RADIUS = 2

cylinder_volume_function = pi * height * CYLINDER_RADIUS**2
cylinder_volume = cylinder_volume_function.subs(height, CYLINDER_HEIGHT)

solved = abs(
    solve(archimedes_law.law, archimedes_law.buoyant_force,
    dict=True)[0][archimedes_law.buoyant_force])
result_buoyant_force_above_liquid = solved.subs({
    archimedes_law.fluid_density: FLUID_DENSITY,
    archimedes_law.displaced_volume: cylinder_volume_function
})
result_buoyant_force_above_liquid = evaluate_expression(result_buoyant_force_above_liquid)

result_buoyant_force_below_liquid = solved.subs({
    archimedes_law.fluid_density: FLUID_DENSITY,
    archimedes_law.displaced_volume: cylinder_volume
})
result_buoyant_force_below_liquid = evaluate_expression(result_buoyant_force_below_liquid)

solved_gravity = solve(gravity_law.law, gravity_law.force, dict=True)[0][gravity_law.force]
result_gravity_force = solved_gravity.subs({
    gravity_law.mass: CYLINDER_MASS,
    gravity_law.acceleration: quantities.acceleration_due_to_gravity,
})
result_gravity_force = evaluate_expression(result_gravity_force)

print(
    f"Buoyant force above liquid function is:\n{print_expression(result_buoyant_force_above_liquid)}"
)
print(
    f"Buoyant force below liquid function is:\n{print_expression(result_buoyant_force_below_liquid)}"
)

p1 = plot(result_buoyant_force_above_liquid, (height, 0, CYLINDER_HEIGHT),
    line_color="blue",
    title="Floating body",
    xlabel="Height below water",
    ylabel="Force",
    label="Buoyant",
    backend=MatplotlibBackend,
    legend=True,
    show=False)

p2 = plot(result_gravity_force, (height, 0, 8),
    line_color="red",
    label="Gravity",
    backend=MatplotlibBackend,
    show=False)

p3 = plot(result_buoyant_force_below_liquid, (height, CYLINDER_HEIGHT, 8),
    line_color="blue",
    label="",
    backend=MatplotlibBackend,
    show=False)

p1.append(p2[0])
p1.append(p3[0])
p1[2].label = ""
p1.show()
