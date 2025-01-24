#!/usr/bin/env python3

from sympy import Symbol, Eq, solve, oo, plot
from symplyphysics import print_expression
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.quantities import speed_of_light
from symplyphysics.laws.conservation import (
    mechanical_energy_after_equals_to_mechanical_energy_before as energy_conservation_law,
    momentum_after_collision_equals_to_momentum_before as momentum_conservation_law,
)
from symplyphysics.laws.relativistic import (
    total_energy_via_momentum_and_rest_mass as energy_momentum_law,
    total_energy_via_relativistic_mass as energy_is_mass,
    relativistic_mass as moving_mass_law,
    relativistic_momentum as moving_momentum_law,
)

# Description
## A flat light wave hits an ideal flat mirror at rest, perpendicular to its surface. Under the influence of the force
## of light pressure, the mirror starts moving. What is the final speed of the mirror and what is the energy of the
## reflected wave? The rest mass of the mirror `m_0` and the energy of the incident wave `W_0` are known.

incident_wave_energy = Symbol("W_0", positive=True)
reflected_wave_energy = Symbol("W_1", positive=True)
mirror_rest_mass = Symbol("m_0", positive=True)
reduced_mirror_speed = Symbol("u", real=True)

mirror_speed_expr_ = reduced_mirror_speed * speed_of_light

# Conservation of energy

mirror_rest_energy = energy_is_mass.law.rhs.subs(energy_is_mass.relativistic_mass, mirror_rest_mass)

total_energy_before = incident_wave_energy + mirror_rest_energy

mirror_mass_expr = moving_mass_law.law.rhs.subs({
    moving_mass_law.rest_mass: mirror_rest_mass,
    moving_mass_law.speed: mirror_speed_expr_,
})

moving_mirror_energy = energy_is_mass.law.rhs.subs(energy_is_mass.relativistic_mass, mirror_mass_expr)

total_energy_after = reflected_wave_energy + moving_mirror_energy

energy_conservation_eqn = energy_conservation_law.law.subs({
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_before):
        total_energy_before,
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_after):
        total_energy_after,
})

# Conservation of momentum

# Electromagnetic waves have zero rest mass.
electromagnetic_wave_momentum_expr = solve(
    energy_momentum_law.law.subs(energy_momentum_law.invariant_mass, 0),
    energy_momentum_law.relativistic_momentum,
)[1]

incident_wave_momentum = electromagnetic_wave_momentum_expr.subs(
    energy_momentum_law.relativistic_energy,
    incident_wave_energy,
)

total_momentum_before = incident_wave_momentum

reflected_wave_momentum = electromagnetic_wave_momentum_expr.subs(
    energy_momentum_law.relativistic_energy,
    reflected_wave_energy,
)

moving_mirror_momentum = moving_momentum_law.law.rhs.subs({
    moving_momentum_law.rest_mass: mirror_rest_mass,
    moving_momentum_law.speed: mirror_speed_expr_,
})

# The reflected wave moves in the opposite direction
total_momentum_after = moving_mirror_momentum - reflected_wave_momentum

momentum_conservation_eqn = momentum_conservation_law.law.subs({
    momentum_conservation_law.momentum(momentum_conservation_law.time_before):
        total_momentum_before,
    momentum_conservation_law.momentum(momentum_conservation_law.time_after):
        total_momentum_after,
})

solved = solve(
    (energy_conservation_eqn, momentum_conservation_eqn),
    (reduced_mirror_speed, reflected_wave_energy),
    dict=True,
)[0]

reduced_mirror_speed_expr = solved[reduced_mirror_speed]
reflected_wave_energy_expr = solved[reflected_wave_energy]

# Simplify the expressions by introducing a reduced energy symbol

reduced_incident_wave_energy = Symbol("w_0", positive=True)

reduced_incident_wave_energy_eqn = Eq(
    reduced_incident_wave_energy,
    incident_wave_energy / (mirror_rest_mass * speed_of_light**2),
)

reduced_mirror_speed_expr = solve(
    (Eq(reduced_mirror_speed, reduced_mirror_speed_expr), reduced_incident_wave_energy_eqn),
    (reduced_mirror_speed, incident_wave_energy),
    dict=True,
)[0][reduced_mirror_speed]

reflected_wave_energy_expr = solve(
    (Eq(reflected_wave_energy, reflected_wave_energy_expr), reduced_incident_wave_energy_eqn),
    (reflected_wave_energy, incident_wave_energy),
    dict=True,
)[0][reflected_wave_energy].replace(
    4 * reduced_incident_wave_energy**2 + 4 * reduced_incident_wave_energy + 1,
    (2 * reduced_incident_wave_energy + 1)**2,
)

# Print formulas

print("Definition of reduced energy of incident wave:\n")
print(print_expression(reduced_incident_wave_energy_eqn), end="\n\n\n")

print("Formula of reduced speed of mirror:\n")
print(print_expression(reduced_mirror_speed_expr), end="\n\n\n")

reduced_reflected_wave_energy_expr = evaluate_expression(reflected_wave_energy_expr /
    (mirror_rest_mass * speed_of_light**2))

print("Formula of reduced energy of reflected wave:\n")
print(print_expression(reduced_reflected_wave_energy_expr), end="\n\n\n")

reduced_reflected_wave_energy_upper_limit = evaluate_expression(reduced_reflected_wave_energy_expr.limit(
    reduced_incident_wave_energy, oo))

print("Upper limit of reduced energy of reflected wave:")
print(print_expression(reduced_reflected_wave_energy_upper_limit))

# Plot formulas

# Notation
## - `W_0, W_1` - energy of incident and reflected waves respectively
## - `w_0, w_1` - reduced energy of incident and reflected waves respectedly
## - `m_0` - rest mass of mirror
## - `v` - speed of mirror after interaction
## - `beta` - reduced speed of mirror after interaction
## - `c` - speed of light

base_plot = plot(
    title=r"Mirror speed and reflected wave energy as functions of incident wave energy",
    xlabel=r"reduced incident wave energy $w_0 = \frac{W_0}{m_0 c^2}$",
    ylabel="reduced quantities",
    legend=True,
    show=False,
)

plot_range = (reduced_incident_wave_energy, 0, 10)

speed_plot = plot(
    reduced_mirror_speed_expr,
    plot_range,
    label=r"reduced mirror speed $\beta = \frac{v}{c}$",
    line_color="blue",
    show=False,
)
base_plot.extend(speed_plot)

energy_plot = plot(
    reduced_reflected_wave_energy_expr,
    plot_range,
    label=r"reduced reflected wave energy $w_1 = \frac{W_1}{m_0 c^2}$",
    line_color="red",
    show=False,
)
base_plot.extend(energy_plot)

energy_asymptote_plot = plot(
    reduced_reflected_wave_energy_upper_limit,
    plot_range,
    label=r"$w_1$ in the limit $w_0 \to \infty$",
    line_color="pink",
    show=False,
)
base_plot.extend(energy_asymptote_plot)

base_plot.show()
