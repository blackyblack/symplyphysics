#!/usr/bin/env python3

from sympy import symbols, Eq, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.laws.waves import frequency_shift_from_velocity as shift_law

reduced_frequency = symbols("reduced_frequency")

reduced_frequency_eqn = Eq(
    reduced_frequency,
    shift_law.observed_frequency / shift_law.real_frequency,
)

reduced_observer_speed = symbols("reduced_observer_speed")

reduced_observer_speed_eqn = Eq(
    reduced_observer_speed,
    shift_law.observer_velocity / shift_law.wave_velocity,
)

reduced_source_speed = symbols("reduced_source_speed")

reduced_source_speed_eqn = Eq(
    reduced_source_speed,
    shift_law.source_velocity / shift_law.wave_velocity,
)

reduced_frequency_expr = solve(
    (
    shift_law.law,
    reduced_frequency_eqn,
    reduced_observer_speed_eqn,
    reduced_source_speed_eqn,
    ),
    (
    shift_law.real_frequency,
    reduced_frequency,
    shift_law.observer_velocity,
    shift_law.source_velocity,
    ),
    dict=True,
)[0][reduced_frequency]

print("Definition of reduced frequency:\n")
print(print_expression(reduced_frequency_eqn))

print("\nDefinition of reduced observer speed:\n")
print(print_expression(reduced_observer_speed_eqn))

print("\nDefinition of reduced source speed:\n")
print(print_expression(reduced_source_speed_eqn))

print("\nEquation for reduced frequency:\n")
print(print_expression(Eq(reduced_frequency, reduced_frequency_expr)))

# `v` is wave speed, `v_S` is source speed, and `v_O` is observer speed
# `f_O` is observed frequency, and `f_S` is source frequency

base_plot = plot(
    title="Reduced observed frequency as a function of reduced observer speed",
    xlabel=r"reduced observer speed $\frac{v_O}{v}$",
    ylabel=r"reduced observed frequency $\frac{f_O}{f_S}$",
    legend=True,
    backend=MatplotlibBackend,
    show=False,
)

for reduced_source_speed_ in (-0.8, -0.5, 0, 0.5, 1):
    expr = reduced_frequency_expr.subs(reduced_source_speed, reduced_source_speed_)
    sub_plot = plot(
        expr,
        (reduced_observer_speed, -2, 1),
        label=rf"$\frac{{v_S}}{{v}} = {reduced_source_speed_}$",
        show=False,
    )
    base_plot.extend(sub_plot)

base_plot.show()
