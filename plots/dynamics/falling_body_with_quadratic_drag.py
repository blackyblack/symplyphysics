#!/usr/bin/env python3
"""
A stationary body begins falling in a medium where the drag force exerted on it is proportional to
its speed squared. Plot the body's speed as a function of time.
"""

from sympy import solve, Symbol, Idx, Eq, dsolve, S, sqrt
from sympy.plotting import plot
from symplyphysics import quantities, global_index, print_expression
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.definitions import (
    acceleration_is_speed_derivative as acceleration_def,
    net_force_is_sum_of_individual_forces as force_superposition_law,
)
from symplyphysics.laws.dynamics import (
    acceleration_is_force_over_mass as newtons_second_law,)

mass = newtons_second_law.mass
time = acceleration_def.time
speed = acceleration_def.speed

# projection on -z:
gravity_force_expr = solve(
    newtons_second_law.law,
    newtons_second_law.force,
)[0].subs(
    newtons_second_law.acceleration,
    quantities.acceleration_due_to_gravity,
)

drag_constant = Symbol("b", positive=True)
# TODO: add law for the force of drag, quadratic w.r.t velocity
drag_force_expr = -1 * drag_constant * speed(time)**2  # projection on -z

total_force = force_superposition_law.definition.rhs.subs({
    global_index: Idx("i", (1, 2)),
}).doit().subs({
    force_superposition_law.force[1]: gravity_force_expr,
    force_superposition_law.force[2]: drag_force_expr,
})

newtons_second_eqn = newtons_second_law.law.subs({
    newtons_second_law.acceleration: acceleration_def.acceleration(time),
    newtons_second_law.force: total_force,
})

tau = Symbol("tau", positive=True)
tau_eqn = Eq(tau, sqrt(mass / (quantities.acceleration_due_to_gravity * drag_constant)))

print(
    "Parameters:",
    f"1. `{tau}`\tcharacteristic time",
    f"2. `{mass}`\tmass",
    f"3. `{quantities.acceleration_due_to_gravity}`\tacceleration due to gravity",
    f"4. `{drag_constant}`\tdrag constant",
    sep="\n",
    end="\n\n",
)

print(
    "Relationship between parameters:",
    print_expression(tau_eqn),
    sep="\n",
    end="\n\n",
)

speed_expr = solve(
    (tau_eqn, acceleration_def.definition, newtons_second_eqn),
    (drag_constant, speed(time), acceleration_def.acceleration(time)),
    dict=True,
)[0][speed(time)]
deqn = Eq(speed(time), speed_expr)

speed_expr = dsolve(deqn, speed(time)).rhs

(c1,) = speed_expr.free_symbols - {quantities.acceleration_due_to_gravity, time, tau}
# Assume the body's initial speed is zero
c1_expr = solve(Eq(speed_expr.subs(time, 0), 0), c1)[0]

speed_expr = speed_expr.subs(c1, c1_expr).simplify()

print(
    f"Speed of a falling body as a function of time `{time}`:",
    print_expression(speed_expr),
    sep="\n",
    end="\n\n",
)

asymptote_expr = speed_expr.limit(time, S.Infinity)

print(
    "Maximum attainable speed (at infinite time):",
    print_expression(asymptote_expr),
    sep="\n",
    end="\n\n",
)

base_plot = plot(
    title="Speed of a body falling in a medium with quadratic resistance",
    legend=True,
    show=False,
    xlabel="time, $s$",
    ylabel="speed, $m/s$",
)

taus = 0.1, 0.5, 1.0, 1.5
colors = "blue", "orange", "green", "red"
tmin, tmax = 0, 4

has_asymptote_label = False  # pylint: disable=invalid-name

for tau_, color in zip(taus, colors, strict=True):
    speed_expr_ = evaluate_expression(speed_expr.subs(tau, tau_))
    asymptote_expr_ = speed_expr_.limit(time, S.Infinity)

    asymptote_label = ""  # pylint: disable=invalid-name
    if not has_asymptote_label:
        has_asymptote_label = True  # pylint: disable=invalid-name
        asymptote_label = "asymptote"  # pylint: disable=invalid-name

    subplot = plot(
        asymptote_expr_,
        (time, tmin, tmax),
        label=asymptote_label,
        line_color="gray",
        show=False,
    )
    base_plot.extend(subplot)

    subplot = plot(
        speed_expr_,
        (time, 0, 4),
        label=rf"$\sqrt{{\frac{{m}}{{g b}}}} = {tau_}$",
        line_color=color,
        show=False,
    )
    base_plot.extend(subplot)

base_plot.show()
