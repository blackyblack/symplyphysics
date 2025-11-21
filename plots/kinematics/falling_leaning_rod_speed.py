#!/usr/bin/env python3
"""
A rigid rod :math:`AB` of length :math:`l` leans with both ends onto the floor (with end
:math:`A`) and the wall (with end :math:`B`). The initial coordinate of end :math:`A` is
:math:`x_0`. Find the expression for the position of end :math:`B` as a function of time :math:`t`
as end :math:`A` moves along the floor with constant speed :math:`v` in the positive direction of
the :math:`x`-axis.
"""

from sympy import Eq, solve, ask, Q
from sympy.plotting import plot
from symplyphysics import symbols, clone_as_symbol, print_expression
from symplyphysics.laws.kinematics import position_via_constant_speed_and_time as position_law

a_initial_position = clone_as_symbol(symbols.position, subscript="0")
a_speed = clone_as_symbol(symbols.speed, positive=True)
time = clone_as_symbol(position_law.time, positive=True)

a_position_expr = position_law.law.rhs.subs({
    position_law.initial_position: a_initial_position,
    position_law.speed: a_speed,
    position_law.time: time,
})

b_position = clone_as_symbol(symbols.position, display_symbol="y", positive=True)

# rod length
length = clone_as_symbol(symbols.length, positive=True)

# From the Pythagorean theorem:
length_eqn = Eq(length**2, a_position_expr**2 + b_position**2)

(b_position_expr,) = filter(
    lambda e: not e.could_extract_minus_sign(),
    solve(length_eqn, b_position),
)

(t_max_expr,) = filter(
    lambda e: ask(Q.positive(e), Q.positive(length - a_initial_position)),
    solve(Eq(b_position_expr, 0), time),
)

b_speed_expr = b_position_expr.diff(time).simplify()

print(
    "Position of end B:",
    print_expression(b_position_expr),
    sep="\n",
    end="\n\n",
)
print(
    "Speed of end B:",
    print_expression(b_speed_expr),
    sep="\n",
    end="\n\n",
)
print(
    "Total time of fall:",
    print_expression(t_max_expr),
    sep="\n",
)

subs = {
    length: 1,
    a_initial_position: 0.5,
    a_speed: 0.5,
}

t_max_ = t_max_expr.subs(subs)
b_position_ = b_position_expr.subs(subs)
b_speed_ = b_speed_expr.subs(subs)

plot(
    b_position_,
    (time, 0, t_max_),
    title="Position of end B as a function of time",
    xlabel="time, s",
    ylabel="position, m",
)

# NOTE that the speed of end B explodes the closer it gets to the ground. In real life there would
# be friction that prevents this behaviour from occurring. This corresponds to end A abruptly
# stopping (and its speed changing to 0) when the rod fully lies on the floor.
plot(
    b_speed_,
    (time, 0, t_max_),
    title="Speed of end B as a function of time",
    xlabel="time, s",
    ylabel="speed, m/s",
)
