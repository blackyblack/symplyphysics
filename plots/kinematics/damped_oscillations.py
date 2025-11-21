#!/usr/bin/env python3
r"""
Plot the displacement during damped oscillations for different values of the damping ratio
:math:`\zeta`.
"""

from collections import namedtuple
from sympy import dsolve, Expr, symbols, Function as SymFunction
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import units, convert_to_si
from symplyphysics.definitions import damped_harmonic_oscillator_equation as damped_eqn

displacement = symbols("q", cls=SymFunction, real=True)
time = symbols("t", positive=True)

Datum = namedtuple("Datum", "zeta label color")

DATA = (
    Datum(zeta=0.3, label="underdamping", color="green"),
    Datum(zeta=0.999999, label="critical damping",
    color="red"),  # cannot put 1.0 due to rounding errors
    Datum(zeta=2, label="overdamping", color="blue"),
)

OMEGA = convert_to_si(1.5 * units.radian / units.second)
INITIAL_POSITION = convert_to_si(1.0 * units.meter)
INITIAL_VELOCITY = convert_to_si(-3.0 * units.meter / units.second)

initial_conditions = {
    displacement(0): INITIAL_POSITION,
    displacement(time).diff(time).subs(time, 0): INITIAL_VELOCITY,
}

eqn = damped_eqn.definition.subs(damped_eqn.time, time).subs({
    damped_eqn.displacement(time): displacement(time),
    damped_eqn.undamped_angular_frequency: OMEGA,
})


def get_solution(zeta: float) -> Expr:
    eqn_subs = eqn.subs(damped_eqn.damping_ratio, zeta)
    dsolved = dsolve(
        eqn_subs,
        displacement(time),
        ics=initial_conditions,
    ).rhs
    return dsolved


p = plot(
    title="Damped oscillations for various values of damping ratio",
    xlabel="time, s",
    ylabel="displacement, m",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
)

for datum in DATA:
    sol = get_solution(datum.zeta)
    sub_p = plot(sol, (time, 0, 10), label=datum.label, line_color=datum.color, show=False)
    p.append(sub_p[0])

p.show()
