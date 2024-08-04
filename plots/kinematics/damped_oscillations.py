#!/usr/bin/env python3

from collections import namedtuple
from sympy import dsolve, Expr, symbols, Function as SymFunction
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.definitions import damped_harmonic_oscillator_equation as damped_eqn

displacement = symbols("displacement", cls=SymFunction, real=True)
time = symbols("time", positive=True)

Data = namedtuple("Data", "zeta label color")

DATAS = (
    Data(zeta=0.3, label="underdamping", color="green"),
    Data(zeta=0.999999, label="critical damping",
    color="red"),  # cannot put 1.0 due to rounding errors
    Data(zeta=2, label="overdamping", color="blue"),
)

OMEGA = 1.5
INITIAL_POSITION = 1.0
INITIAL_VELOCITY = -3.0

initial_conditions = {
    displacement(0): INITIAL_POSITION,
    displacement(time).diff(time).subs(time, 0): INITIAL_VELOCITY,
}

eqn = damped_eqn.definition.subs(
    damped_eqn.time,
    time,
).subs({
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

for data in DATAS:
    sol = get_solution(data.zeta)
    sub_p = plot(sol, (time, 0, 10), label=data.label, line_color=data.color, show=False)
    p.append(sub_p[0])

p.show()
