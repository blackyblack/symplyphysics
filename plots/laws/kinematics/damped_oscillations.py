#!/usr/bin/env python3

from collections import namedtuple
from sympy import Eq, solve, dsolve, Expr
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.definitions import damped_harmonic_oscillator_equation as damped_eqn

Data = namedtuple("Data", "zeta label color")

DATAS = (
    Data(zeta=0.5, label="underdamping", color="red"),
    Data(zeta=1.0, label="critical damping", color="blue"),
    Data(zeta=2.0, label="overdamping", color="magenta"),
)

OMEGA = 1.0
INITIAL_POSITION = 1.0
INITIAL_VELOCITY = 0.0

def get_solution(data: Data) -> Expr:
    # eqn = damped_eqn.definition.subs({
    #     damped_eqn.damping_ratio: data.zeta,
    #     damped_eqn.undamped_angular_frequency: OMEGA,
    # })
    # dsolved = dsolve(eqn, damped_eqn.displacement(damped_eqn.time)).rhs
    coefs = solve(
        [
            Eq(INITIAL_POSITION, dsolved.subs(damped_eqn.time, 0)),
            Eq(INITIAL_VELOCITY, dsolved.diff(damped_eqn.time).subs(damped_eqn.time, 0)),
        ],
        ("C1", "C2"),
        dict=True,
    )[0]
    dsolved_subs = dsolved.subs(coefs)
    print(dsolved_subs.free_symbols)
    return dsolved_subs


p = plot(
    title="Damped oscillations for various values of damping ratio",
    xlabel="time, s",
    ylabel="position, m",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
    annotations=None,
)

for data in DATAS:
    sub_p = plot(
        get_solution(data),
        (damped_eqn.time, 0, 5),
        label=data.label,
        line_color=data.color,
        show=False
    )
    p.append(sub_p[0])

p.show()
