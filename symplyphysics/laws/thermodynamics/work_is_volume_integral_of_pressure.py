from sympy import Eq, Integral, Point2D, solve, symbols, Function as SymFunction
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.geometry.line import two_point_function
from symplyphysics.laws.thermodynamics import infinitesimal_work_in_quasistatic_process as work_law

# Description
## The pressure-volume work occurs when the volume of a system changes.

# Law: W = Integral(p(V), V)
## W - work done by gas
## p - gas pressure
## V - gas volume

# Conditions
## - The process is reversible or quasi-static.
## - The system is closed.

work = Symbol("work", units.energy)
pressure = Function("pressure", units.pressure)
volume = Symbol("volume", units.volume)
volume_before = Symbol("volume_before", units.volume)
volume_after = Symbol("volume_after", units.volume)

law = Eq(work, Integral(pressure(volume), (volume, volume_before, volume_after)))

# Derive law from infinitesimal counterpart

_work_function = symbols("work", cls=SymFunction)

# Prepare law for integration
# Note that `dW = (dW/dV) * dV` and we set `dV` to 1
_infinitesimal_law = work_law.law.subs({
    work_law.infinitesimal_work_done: _work_function(volume).diff(volume),
    work_law.infinitesimal_volume_change: 1,
    work_law.pressure_inside_system: pressure(volume),
})

_integrated_law = Eq(
    _infinitesimal_law.lhs.integrate((volume, volume_before, volume_after)),
    _infinitesimal_law.rhs.integrate((volume, volume_before, volume_after))
).subs({
    _work_function(volume_before): 0,  # there was no work done at the start of the process
    _work_function(volume_after): work,  # the work at the end of the process is the work done during the process
})

_work_derived = solve(_integrated_law, work)[0]

assert expr_equals(_work_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    volume_before_=volume_before,
    volume_after_=volume_after,
    pressure_before_=pressure,
    pressure_after_=pressure,
)
@validate_output(work)
def calculate_work(
    volume_before_: Quantity,
    volume_after_: Quantity,
    pressure_before_: Quantity,
    pressure_after_: Quantity,
) -> Quantity:
    pressure_ = two_point_function(
        Point2D(volume_before_, pressure_before_),
        Point2D(volume_after_, pressure_after_),
        volume,
    )
    result = law.rhs.subs({
        pressure(volume): pressure_,
        volume_before: volume_before_,
        volume_after: volume_after_,
    }).doit()
    return Quantity(result)
