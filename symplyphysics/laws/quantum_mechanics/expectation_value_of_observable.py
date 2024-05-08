from sympy import (
    Eq,
    Integral,
    Function as SymFunction,
    symbols,
    sqrt,
    conjugate,
    S,
    Expr,
)
from symplyphysics import (
    units,
    Symbol,
    Function,
)
from symplyphysics.core.type_aliases import WaveFunction, Observable

# Description
## Roughly speaking, an observable is a measurable property of a physical system, such as spin, position,
## energy, momentum, etc. In terms of Quantum Mechanics, every physical observable corresponds to an operator
## which acts on the wave function.

# Law: <O> = Integral(conj(psi(x, t)) * O[psi(x, t)], (x, -oo, oo))
## O - observable operator
## psi - wave function
## x - position
## t - time
## oo - infinity
## <O> - expectation value of O
## conj - complex conjugate
## O[f] - operator `O` applied to function `f`

mean_observable_value = symbols("mean_observable_value")
observable = symbols("observable", cls=SymFunction)
wave_function = Function("wave_function", 1 / sqrt(units.length))
position = Symbol("position", units.length, real=True)
time = Symbol("time", units.time)

law = Eq(
    mean_observable_value,
    Integral(
        conjugate(wave_function(position, time))
        * observable(wave_function(position, time)),
        (position, S.NegativeInfinity, S.Infinity)
    ),
)

def calculate_mean_observable_value(
    observable_: Observable,
    wave_function_: WaveFunction,
) -> Expr:
    result = law.rhs.replace(
        wave_function, wave_function_,
    ).replace(
        observable, lambda _wave_function: observable_(wave_function_)(position, time),
    ).doit()
    return result
