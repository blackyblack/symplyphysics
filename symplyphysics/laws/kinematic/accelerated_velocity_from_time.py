from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Accelerated velocity is time dependent and increases with time if acceleration is co-directed with velocity and decreases if they are counter-directed.
## Initial velocity is velocity at the start of observation, when time is 0.

## Law: V = V0 + at
## Where:
## V is the velocity in the moment of time t
## V0 is the initial velocity (in the moment of time 0)
## a is the acceleration.

velocity = Symbol("velocity", units.velocity)
time = Symbol("time", units.time)
acceleration = Symbol("acceleration", units.acceleration)
initial_velocity = Symbol("initial_velocity", units.velocity)

law = Eq(velocity, initial_velocity + acceleration * time)

def print(expr: Expr) -> str:
    symbols = [velocity, time, acceleration, initial_velocity]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(initial_velocity_=initial_velocity, acceleration_=acceleration, time_=time)
@validate_output_symbol(velocity)
def calculate_velocity(initial_velocity_: Quantity, acceleration_: Quantity, time_: Quantity) -> Quantity:
    result_velocity_expression = solve(law, velocity, dict=True)[0][velocity]
    result_expr = result_velocity_expression.subs({initial_velocity: initial_velocity_, acceleration: acceleration_, time: time_})
    return expr_to_quantity(result_expr)
