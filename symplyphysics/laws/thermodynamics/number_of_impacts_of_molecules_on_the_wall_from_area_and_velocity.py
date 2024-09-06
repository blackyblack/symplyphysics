from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to_float)

# Description
## Consider the volume of gas in a vessel in the form of a cube and determine the number of impacts on its wall perpendicular to the x-axis.
## Denote by n the number of molecules per unit volume whose velocity projections on the x-axis are ±v_x.
## In chaotic motion, the number of molecules with positive and negative projections on a given axis can be considered the same and equal to n / 2, since movement in all directions is equally likely.
## Of the selected group of molecules, only those molecules whose velocities are directed towards the wall and which are located from the wall at a distance not
## exceeding v_x * t, or those molecules that are inside the volume V = S * v_x * t, will reach the wall area S over a period of time t.

## Law: N = 0.5 * n * S * v_x * t
## Where:
## N is number of impacts
## n is the number of molecules per unit volume
## S is wall area
## v_x is velocity projection on the x-axis
## t is time

number_of_impacts = Symbol("number_of_impacts", dimensionless)
molecules_concentration = Symbol("molecules_concentration", 1 / units.volume)
area = Symbol("area", units.area)
velocity_projection = Symbol("velocity_projection", units.velocity)
time = Symbol("time", units.time)

law = Eq(number_of_impacts, (molecules_concentration * area * velocity_projection * time) / 2)

## Conditions
## Gas is ideal
## Wall is flat and perpendicular to X-axis


def print_law() -> str:
    return print_expression(law)


@validate_input(molecules_concentration_=molecules_concentration,
    area_=area,
    velocity_projection_=velocity_projection,
    time_=time)
@validate_output(number_of_impacts)
def calculate_number_of_impacts(molecules_concentration_: Quantity, area_: Quantity,
    velocity_projection_: Quantity, time_: Quantity) -> int:
    result_expr = solve(law, number_of_impacts, dict=True)[0][number_of_impacts]
    result_number_of_impacts = result_expr.subs({
        molecules_concentration: molecules_concentration_,
        area: area_,
        velocity_projection: velocity_projection_,
        time: time_
    })
    return int(convert_to_float(result_number_of_impacts))
