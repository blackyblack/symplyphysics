from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    Dimensionless, validate_input, validate_output)
from symplyphysics.laws.thermodynamics import pressure_from_temperature_and_volume as thermodynamics_law

# Description
## Adiabatic process: Q = 0, P * V^y = const
## Where:
## Q is amount of transferred heat between system and environment
## P is pressure,
## V is volume,
## y is the ratio of specific heats (also known as heat capacity
##   ratio) (https://en.wikipedia.org/wiki/Heat_capacity_ratio)

specific_heats_ratio = Symbol("specific_heats_ratio", Dimensionless)
# Some of these parameters depend on each other. It is up to user, which of these parameters to choose
# as known.
temperature_start = Symbol("temperature_start", units.temperature)
temperature_end = Symbol("temperature_end", units.temperature)
volume_start = Symbol("volume_start", units.volume)
volume_end = Symbol("volume_end", units.volume)
pressure_start = Symbol("pressure_start", units.pressure)
pressure_end = Symbol("pressure_end", units.pressure)

adiabatic_condition = Eq(pressure_start * (volume_start**specific_heats_ratio),
    pressure_end * (volume_end**specific_heats_ratio))

eq_start = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: temperature_start,
    thermodynamics_law.volume: volume_start,
    thermodynamics_law.pressure: pressure_start
})

eq_end = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: temperature_end,
    thermodynamics_law.volume: volume_end,
    thermodynamics_law.pressure: pressure_end
})

law = [eq_start, eq_end, adiabatic_condition]


def print() -> str:
    return print_expression(law)


@validate_input(mole_count_=thermodynamics_law.mole_count,
    temperature_start_=temperature_start,
    volume_start_=volume_start,
    volume_end_=volume_end,
    specific_heats_ratio_=specific_heats_ratio)
@validate_output(pressure_end)
def calculate_pressure(mole_count_: Quantity, temperature_start_: Quantity, volume_start_: Quantity,
    volume_end_: Quantity, specific_heats_ratio_: float) -> Quantity:

    solved = solve(law, (pressure_start, temperature_end, pressure_end), dict=True)[0][pressure_end]
    result_pressure = solved.subs({
        thermodynamics_law.mole_count: mole_count_,
        temperature_start: temperature_start_,
        volume_start: volume_start_,
        volume_end: volume_end_,
        specific_heats_ratio: specific_heats_ratio_
    })
    return expr_to_quantity(result_pressure)
