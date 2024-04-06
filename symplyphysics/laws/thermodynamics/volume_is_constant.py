from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation as thermodynamics_law

# Description
## Isochoric process: V = const, P1 * T2 = P2 * T1
## Where:
## P is pressure,
## V is volume,
## T is temperature
## T1 is initial temperature, T2 is resulting temperature
## P1 is initial pressure, P1 is resulting pressure

pressure_start = Symbol("pressure_start", units.pressure)
pressure_end = Symbol("pressure_end", units.pressure)
temperature_start = Symbol("temperature_start", units.temperature)
temperature_end = Symbol("temperature_end", units.temperature)

law = Eq(pressure_start * temperature_end, pressure_end * temperature_start)

## Derive the same law from the general ideal gas law

volume_start = Symbol("volume_start", units.volume)
volume_end = Symbol("volume_end", units.volume)

isochoric_condition = Eq(volume_start, volume_end)

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

derived_law = [eq_start, eq_end, isochoric_condition]

## Check the equivalence of 'law' and 'derived_law'
derived_pressure_end = solve(derived_law, (volume_start, volume_end, pressure_end),
    dict=True)[0][pressure_end]
assert solve(law, pressure_end, dict=True)[0][pressure_end] == derived_pressure_end


def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_start_=temperature_start,
    pressure_start_=pressure_start,
    temperature_end_=temperature_end)
@validate_output(pressure_end)
def calculate_pressure(temperature_start_: Quantity, pressure_start_: Quantity,
    temperature_end_: Quantity) -> Quantity:
    solved = solve(law, pressure_end, dict=True)[0][pressure_end]
    result_expr = solved.subs({
        pressure_start: pressure_start_,
        temperature_start: temperature_start_,
        temperature_end: temperature_end_
    })
    return Quantity(result_expr)
