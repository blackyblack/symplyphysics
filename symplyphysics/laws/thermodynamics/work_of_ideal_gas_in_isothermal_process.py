from sympy import Eq, solve, log, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import work_is_volume_integral_of_pressure as work_law
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation
from symplyphysics.laws.quantities import (
    quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law,
)

## Description
## The isothermal process of expansion (or compression) of a gas can occur under conditions where heat exchange between the gas and the external environment is carried out at a constant temperature difference.
## To do this, the heat capacity of the external environment must be large enough, and the expansion (or compression) process must be slow enough.

## Law: A = (m / mu) * R * T * ln(V_2 / V_1)
## Where:
## A is work of gas
## m is mass
## mu is molar mass
## V_1 is start volume
## V_2 is final volume
## R is ideal gas constant,
## T is temperature

## Note
## Since the internal energy of an ideal gas does not change during the isothermal process, all the heat received from the source is converted into work.
## When gas is compressed, the work of external forces is positive, and when expanding, the gas does the positive work.

work = Symbol("work", units.energy)
molar_mass = Symbol("molar_mass", units.mass / units.amount_of_substance)
start_volume = Symbol("start_volume", units.volume)
final_volume = Symbol("final_volume", units.volume)
temperature = symbols.thermodynamics.temperature
gas_mass = clone_symbol(symbols.basic.mass, "gas_mass")

law = Eq(work, (gas_mass / molar_mass) * units.molar_gas_constant * temperature *
    log(final_volume / start_volume))

# Derive from ideal gas equation

_mole_count = solve(
    molar_qty_law.law, molar_qty_law.amount_of_substance,
)[0].subs({
    molar_qty_law.extensive_quantity: gas_mass,
    molar_qty_law.molar_quantity: molar_mass,
})

_ideal_gas_eqn = ideal_gas_equation.law.subs({
    ideal_gas_equation.mole_count: _mole_count,
    ideal_gas_equation.temperature: temperature,
})

_pressure_expr = solve(
    _ideal_gas_eqn, ideal_gas_equation.pressure
)[0].subs(
    ideal_gas_equation.volume, work_law.volume
)

_volume_before = SymSymbol("volume_before", positive=True)
_volume_after = SymSymbol("volume_after", positive=True)

_work_expr = work_law.law.rhs.subs({
    work_law.volume_before: _volume_before,
    work_law.volume_after: _volume_after,
    work_law.pressure(work_law.volume): _pressure_expr
}).doit().simplify()

_work_from_law = law.rhs.subs({
    start_volume: _volume_before,
    final_volume: _volume_after,
})

assert expr_equals(_work_expr, _work_from_law)


@validate_input(mass_=gas_mass,
    molar_mass_=molar_mass,
    start_volume_=start_volume,
    final_volume_=final_volume,
    temperature_=temperature)
@validate_output(work)
def calculate_work(mass_: Quantity, molar_mass_: Quantity, start_volume_: Quantity,
    final_volume_: Quantity, temperature_: Quantity) -> Quantity:
    solved = solve(law, work, dict=True)[0][work]
    result_expr = solved.subs({
        gas_mass: mass_,
        molar_mass: molar_mass_,
        start_volume: start_volume_,
        final_volume: final_volume_,
        temperature: temperature_
    })
    return Quantity(result_expr)
