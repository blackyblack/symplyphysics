from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol,
                           print_expression, validate_input_symbols,
                           validate_output_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import power_from_energy_time as power_and_time
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohm_law
from symplyphysics.laws.electricity import amount_energy_from_voltage_time_resistance as joule_lenz_law

# Description
# Dissipated heat power is proportional to current square and resistance
# P = I**2 * R
# where :
# P - is dissipated heat power of the element
# I is current flowing through this element
# R is impedance of this element

heat_power = Symbol("heat_power", units.power)
current = Symbol("current", units.current)
resistance = Symbol("resistance", units.impedance)

law = Eq(heat_power, current**2 * resistance)

# This law might be easily derived via Joule-Lenz law and dependence of power from energy and time

ohm_law_applied = ohm_law.law.subs({
    ohm_law.voltage: joule_lenz_law.voltage,
    ohm_law.current: current,
    ohm_law.resistance: resistance
})
power_and_time_applied = power_and_time.law.subs({
    power_and_time.energy: joule_lenz_law.amount_energy,
    power_and_time.time: joule_lenz_law.time
})
joule_lenz_law_applied = joule_lenz_law.law.subs(
    {joule_lenz_law.resistance: resistance})

law_derived = [ohm_law_applied, power_and_time_applied, joule_lenz_law_applied]
power_derived = solve(law_derived,
                      (joule_lenz_law.voltage, joule_lenz_law.amount_energy,
                       power_and_time.power),
                      dict=True)[0][power_and_time.power]

# Check if derived power is same as declared
assert (expr_equals(power_derived, law.rhs))


def print() -> str:
    return print_expression(law)


@validate_input_symbols(current_=current, resistance_=resistance)
@validate_output_symbol(heat_power)
def calculate_heat_power(current_: Quantity, resistance_: Quantity) -> Quantity:
    result_power_expr = solve(law, heat_power, dict=True)[0][heat_power]
    result_expr = result_power_expr.subs({
        current: current_,
        resistance: resistance_
    })
    return expr_to_quantity(result_expr)
