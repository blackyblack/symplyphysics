from symplyphysics import (
    units, convert_to,solve, SI,Quantity, expr_to_quantity
)
from symplyphysics import (
    symbols, Eq, pretty
)
from symplyphysics.laws.electricity import amount_energy_from_voltage_time_resistance as joule_lenz_law
from symplyphysics.definitions import amount_energy_from_mass_and_temperature as operate_energy
from symplyphysics.definitions import density_from_mass_volume as density
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohm_law
from symplyphysics.definitions import power_is_proportional_voltage_and_current as operate_power

# The household electric kettle heated 0.5 liters of water to boiling.
# The power of the kettle is 1500 watts. How long will it take for the
# temperature of the water to drop to a comfortable 35 degrees Celsius?
# Since physics uses Kelvin temperature units, we must consider that 0 degrees Celsius = 273 Kelvin.
# Density of water 1000 kilogram / meter**3
density_example = units.Quantity('density_example')
SI.set_quantity_dimension(density_example, units.mass / units.volume)
SI.set_quantity_scale_factor(density_example, 1000 * units.kilogram / units.meter**3)
P_example = units.Quantity('P_example')
SI.set_quantity_dimension(P_example, units.power)
SI.set_quantity_scale_factor(P_example, 1500 * units.watt)
v_example = units.Quantity('v_example')
SI.set_quantity_dimension(v_example, units.volume)
SI.set_quantity_scale_factor(v_example, 0.5 * units.liter)
C_example = units.Quantity('C_example')
SI.set_quantity_dimension(C_example, units.energy / (units.mass * units.temperature))
SI.set_quantity_scale_factor(C_example, 4200 * units.joule / (units.kilogram * units.kelvin))
t1_example = units.Quantity('t1_example')
SI.set_quantity_dimension(t1_example, units.temperature)
SI.set_quantity_scale_factor(t1_example, (100 + 273) * units.kelvin)
t2_example = units.Quantity('t2_example')
SI.set_quantity_dimension(t2_example, units.temperature)
SI.set_quantity_scale_factor(t2_example, (35 + 273) * units.kelvin)
Celsius = symbols('Celsius')
body_mass = solve(density.definition, density.mass, dict=True)[0][density.mass]
resistance = solve(ohm_law.law, ohm_law.resistance, dict=True)[0][ohm_law.resistance]
current_1 = solve(ohm_law.law, ohm_law.current, dict=True)[0][ohm_law.current]
current_2 = solve(operate_power.law,operate_power.current, dict=True)[0][operate_power.current]
current = Eq(current_1, current_2)
resistance_solved = solve(current, ohm_law.resistance, dict=True)[0][ohm_law.resistance]
law = Eq(operate_energy.law.subs({operate_energy.body_mass: body_mass}),
    joule_lenz_law.law.subs({joule_lenz_law.resistance: resistance_solved}))
solved = solve(law, joule_lenz_law.time, dict=True)[0][joule_lenz_law.time]
time = Eq(joule_lenz_law.time, solved)
print("\nFormula is:\n\n {}".format(pretty(time, use_unicode=False)))
law = law.subs({operate_power.power: P_example, operate_energy.body_mass: body_mass,
    operate_energy.specific_heat_capacity:C_example,
    operate_energy.temperature_begin: t1_example, operate_energy.temperature_end: t2_example,
    density.density: density_example, density.volume: v_example
})
solved_example = solve(law, joule_lenz_law.time, dict=True)[0][joule_lenz_law.time]
print("\n For Power = {} {}, water volume = {} {}, temperature_begin = {} {}, temperature_end = {} {}: cooling time = {} {}\n"
    .format(
        convert_to(P_example, units.watt).subs(units.watt, 1).evalf(5),
        units.watt,
        convert_to(v_example, units.liter).subs(units.liter, 1).evalf(3),
        units.liter,
        convert_to(t1_example - 273, units.kelvin).subs(units.kelvin, 1).evalf(3),
        Celsius,
        convert_to(t2_example - 273, units.kelvin).subs(units.kelvin, 1).evalf(3),
        Celsius,
        convert_to(solved_example, units.second).subs(units.second, 1).evalf(3),
        units.second
))