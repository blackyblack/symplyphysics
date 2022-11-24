from symplyphysics import (
    units, convert_to,solve, SI,Quantity, expr_to_quantity
)
from symplyphysics import (
    symbols, Eq, pretty
)
from symplyphysics.laws.electricity import amount_energy_from_voltage_time_resistance as joule_lenz_law
from symplyphysics.definitions import power_is_proportional_voltage_and_current as operate_power
from symplyphysics.definitions import amount_energy_from_mass_and_temperature as operate_energy
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohm_law

# The household electric kettle heated 0.5 liters of water to boiling.
# The power of the kettle is 1500 watts. How long will it take for the
# temperature of the water to drop to a comfortable 35 degrees?
P_example = units.Quantity('P_example')
SI.set_quantity_dimension(P_example, units.power)
SI.set_quantity_scale_factor(P_example, 1500 * units.watt)
m_example = units.Quantity('m_example')
SI.set_quantity_dimension(m_example, units.mass)
SI.set_quantity_scale_factor(m_example, 0.5 * units.kilogram)
C_example = units.Quantity('C_example')
SI.set_quantity_dimension(C_example, units.energy / (units.mass * units.temperature))
SI.set_quantity_scale_factor(C_example, 4200 * units.joule / (units.kilogram * units.kelvin))
t1_example = units.Quantity('t1_example')
SI.set_quantity_dimension(t1_example, units.temperature)
SI.set_quantity_scale_factor(t1_example, 373 * units.kelvin)
t2_example = units.Quantity('t2_example')
SI.set_quantity_dimension(t2_example, units.temperature)
SI.set_quantity_scale_factor(t2_example, 308 * units.kelvin)
P, t, m , C, t1, t2 = symbols('P t m C t2 t1')
Q1 = joule_lenz_law.law.subs({
    joule_lenz_law.operating_voltage / joule_lenz_law.resistance_heater: ohm_law.current,
    ohm_law.current * joule_lenz_law.operating_voltage: P, joule_lenz_law.operating_time: t
})
Q2 = operate_energy.law.subs({operate_energy.body_mass: m, operate_energy.specific_heat_capacity: C,
     operate_energy.final_temperature: t2, operate_energy.initial_temperature: t1
})
law = Eq(Q1, Q2)
solved = solve(law, t, dict=True)[0][t]
answer = Eq(t, solved)
print("\nFormula is:\n\n {}".format(pretty(answer, use_unicode=False)))
law_example = Eq(Q1, Q2).subs({
    P: convert_to(P_example, units.watt).subs(units.watt, 1).evalf(5),
    m: convert_to(m_example, units.kilogram).subs(units.kilogram, 1).evalf(4),
    C: convert_to(C_example, units.joule / (units.kilogram * units.kelvin)).subs(units.joule / (units.kilogram * units.kelvin), 1).evalf(6),
    t2: convert_to(t1_example, units.kelvin).subs(units.kelvin, 1).evalf(3),
    t1: convert_to(t2_example, units.kelvin).subs(units.kelvin, 1).evalf(3)
})
solved_example = solve(law_example, t, dict=True)[0][t]
print("\n For P = {} {}, water mass = {} {}, t2 = {} {}, t1 = {} {}, cooling time = {} {}\n"
    .format(
        convert_to(P_example, units.watt).subs(units.watt, 1).evalf(5),
        units.watt,
        convert_to(m_example, units.kilogram).subs(units.kilogram, 1).evalf(4),
        units.kilogram,
        convert_to(t1_example, units.kelvin).subs(units.kelvin, 1).evalf(3),
        units.kelvin,
        convert_to(t2_example, units.kelvin).subs(units.kelvin, 1).evalf(3),
        units.kelvin,
        convert_to(solved_example, units.second).subs(units.second, 1).evalf(5),
        units.second
))