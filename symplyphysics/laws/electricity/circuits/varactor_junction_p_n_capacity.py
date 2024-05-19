from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output,
                           dimensionless,)

# Description
## The principle of operation of the variator is based on the dependence of the barrier capacitance of the p-n junction on the voltage value.
## Knowing the doping coefficient, the material parameter and the capacity without bias voltage, it is possible to calculate the capacity
## of the varactor at a given voltage.

## Law is: C = C0 / (1 - V / V0)^y, where
## C - barrier capacity of the p-n junction,
## C0 - barrier capacitance of the p-n junction without bias voltage,
## V - voltage,
## V0 - material parameter (a certain constant having the dimension of the voltage and depending on the material),
## y - doping coefficient.

junction_capacitance = Symbol("junction_capacitance", units.capacitance)

junction_capacitance_without_bias_voltage = Symbol("junction_capacitance_without_bias_voltage", units.capacitance)
voltage = Symbol("voltage", units.voltage)
material_parameter = Symbol("material_parameter", units.voltage)
doping_coefficient = Symbol("doping_coefficient", dimensionless)

law = Eq(junction_capacitance, junction_capacitance_without_bias_voltage / (1 - voltage / material_parameter)**doping_coefficient)


@validate_input(junction_capacitance_without_bias_voltage_=junction_capacitance_without_bias_voltage,
                voltage_=voltage,
                material_parameter_=material_parameter,
                doping_coefficient_=doping_coefficient)
@validate_output(junction_capacitance)
def calculate_junction_capacitance(junction_capacitance_without_bias_voltage_: Quantity, voltage_: Quantity,
                                                  material_parameter_: Quantity, doping_coefficient_: float) -> Quantity:
    result_expr = solve(law, junction_capacitance, dict=True)[0][junction_capacitance]
    result_expr = result_expr.subs({
        junction_capacitance_without_bias_voltage: junction_capacitance_without_bias_voltage_,
        voltage: voltage_,
        material_parameter: material_parameter_,
        doping_coefficient: doping_coefficient_,
    })
    return Quantity(result_expr)
