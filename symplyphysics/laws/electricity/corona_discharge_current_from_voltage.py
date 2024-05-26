from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The corona discharge is an independent discharge in a relatively dense gas. If an electric field is applied to two electrodes between
## which there is a gas gap, then a corona discharge occurs at a certain potential difference between the electrodes.

## Law is: I = A * k * U * (U - U0), where
## I - corona discharge current,
## A - experimentally determined gas coefficient,
## k - mobility of charged particles,
## U - voltage,
## U0 - corona discharge occurrence voltage.

current = Symbol("current", units.current)

gas_coefficient = Symbol("gas_coefficient", units.charge / (units.area * units.voltage))
mobility_of_charged_particles = Symbol("mobility_of_charged_particles", units.area / (units.voltage * units.time))
voltage = Symbol("voltage", units.voltage)
corona_discharge_occurrence_voltage = Symbol("corona_discharge_occurrence_voltage", units.voltage)

law = Eq(current, gas_coefficient * mobility_of_charged_particles * voltage * (voltage - corona_discharge_occurrence_voltage))


@validate_input(gas_coefficient_=gas_coefficient,
    mobility_of_charged_particles_=mobility_of_charged_particles,
    voltage_=voltage,
    corona_discharge_occurrence_voltage_=corona_discharge_occurrence_voltage)
@validate_output(current)
def calculate_current(gas_coefficient_: Quantity,
    mobility_of_charged_particles_: Quantity, voltage_: Quantity,
    corona_discharge_occurrence_voltage_: Quantity) -> Quantity:
    result_expr = solve(law, current,
        dict=True)[0][current]
    result_expr = result_expr.subs({
        gas_coefficient: gas_coefficient_,
        mobility_of_charged_particles: mobility_of_charged_particles_,
        voltage: voltage_,
        corona_discharge_occurrence_voltage: corona_discharge_occurrence_voltage_,
    })
    return Quantity(result_expr)
