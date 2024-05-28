from sympy import Eq, solve, pi, exp, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol
)
from sympy.physics.units import elementary_charge, electron_rest_mass, boltzmann_constant

# Description
## A gas–discharge plasma is an ionized gas in which the concentrations of positively and negatively charged particles are approximately
## equal to each other, and the Debye shielding radius is significantly smaller than the characteristic size of the volume in which the
## ionized gas is located.
## A probe is an auxiliary metal electrode that is inserted into a plasma volume for its examination. Flat, cylindrical and spherical probes
## are used. A potential is set relative to one of the electrodes on the probe and the dependence of the incoming current on the value of
## this potential is removed.

## Law is: I = 0.25 * S * e * n * sqrt(8 * k * T / (pi * m)) * exp(-e * (Uf - U) / (k * T)), where
## I - probe current,
## S - probe surface area,
## e - elementary charge,
## n - concentration of electrons in plasma,
## k - boltzmann constant,
## T - temperature of electrons in plasma,
## m - electron rest mass,
## Uf - floating plasma potential,
## U - potential at the location of the probe.

current = Symbol("current", units.current)

area_probe_surface = Symbol("area_probe_surface", units.area)
electron_concentration = Symbol("electron_concentration", 1 / units.volume)
electron_temperature = clone_symbol(symbols.thermodynamics.temperature, "electron_temperature")
floating_plasma_potential = Symbol("floating_plasma_potential", units.voltage)
potential_at_location_of_probe = Symbol("potential_at_location_of_probe", units.voltage)

expression_1 = 0.25 * area_probe_surface * elementary_charge * electron_concentration
expression_2 = sqrt(8 * boltzmann_constant * electron_temperature / (pi * electron_rest_mass))
expression_3 = exp(-elementary_charge * (floating_plasma_potential - potential_at_location_of_probe) / (boltzmann_constant * electron_temperature))

law = Eq(current, expression_1 * expression_2 * expression_3)


@validate_input(area_probe_surface_=area_probe_surface,
    electron_concentration_=electron_concentration,
    electron_temperature_=electron_temperature,
    floating_plasma_potential_=floating_plasma_potential,
    potential_at_location_of_probe_=potential_at_location_of_probe)
@validate_output(current)
def calculate_current(area_probe_surface_: Quantity,
    electron_concentration_: Quantity, electron_temperature_: Quantity,
    floating_plasma_potential_: Quantity, potential_at_location_of_probe_: Quantity) -> Quantity:
    result_expr = solve(law, current,
        dict=True)[0][current]
    result_expr = result_expr.subs({
        area_probe_surface: area_probe_surface_,
        electron_concentration: electron_concentration_,
        electron_temperature: electron_temperature_,
        floating_plasma_potential: floating_plasma_potential_,
        potential_at_location_of_probe: potential_at_location_of_probe_,
    })
    return Quantity(result_expr)
