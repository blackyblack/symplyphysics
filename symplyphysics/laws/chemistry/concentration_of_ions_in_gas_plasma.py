from sympy import Eq, solve, sqrt, pi
from sympy.physics.units import elementary_charge, electron_rest_mass, boltzmann_constant
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, clone_symbol, symbols)

# Description
## A gas–discharge plasma is an ionized gas in which the concentrations of positively and negatively charged particles are approximately
## equal to each other, and the Debye shielding radius is significantly smaller than the characteristic size of the volume in which the
## ionized gas is located.
## A probe is an auxiliary metal electrode that is inserted into a plasma volume for its examination. Flat, cylindrical and spherical probes
## are used. A potential is set relative to one of the electrodes on the probe and the dependence of the incoming current on the value of
## this potential is removed.
## By measuring the probe current, the concentration of ions in the plasma can be determined.

## Law is: n = I0 / (S * e * sqrt(k * T / (2 * pi * m))), where
## n - concentration of ions in plasma,
## I - probe current (at a probe potential equal to the plasma potential),
## S - probe surface area,
## e - elementary charge,
## k - boltzmann constant,
## T - temperature of electrons in plasma,
## m - electron rest mass.

ion_concentration = Symbol("ion_concentration", 1 / units.volume)

probe_current = Symbol("probe_current", units.current)
electron_temperature = clone_symbol(symbols.thermodynamics.temperature, "electron_temperature")
area_probe_surface = Symbol("area_probe_surface", units.area)

law = Eq(ion_concentration, probe_current / (area_probe_surface * elementary_charge * sqrt(boltzmann_constant * electron_temperature / (2 * pi * electron_rest_mass))))


@validate_input(probe_current_=probe_current,
    electron_temperature_=electron_temperature,
    area_probe_surface_=area_probe_surface)
@validate_output(ion_concentration)
def calculate_ion_concentration(probe_current_: Quantity, electron_temperature_: Quantity,
    area_probe_surface_: Quantity) -> Quantity:
    result_expr = solve(law, ion_concentration, dict=True)[0][ion_concentration]
    result_expr = result_expr.subs({
        probe_current: probe_current_,
        electron_temperature: electron_temperature_,
        area_probe_surface: area_probe_surface_,
    })
    return Quantity(result_expr)
