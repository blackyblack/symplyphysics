from sympy import Eq
from symplyphysics import (
    Quantity,
    Symbol,
    units,
    dimensionless,
    validate_input,
    validate_output,
)

# Description
## The kinetic energy of an object is the form of energy that is possesses due to its motion.

# Law: E_kin = (gamma - 1) * m_rest * c**2
## E_kin - kinetic energy of body
## gamma - [Lorentz factor](../../definitions/lorentz_factor.py)
## m_rest - rest (invariant) mass of body
## c - speed of light

# Notes
## - The work expended accelerating an object from rest approaches infinity as the velocity approaches the
##   speed of light. Thus it is impossible to accelerate an object across this boundary.

kinetic_energy = Symbol("kinetic_energy", units.energy)
lorentz_factor = Symbol("lorentz_factor", dimensionless)
rest_mass = Symbol("rest_mass", units.mass)

law = Eq(kinetic_energy, (lorentz_factor - 1) * rest_mass * units.speed_of_light**2)

# TODO: derive from relativistic momentum and expression for kinetic energy


@validate_input(
    lorentz_factor_=lorentz_factor,
    rest_mass_=rest_mass,
)
@validate_output(kinetic_energy)
def calculate_kinetic_energy(
    lorentz_factor_: float,
    rest_mass_: Quantity,
) -> Quantity:
    if lorentz_factor_ < 1:
        raise ValueError("Lorentz factor must be greater or equal to 1")

    result = law.rhs.subs({
        lorentz_factor: lorentz_factor_,
        rest_mass: rest_mass_,
    })
    return Quantity(result)
