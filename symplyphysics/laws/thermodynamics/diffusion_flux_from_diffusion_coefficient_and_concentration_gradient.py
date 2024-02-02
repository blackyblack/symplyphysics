from sympy import (Eq, solve, Derivative)
from symplyphysics import (units, Function, Quantity, Symbol, print_expression, validate_input,
    validate_output)


# Description
## Fick's first law relates the diffusive flux to the gradient of the concentration. It postulates that the flux goes from regions of high concentration to regions of low concentration,
## with a magnitude that is proportional to the concentration gradient (spatial derivative), or in simplistic terms the concept that a solute will move from a region of high concentration
## to a region of low concentration across a concentration gradient.
##The sign "-" in Fick's law indicates that the directions of the diffusion flow and the gradient are opposite: the gas diffuses towards a lower concentration, and the gradient is directed towards a higher concentration.

## Law: J = -D * dn / dx
## Where:
## J is the diffusion flux
## D is the diffusion coefficient
## n is the concentration, or amount of the substance per unit volume
## x is position

## Conditions
## Mixture is ideal

diffusion_flux = Function("diffusion_flux", units.amount_of_substance / (units.area * units.time))
diffusion_coefficient = Symbol("diffusion_coefficient", units.area / units.time)
concentration = Function("concentration", units.amount_of_substance / units.volume)
position = Symbol("position", units.length)

law = Eq(diffusion_flux(position), -diffusion_coefficient * Derivative(concentration(position), position))


def print_law() -> str:
    return print_expression(law)


@validate_input(diffusion_coefficient_=diffusion_coefficient, concentration_start_=concentration, concentration_end_=concentration, position_=position)
@validate_output(diffusion_flux)
def calculate_diffusion_flux(diffusion_coefficient_: Quantity, concentration_start_: Quantity,  concentration_end_: Quantity, position_: Quantity) -> Quantity:
    diffusion_flux_function_ = position * (concentration_end_ - concentration_start_) / position_
    applied_definition = law.subs({
        concentration(position): diffusion_flux_function_,
        diffusion_coefficient: diffusion_coefficient_
    })
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
