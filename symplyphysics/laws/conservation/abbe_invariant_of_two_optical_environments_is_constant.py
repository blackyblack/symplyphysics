from sympy import Eq, solve
from symplyphysics import (Symbol, units, print_expression, Quantity,
    validate_input, validate_output)

# Description
## The Abbe's invariant connects the front and back segments S and S',
## allowing one of them to be determined if the second one is known

# Law: Q = Q'
# Where:
## Q -  the Abbe's invariant of the medium before refraction
## Q' - the Abbe invariant of the refractive surface

# Conditions:
## - Abbe's formula is valid only for paraxial rays;
## - All rays emanating from point S and forming different but necessarily small angles with the axis will pass through the same point S' after refraction;
## - For real systems of only the paraxial region, formulas and positions that are valid for an ideal optical system can be applied.
## - The segment S does not depend on the angle, i.e. the homocentric beam of paraxial rays remains homocentric after passing through the refractive surface.

# NOTE:
## proofs: https://studme.org/341451/matematika_himiya_fizik/prelomlenie_otrazhenie_sveta_sfericheskoy_poverhnosti

invariant_abbe_before = Symbol("invariant_abbe_before", 1 / units.length)
invariant_abbe_after = Symbol("invariant_abbe_after", 1 / units.length)

law = Eq(invariant_abbe_before, invariant_abbe_after)


def print_law() -> str:
    return print_expression(law)


@validate_input(invariant_abbe_medium_=invariant_abbe_before)
@validate_output(invariant_abbe_after)
def calculate_surface_abbe_invariant(invariant_abbe_before_: Quantity) -> Quantity:
    solved = solve(law, invariant_abbe_after, dict=True)[0][invariant_abbe_after]
    result_expr = solved.subs({
        invariant_abbe_before: invariant_abbe_before_
    })
    result = Quantity(result_expr)
    return result
