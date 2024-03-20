from sympy import Eq, sqrt, pi
from symplyphysics import (
    Quantity,
    Symbol,
    print_expression,
    units,
    dimensionless,
    validate_input,
    validate_output,
)

# Description
## The mean free path of a molecule in random motion is its average path length between collisions.

# Law: lambda = 1 / (sqrt(2) * pi * d**2 * (N / V))
## lambda - mean free path estimate
## d - molecular diameter
## N - number of molecules found in volume V
## V - volume

mean_free_path = Symbol("mean_free_path", units.length, positive=True)
molecular_diameter = Symbol("molecular_diameter", units.length, positive=True)
molecule_count = Symbol("molecule_count", dimensionless, integer=True, positive=True)
volume = Symbol("volume", units.volume, positive=True)

law = Eq(mean_free_path, 1 / (sqrt(2) * pi * molecular_diameter**2 * (molecule_count / volume)))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    molecular_diameter_=molecular_diameter,
    molecule_count_=molecule_count,
    volume_=volume,
)
@validate_output(mean_free_path)
def calculate_mean_free_path(
    molecular_diameter_: Quantity,
    molecule_count_: int,
    volume_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        molecular_diameter: molecular_diameter_,
        molecule_count: molecule_count_,
        volume: volume_,
    })
    return Quantity(result)
