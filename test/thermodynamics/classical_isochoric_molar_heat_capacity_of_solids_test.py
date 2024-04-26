from symplyphysics import (
    assert_equal,
    units,
)
from symplyphysics.laws.thermodynamics import (
    classical_isochoric_molar_heat_capacity_of_solids as dulong_petit_law,
)

def test_law() -> None:
    result = dulong_petit_law.calculate_isochoric_molar_heat_capacity()
    assert_equal(result, 24.9 * units.joule / (units.kelvin * units.mole), tolerance=2e-3)
