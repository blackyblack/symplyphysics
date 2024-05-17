from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    isochoric_molar_heat_capacity_of_ideal_gas_via_degrees_of_freedom as heat_capacity_law,)

# Description
## The isochoric molar heat capacity of an ideal gas consisting of molecules with 6 degrees of freedom
## is 24.9 J/(K*mol).

Args = namedtuple("Args", "f")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = 6
    return Args(f=f)


def test_law(test_args: Args) -> None:
    result = heat_capacity_law.calculate_isochoric_molar_heat_capacity(test_args.f)
    assert_equal(result, 24.9 * units.joule / (units.kelvin * units.mole), tolerance=2e-3)


def test_bad_degree_count() -> None:
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heat_capacity_law.calculate_isochoric_molar_heat_capacity(fb)
