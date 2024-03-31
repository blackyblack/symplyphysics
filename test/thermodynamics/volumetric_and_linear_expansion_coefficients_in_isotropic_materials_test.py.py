from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    volumetric_and_linear_expansion_coefficients_in_isotropic_materials as coefficients_law,
)

# Description
## The volumetric coefficient of thermal expansion of a material with linear coefficient of thermal
## expansion alpha = 0.2 1/K is beta = 0.6 1/K.

Args = namedtuple("Args", "a")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(0.2 / units.kelvin)
    return Args(a=a)


def test_law(test_args: Args) -> None:
    result = coefficients_law.calculate_volumetric_expansion_coefficient(test_args.a)
    assert_equal(result, 0.6 / units.kelvin)


def test_bad_linear_coefficient() -> None:
    ab = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        coefficients_law.calculate_volumetric_expansion_coefficient(ab)
    with raises(TypeError):
        coefficients_law.calculate_volumetric_expansion_coefficient(100)
