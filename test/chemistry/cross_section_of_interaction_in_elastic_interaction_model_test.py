from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import cross_section_of_interaction_in_elastic_interaction_model as cross_section_law

# Description
## The diameter of the atom is 1.31e-10 meter. Sutherland's constant is 133.5 kelvin. The temperature is 573 kelvin.
## Then the cross-sectional area of the interaction is 6.65e-20 meter^2.

Args = namedtuple("Args", ["diameter_of_atom", "constant_sutherland", "temperature"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    diameter_of_atom = Quantity(1.31e-10 * units.meter)
    constant_sutherland = Quantity(133.5 * units.kelvin)
    temperature = Quantity(573 * units.kelvin)

    return Args(diameter_of_atom=diameter_of_atom,
                constant_sutherland=constant_sutherland,
                temperature=temperature,)


def test_basic_cross_sectional_area_of_interaction(test_args: Args) -> None:
    result = cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.diameter_of_atom, test_args.constant_sutherland, test_args.temperature,)
    assert_equal(result, 6.65e-20 * units.meter**2)


def test_bad_diameter_of_atom(test_args: Args) -> None:
    diameter_of_atom = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(diameter_of_atom, test_args.constant_sutherland, test_args.temperature,)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(100, test_args.constant_sutherland, test_args.temperature,)


def test_bad_constant_sutherland(test_args: Args) -> None:
    constant_sutherland = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.diameter_of_atom, constant_sutherland, test_args.temperature,)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.diameter_of_atom, 100, test_args.temperature,)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.diameter_of_atom, test_args.constant_sutherland, temperature,)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.diameter_of_atom, test_args.constant_sutherland, 100,)
