from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes,)
from symplyphysics.laws.electricity.circuits import resonant_frequency_of_rectangular_resonator as frequency_law

# Description
## The height, width and length of the resonator are 2 centimeter, respectively. The dielectric constant and magnetic permeability of the material
## filling the resonator are 2.2. and 1, respectively. The first, second and third indexes are equal to 1.
## Then the resonant frequency is 8.752 gigahertz.
## ## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", [
    "first_index", "second_index", "third_index", "resonator_width",
    "resonator_height", "resonator_length", "relative_permittivity", "relative_permeability"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    first_index = 1
    second_index = 1
    third_index = 1
    resonator_width = Quantity(2 * units.centimeter)
    resonator_height = Quantity(2 * units.centimeter)
    resonator_length = Quantity(2 * units.centimeter)
    relative_permittivity = 2.2
    relative_permeability = 1

    return Args(first_index=first_index,
        second_index=second_index,
        third_index=third_index,
        resonator_width=resonator_width,
        resonator_height=resonator_height,
        resonator_length=resonator_length,
        relative_permittivity=relative_permittivity,
        relative_permeability=relative_permeability
        )


def test_basic_resonant_frequency(test_args: Args) -> None:
    result = frequency_law.calculate_resonant_frequency(test_args.first_index,
        test_args.second_index, test_args.third_index, test_args.resonator_width,
        test_args.resonator_height, test_args.resonator_length, test_args.relative_permittivity, test_args.relative_permeability)
    assert_equal(result, 8.752 * prefixes.giga * units.hertz)


def test_bad_index(test_args: Args) -> None:
    bad_index = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_resonant_frequency(bad_index,
            test_args.second_index, test_args.third_index, test_args.resonator_width,
            test_args.resonator_height, test_args.resonator_length, test_args.relative_permittivity, test_args.relative_permeability)
    with raises(errors.UnitsError):
        frequency_law.calculate_resonant_frequency(test_args.first_index,
            bad_index, test_args.third_index, test_args.resonator_width,
            test_args.resonator_height, test_args.resonator_length, test_args.relative_permittivity, test_args.relative_permeability)
    with raises(errors.UnitsError):
        frequency_law.calculate_resonant_frequency(test_args.first_index,
            test_args.second_index, bad_index, test_args.resonator_width,
            test_args.resonator_height, test_args.resonator_length, test_args.relative_permittivity, test_args.relative_permeability)


def test_bad_resonator_width(test_args: Args) -> None:
    resonator_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_resonant_frequency(test_args.first_index, test_args.second_index, test_args.third_index, resonator_width, test_args.resonator_height, test_args.resonator_length, test_args.relative_permittivity, test_args.relative_permeability)
    with raises(TypeError):
        frequency_law.calculate_resonant_frequency(test_args.first_index, test_args.second_index, test_args.third_index, 100, test_args.resonator_height, test_args.resonator_length, test_args.relative_permittivity, test_args.relative_permeability)


def test_bad_resonator_height(test_args: Args) -> None:
    resonator_height = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_resonant_frequency(test_args.first_index, test_args.second_index, test_args.third_index, test_args.resonator_width, resonator_height, test_args.resonator_length, test_args.relative_permittivity, test_args.relative_permeability)
    with raises(TypeError):
        frequency_law.calculate_resonant_frequency(test_args.first_index, test_args.second_index, test_args.third_index, test_args.resonator_width, 100, test_args.resonator_length, test_args.relative_permittivity, test_args.relative_permeability)


def test_bad_resonator_length(test_args: Args) -> None:
    resonator_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_resonant_frequency(test_args.first_index,
            test_args.second_index, test_args.third_index, test_args.resonator_width,
            test_args.resonator_height, resonator_length, test_args.relative_permittivity, test_args.relative_permeability)
    with raises(TypeError):
        frequency_law.calculate_resonant_frequency(test_args.first_index,
            test_args.second_index, test_args.third_index, test_args.resonator_width,
            test_args.resonator_height, 100, test_args.relative_permittivity, test_args.relative_permeability)


def test_bad_relative_permittivity(test_args: Args) -> None:
    relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_resonant_frequency(test_args.first_index, test_args.second_index, test_args.third_index, test_args.resonator_width, test_args.resonator_height, test_args.resonator_length, relative_permittivity, test_args.relative_permeability)


def test_bad_relative_permeability(test_args: Args) -> None:
    relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_resonant_frequency(test_args.first_index, test_args.second_index, test_args.third_index, test_args.resonator_width, test_args.resonator_height, test_args.resonator_length, test_args.relative_permittivity, relative_permeability)
