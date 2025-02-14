from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, units, Quantity, errors
from symplyphysics.laws.gravity import relative_aperture_of_telescope as relative_aperture_law

# Description
## With a lens diameter of 37.5 centimeters and a focal length of 600 cm, the relative aperture
## of the telescope is 1/16.
## https://znanija.com/task/47215311

Args = namedtuple("Args", ["lens_diameter", "focal_length_lens"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    lens_diameter = Quantity(37.5 * units.centimeter)
    focal_length_lens = Quantity(600 * units.centimeter)

    return Args(lens_diameter=lens_diameter, focal_length_lens=focal_length_lens)


def test_basic_relative_aperture(test_args: Args) -> None:
    result = relative_aperture_law.calculate_relative_aperture(test_args.lens_diameter,
        test_args.focal_length_lens)
    assert_equal(result, 1 / 16)


def test_bad_lens_diameter(test_args: Args) -> None:
    lens_diameter = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relative_aperture_law.calculate_relative_aperture(lens_diameter,
            test_args.focal_length_lens)
    with raises(TypeError):
        relative_aperture_law.calculate_relative_aperture(100, test_args.focal_length_lens)


def test_bad_focal_length_lens(test_args: Args) -> None:
    focal_length_lens = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relative_aperture_law.calculate_relative_aperture(test_args.lens_diameter,
            focal_length_lens)
    with raises(TypeError):
        relative_aperture_law.calculate_relative_aperture(test_args.lens_diameter, 100)
