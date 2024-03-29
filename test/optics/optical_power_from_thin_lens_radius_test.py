from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.optics import optical_power_from_thin_lens_radius as optical_power

# For a glass (n=1.5) biconvex lens with radii R1 = R2 = 50 cm in air (n0=1), the optical power will be 2 diopters.
# Note: R1 > 0, R2 < 0, as the lens is biconvex.

Args = namedtuple("Args", ["lens_index", "medium_index", "front_radius", "back_radius"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    lens_index = 1.5
    medium_index = 1
    front_radius = Quantity(50 * units.cm)
    back_radius = Quantity(-50 * units.cm)
    return Args(lens_index=lens_index,
        medium_index=medium_index,
        front_radius=front_radius,
        back_radius=back_radius)


def test_basic_power(test_args: Args) -> None:
    result = optical_power.calculate_optical_power(test_args.lens_index, test_args.medium_index,
        test_args.front_radius, test_args.back_radius)
    assert_equal(result, 2 * units.dioptre)


def test_bad_radius(test_args: Args) -> None:
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        optical_power.calculate_optical_power(test_args.lens_index, test_args.medium_index, rb,
            test_args.back_radius)
    with raises(TypeError):
        optical_power.calculate_optical_power(test_args.lens_index, test_args.medium_index, 100,
            test_args.back_radius)
    with raises(errors.UnitsError):
        optical_power.calculate_optical_power(test_args.lens_index, test_args.medium_index,
            test_args.front_radius, rb)
    with raises(TypeError):
        optical_power.calculate_optical_power(test_args.lens_index, test_args.medium_index,
            test_args.front_radius, 100)
