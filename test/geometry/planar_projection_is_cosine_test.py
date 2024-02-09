from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law

Args = namedtuple("Args", ["force_vector_amplitude", "angle_between_vector_and_horizontal_axis"])


#We are having force vector of 3 Newtons with angle 60 degrees to horizontal axis. Projection of this vector to horizontal axis should be 2 times less than the vector.
@fixture(name="test_args")
def test_args_fixture() -> Args:
    force_vector_amplitude = Quantity(3 * units.newton)
    angle_between_vector_and_axis = Quantity(60 * units.degree)
    return Args(force_vector_amplitude=force_vector_amplitude,
        angle_between_vector_and_horizontal_axis=angle_between_vector_and_axis)


def test_basic_projection(test_args: Args) -> None:
    result = projection_law.calculate_projection(test_args.force_vector_amplitude,
        test_args.angle_between_vector_and_horizontal_axis)
    assert_equal(result, 1.5 * units.newton)


def test_projection_with_number(test_args: Args) -> None:
    result = projection_law.calculate_projection(test_args.force_vector_amplitude,
        (pi / 3).evalf(4))
    assert_equal(result, 1.5 * units.newton)


def test_bad_angle(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        projection_law.calculate_projection(test_args.force_vector_amplitude, ab)
