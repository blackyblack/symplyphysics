from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.optics import angle_of_light_deflection_in_prism as angle_law

# Description
## Consider a prism with an angle between the faces equal to 45 degree (pi / 4 radian) and a refractive index of 2.5.
## The angle of deviation will be 67.5 degree (1.1781 radian).
## https://www.indigomath.ru//raschety/s3Jncy.html

Args = namedtuple("Args", ["angle_faces", "refractive_index"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    angle_faces = pi / 4
    refractive_index = 2.5
    return Args(angle_faces=angle_faces, refractive_index=refractive_index)


def test_basic_angle_deviation(test_args: Args) -> None:
    result = angle_law.calculate_angle_deviation(test_args.angle_faces, test_args.refractive_index)
    assert_equal(result, 1.1781)


def test_bad_angle_faces(test_args: Args) -> None:
    angle_faces = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        angle_law.calculate_angle_deviation(angle_faces, test_args.refractive_index)


def test_bad_refractive_index(test_args: Args) -> None:
    refractive_index = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        angle_law.calculate_angle_deviation(test_args.angle_faces, refractive_index)
