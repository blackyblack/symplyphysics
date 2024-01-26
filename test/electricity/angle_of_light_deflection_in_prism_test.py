from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors, dimensionless,)
from symplyphysics.laws.electricity import angle_of_light_deflection_in_prism as angle_law

# Description
## Consider a prism with an angle between the faces equal to 45 degree and a refractive index of 15.
## The angle of deviation will be 630 degrees.
## https://www.indigomath.ru//raschety/s3Jncy.html


@fixture(name="test_args")
def test_args_fixture():
    angle_faces = 45
    refractive_index = 15

    Args = namedtuple("Args", ["angle_faces", "refractive_index"])
    return Args(angle_faces=angle_faces, refractive_index=refractive_index)


def test_basic_angle_deviation(test_args):
    result = angle_law.calculate_angle_deviation(test_args.angle_faces, test_args.refractive_index)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, dimensionless)
    result = convert_to(result, dimensionless).evalf(5)
    assert result == approx(630, rel=0.01)


def test_bad_angle_faces(test_args):
    angle_faces = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        angle_law.calculate_angle_deviation(angle_faces, test_args.refractive_index)
    with raises(AttributeError):
        angle_law.calculate_angle_deviation(True, test_args.refractive_index)


def test_bad_refractive_index(test_args):
    refractive_index = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        angle_law.calculate_angle_deviation(test_args.angle_faces, refractive_index)
    with raises(AttributeError):
        angle_law.calculate_angle_deviation(test_args.angle_faces, True)
