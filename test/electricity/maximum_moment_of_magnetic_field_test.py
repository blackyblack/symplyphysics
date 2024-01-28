from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.electricity import maximum_moment_of_magnetic_field as moment_law

# Description
## With a current value of 0.5 [ampere] and a contour area of 0.02 [meter^2],
## the magnetic moment is 0.01 [ampere * meter^2].
## https://online.mephi.ru/courses/physics/electricity/data/course/4/4.1.html


@fixture(name="test_args")
def test_args_fixture():
    current = Quantity(0.5 * units.ampere)
    area = Quantity(0.02 * units.meter**2)

    Args = namedtuple("Args", ["current", "area"])
    return Args(current=current, area=area)


def test_basic_moment(test_args):
    result = moment_law.calculate_moment(test_args.current, test_args.area)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current * units.area)
    result = convert_to(result, units.ampere * units.meter**2).evalf(5)
    assert result == approx(0.01, rel=0.01)


def test_bad_current(test_args):
    current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        moment_law.calculate_moment(current, test_args.area)
    with raises(TypeError):
        moment_law.calculate_moment(100, test_args.area)


def test_bad_area(test_args):
    area = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        moment_law.calculate_moment(test_args.current, area)
    with raises(TypeError):
        moment_law.calculate_moment(test_args.current, 100)
