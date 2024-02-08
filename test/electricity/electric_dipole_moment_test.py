from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import electric_dipole_moment as moment_law

# Description
## It is known that with a charge of 1 coulomb and a distance of 2 meters,
## dipole moment will be equal to 2 coulomb * meters.
## https://www.calculatoratoz.com/ru/electric-dipole-moment-calculator/Calc-579


@fixture(name="test_args")
def test_args_fixture():
    charge = Quantity(1 * units.coulomb)
    distance = Quantity(2 * units.meter)

    Args = namedtuple("Args", ["charge", "distance"])
    return Args(charge=charge, distance=distance)


def test_basic_electric_dipole_moment(test_args):
    result = moment_law.calculate_electric_moment(test_args.charge, test_args.distance)
    assert_equal(result, 2 * units.coulomb * units.meter)


def test_bad_charge(test_args):
    charge = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        moment_law.calculate_electric_moment(charge, test_args.distance)
    with raises(TypeError):
        moment_law.calculate_electric_moment(100, test_args.distance)


def test_bad_distance(test_args):
    distance = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        moment_law.calculate_electric_moment(test_args.charge, distance)
    with raises(TypeError):
        moment_law.calculate_electric_moment(test_args.charge, 100)
