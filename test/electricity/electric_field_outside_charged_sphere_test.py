from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import electric_field_outside_charged_sphere as intensity_law

# Description
## It is known that when the sphere is charged with 1 coulomb and distance from
## sphere is 1 meter, the electric field intensity is 8.992e9 volt / meter.
## https://vrcacademy.com/calculator/electric-field-uniformly-charged-sphere-calculator/


@fixture(name="test_args")
def test_args_fixture():
    charge = Quantity(1 * units.coulomb)
    distance = Quantity(1 * units.meter)

    Args = namedtuple("Args", ["charge", "distance"])
    return Args(charge=charge, distance=distance)


def test_basic_electric_intensity(test_args):
    result = intensity_law.calculate_electric_intensity(test_args.charge, test_args.distance)
    assert_equal(result, 8.992e9 * units.volt / units.meter)


def test_bad_charge(test_args):
    charge = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        intensity_law.calculate_electric_intensity(charge, test_args.distance)
    with raises(TypeError):
        intensity_law.calculate_electric_intensity(100, test_args.distance)


def test_bad_distance(test_args):
    distance = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        intensity_law.calculate_electric_intensity(test_args.charge, distance)
    with raises(TypeError):
        intensity_law.calculate_electric_intensity(test_args.charge, 100)
