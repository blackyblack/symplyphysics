from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)

from symplyphysics.laws.thermodynamics import number_of_impacts_of_molecules_on_the_wall_from_area_and_velocity as number_of_impacts


@fixture(name="test_args")
def test_args_fixture():
    molecules_concentration = Quantity(1e25 * (1 / units.meter**3))
    area = Quantity(0.5 * units.meter**2)
    velocity_projection = Quantity(500 * units.meter / units.second)
    time = Quantity(5 * units.second)
    Args = namedtuple("Args", ["molecules_concentration", "area", "velocity_projection", "time"])
    return Args(molecules_concentration=molecules_concentration,
        area=area,
        velocity_projection=velocity_projection,
        time=time)


def test_basic_number_of_impacts(test_args):
    result = number_of_impacts.calculate_number_of_impacts(test_args.molecules_concentration,
        test_args.area, test_args.velocity_projection, test_args.time)
    assert_equal(result, 6.25e27)


def test_bad_molecules_concentration(test_args):
    bad_molecules_concentration = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        number_of_impacts.calculate_number_of_impacts(bad_molecules_concentration, test_args.area,
            test_args.velocity_projection, test_args.time)
    with raises(TypeError):
        number_of_impacts.calculate_number_of_impacts(100, test_args.area,
            test_args.velocity_projection, test_args.time)


def test_bad_area(test_args):
    bad_area = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        number_of_impacts.calculate_number_of_impacts(test_args.molecules_concentration, bad_area,
            test_args.velocity_projection, test_args.time)
    with raises(TypeError):
        number_of_impacts.calculate_number_of_impacts(test_args.molecules_concentration, 100,
            test_args.velocity_projection, test_args.time)


def test_bad_velocity_projection(test_args):
    bad_velocity_projection = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        number_of_impacts.calculate_number_of_impacts(test_args.molecules_concentration,
            test_args.area, bad_velocity_projection, test_args.time)
    with raises(TypeError):
        number_of_impacts.calculate_number_of_impacts(test_args.molecules_concentration,
            test_args.area, 100, test_args.time)


def test_bad_time(test_args):
    bad_time = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        number_of_impacts.calculate_number_of_impacts(test_args.molecules_concentration,
            test_args.area, test_args.velocity_projection, bad_time)
    with raises(TypeError):
        number_of_impacts.calculate_number_of_impacts(test_args.molecules_concentration,
            test_args.area, test_args.velocity_projection, 100)
