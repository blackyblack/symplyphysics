from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity import free_fall_acceleration_from_height as free_fall_acceleration

# G - universal gravity constant  6.672e-11 N*m^2/kg^2
# M - Earth mass constant         5.9722e24 kg
# R - Earth radius constant       6371 km

Args = namedtuple("Args", ["earth_mass", "earth_radius", "height_from_surface"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    earth_mass = Quantity(5.9722e24 * units.kilogram)
    earth_radius = Quantity(6371 * units.kilometer)
    height_from_surface = Quantity(1 * units.meter)
    return Args(earth_mass=earth_mass,
        earth_radius=earth_radius,
        height_from_surface=height_from_surface)


def test_basic_acceleration(test_args: Args) -> None:
    result = free_fall_acceleration.calculate_acceleration(test_args.earth_mass,
        test_args.earth_radius, test_args.height_from_surface)
    assert_equal(result, 9.823 * units.meter / units.second**2)


def test_bad_earth_mass(test_args: Args) -> None:
    emb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        free_fall_acceleration.calculate_acceleration(emb, test_args.earth_radius,
            test_args.height_from_surface)
    with raises(TypeError):
        free_fall_acceleration.calculate_acceleration(100, test_args.earth_radius,
            test_args.height_from_surface)


def test_bad_earth_radius(test_args: Args) -> None:
    erb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        free_fall_acceleration.calculate_acceleration(test_args.earth_mass, erb,
            test_args.height_from_surface)
    with raises(TypeError):
        free_fall_acceleration.calculate_acceleration(test_args.earth_mass, 100,
            test_args.height_from_surface)


def test_bad_height(test_args: Args) -> None:
    hb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        free_fall_acceleration.calculate_acceleration(test_args.earth_mass, test_args.earth_radius,
            hb)
    with raises(TypeError):
        free_fall_acceleration.calculate_acceleration(test_args.earth_mass, test_args.earth_radius,
            100)
