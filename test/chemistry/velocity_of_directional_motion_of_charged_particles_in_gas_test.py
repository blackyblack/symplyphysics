from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import velocity_of_directional_motion_of_charged_particles_in_gas as velocity_law

# Description
## The mobility of charged particles at a single pressure is 2450 [meter^2 * pascal / (volt * second)].
## Pressure is 100 pascal. The electric intensity is 3000 volt per meter.
## Then the velocity of the directional motion of charged particles in a gas is 73.5e3 meter per second.

Args = namedtuple("Args", ["mobility_at_unit_pressure", "pressure", "electric_intensity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mobility_at_unit_pressure = Quantity(2450 * units.meter**2 * units.pascal / (units.volt * units.second))
    pressure = Quantity(100 * units.pascal)
    electric_intensity = Quantity(3000 * units.volt / units.meter)

    return Args(
        mobility_at_unit_pressure=mobility_at_unit_pressure,
        pressure=pressure,
        electric_intensity=electric_intensity,
    )


def test_basic_velocity(test_args: Args) -> None:
    result = velocity_law.calculate_velocity(
        test_args.mobility_at_unit_pressure,
        test_args.pressure,
        test_args.electric_intensity,
    )
    assert_equal(result, 73.5e3 * units.meter / units.second)


def test_bad_mobility_at_unit_pressure(test_args: Args) -> None:
    mobility_at_unit_pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_velocity(
            mobility_at_unit_pressure,
            test_args.pressure,
            test_args.electric_intensity,
        )
    with raises(TypeError):
        velocity_law.calculate_velocity(
            100,
            test_args.pressure,
            test_args.electric_intensity,
        )


def test_bad_pressure(test_args: Args) -> None:
    pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_velocity(
            test_args.mobility_at_unit_pressure,
            pressure,
            test_args.electric_intensity,
        )
    with raises(TypeError):
        velocity_law.calculate_velocity(
            test_args.mobility_at_unit_pressure,
            100,
            test_args.electric_intensity,
        )


def test_bad_electric_intensity(test_args: Args) -> None:
    electric_intensity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_velocity(
            test_args.mobility_at_unit_pressure,
            test_args.pressure,
            electric_intensity,
        )
    with raises(TypeError):
        velocity_law.calculate_velocity(
            test_args.mobility_at_unit_pressure,
            test_args.pressure,
            100,
        )
