from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import kinetic_energy_of_rotating_body as kinetic_energy_law

# Description
## A rigid body of rotational inertia I = 0.02 kg*m^2 rotates around an axis with the angular
## speed of 3 rad/s. It kinetic energy should amount to 0.09 J.


@fixture(name="test_args")
def test_args_fixture():
    I = Quantity(0.02 * units.kilogram * units.meter**2)
    w = Quantity(3.0 * units.radian / units.second)
    Args = namedtuple("Args", "I w")
    return Args(I=I, w=w)


def test_basic_law(test_args):
    result = kinetic_energy_law.calculate_kinetic_energy(test_args.I, test_args.w)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).evalf(3)
    assert result_energy == approx(0.09, 1e-4)


def test_bad_rotational_inertia(test_args):
    Ib = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        kinetic_energy_law.calculate_kinetic_energy(Ib, test_args.w)
    with raises(TypeError):
        kinetic_energy_law.calculate_kinetic_energy(100, test_args.w)


def test_bad_angular_velocity(test_args):
    wb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        kinetic_energy_law.calculate_kinetic_energy(test_args.I, wb)
    with raises(TypeError):
        kinetic_energy_law.calculate_kinetic_energy(test_args.I, 100)
