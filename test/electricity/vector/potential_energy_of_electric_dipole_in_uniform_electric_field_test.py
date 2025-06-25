from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, units, errors
from symplyphysics.laws.electricity.vector import (
    potential_energy_of_electric_dipole_in_uniform_electric_field as law)

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector

Args = namedtuple("Args", "p e")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = QuantityCoordinateVector([
        1e-10 * units.coulomb * units.meter,
        0,
        -4e-10 * units.coulomb * units.meter,
    ], CARTESIAN)

    e = QuantityCoordinateVector([
        100 * units.volt / units.meter,
        -90 * units.volt / units.meter,
        -1 * units.volt / units.meter,
    ], CARTESIAN)

    return Args(p=p, e=e)


def test_law(test_args: Args) -> None:
    result = law.calculate_potential_energy(test_args.p, test_args.e)
    assert_equal(result, -1.04e-8 * units.joule)


def test_bad_electric_dipole_moment(test_args: Args) -> None:
    pb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_potential_energy(pb_vector, test_args.e)

    pb_scalar = units.meter * units.coulomb
    with raises(ValueError):
        law.calculate_potential_energy(pb_scalar, test_args.e)

    with raises(TypeError):
        law.calculate_potential_energy(100, test_args.e)
    with raises(TypeError):
        law.calculate_potential_energy([100], test_args.e)


def test_bad_electric_field(test_args: Args) -> None:
    eb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_potential_energy(test_args.p, eb_vector)

    eb_scalar = units.volt / units.meter
    with raises(ValueError):
        law.calculate_potential_energy(test_args.p, eb_scalar)

    with raises(TypeError):
        law.calculate_potential_energy(test_args.p, 100)
    with raises(TypeError):
        law.calculate_potential_energy(test_args.p, [100])
