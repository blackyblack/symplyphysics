from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors, Quantity
from symplyphysics.laws.electricity.vector import (
    electric_dipole_moment_of_electrically_neutral_system as law)

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector, CARTESIAN
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "qs rs")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q1 = Quantity(1e-10 * units.coulomb)
    q2 = Quantity(-7e-11 * units.coulomb)
    q3 = Quantity(-3e-11 * units.coulomb)
    q4 = Quantity(5e-10 * units.coulomb)

    r1 = QuantityCoordinateVector([1.0e-3 * units.meter, 0, 0], CARTESIAN)
    r2 = QuantityCoordinateVector([0, -2.0e-3 * units.meter, 0], CARTESIAN)
    r3 = QuantityCoordinateVector([
        1.2e-3 * units.meter,
        -2.0e-2 * units.meter,
        3.0e-3 * units.meter,
    ], CARTESIAN)
    r4 = QuantityCoordinateVector([0, 0, 4.2e-2 * units.meter], CARTESIAN)

    return Args(qs=(q1, q2, q3, q4), rs=(r1, r2, r3, r4))


def test_law(test_args: Args) -> None:
    result = law.calculate_electric_dipole_moment(test_args.qs[:-1], test_args.rs[:-1])

    expected = QuantityCoordinateVector([
        6.4e-14 * units.coulomb * units.meter,
        7.4e-13 * units.coulomb * units.meter,
        -9.0e-14 * units.coulomb * units.meter,
    ], CARTESIAN)

    assert_equal_vectors(expected, result)


def test_bad_inputs(test_args: Args) -> None:
    # input lengths don't match

    with raises(ValueError):
        law.calculate_electric_dipole_moment(test_args.qs, test_args.rs[:-1])

    with raises(ValueError):
        law.calculate_electric_dipole_moment(test_args.qs[:-1], test_args.rs)

    # the system is not electrically neutral

    with raises(AssertionError):
        law.calculate_electric_dipole_moment(test_args.qs, test_args.rs)


def test_bad_charge(test_args: Args) -> None:
    qb = Quantity(units.candela)

    with raises(errors.UnitsError):
        law.calculate_electric_dipole_moment((qb, test_args.qs[0]), test_args.rs[:2])
    with raises(TypeError):
        law.calculate_electric_dipole_moment((100, test_args.qs[0]), test_args.rs[:2])

    with raises(errors.UnitsError):
        law.calculate_electric_dipole_moment((test_args.qs[0], qb), test_args.rs[:2])
    with raises(TypeError):
        law.calculate_electric_dipole_moment((test_args.qs[0], 100), test_args.rs[:2])

    with raises(errors.UnitsError):
        law.calculate_electric_dipole_moment(qb, test_args.rs)
    with raises(TypeError):
        law.calculate_electric_dipole_moment(100, test_args.rs)


def test_bad_position_vector(test_args: Args) -> None:
    rb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    rb_scalar = units.meter

    with raises(errors.UnitsError):
        law.calculate_electric_dipole_moment(test_args.qs[:-1], (rb_vector,) + test_args.rs[:2])
    with raises(ValueError):
        law.calculate_electric_dipole_moment(test_args.qs[:-1], (rb_scalar,) + test_args.rs[:2])
    with raises(TypeError):
        law.calculate_electric_dipole_moment(test_args.qs[:-1], (100,) + test_args.rs[:2])

    with raises(errors.UnitsError):
        law.calculate_electric_dipole_moment(test_args.qs[:-1], test_args.rs[:2] + (rb_vector,))
    with raises(ValueError):
        law.calculate_electric_dipole_moment(test_args.qs[:-1], test_args.rs[:2] + (rb_scalar,))
    with raises(TypeError):
        law.calculate_electric_dipole_moment(test_args.qs[:-1], test_args.rs[:2] + (100,))

    with raises(errors.UnitsError):
        law.calculate_electric_dipole_moment(test_args.qs, rb_vector)
    with raises(TypeError):
        law.calculate_electric_dipole_moment(test_args.qs, rb_scalar)
