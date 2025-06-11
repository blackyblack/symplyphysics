from collections import namedtuple
from pytest import fixture, raises
from sympy.physics.units import speed_of_light
from symplyphysics import assert_equal, units, Quantity, errors
from symplyphysics.laws.relativistic.vector import energy_momentum_relation as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.solvers import solve_for_vector

Args = namedtuple("Args", "p e v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    momentum_unit = units.kilogram * speed_of_light

    p = QuantityCoordinateVector([
        Quantity(5.0e-5 * momentum_unit),
        Quantity(1e-4 * momentum_unit),
        Quantity(-5e-4 * momentum_unit),
    ], CARTESIAN)

    e = Quantity(1e-3 * units.kilogram * speed_of_light**2)

    v = QuantityCoordinateVector([
        Quantity(0.05 * speed_of_light),
        Quantity(0.1 * speed_of_light),
        Quantity(-0.5 * speed_of_light),
    ], CARTESIAN)

    return Args(p=p, e=e, v=v)


def test_momentum_law(test_args: Args) -> None:
    result = law.calculate_momentum(test_args.e, test_args.v)
    assert_equal_vectors(result, test_args.p)


def test_velocity_law(test_args: Args) -> None:
    result = solve_for_vector(law.vector_law, law.velocity).subs({
        law.momentum: test_args.p,
        law.total_energy: test_args.e,
    })
    result = QuantityCoordinateVector.from_expr(result)

    assert_equal_vectors(result, test_args.v)


def test_energy_law(test_args: Args) -> None:
    result = law.energy_law.rhs.subs({
        law.momentum: test_args.p,
        law.velocity: test_args.v,
    })

    assert_equal(result, test_args.e)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_momentum(eb, test_args.v)
    with raises(TypeError):
        law.calculate_momentum(100, test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityCoordinateVector([
        Quantity(1 * units.coulomb),
        Quantity(1 * units.coulomb),
        Quantity(1 * units.coulomb),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_momentum(test_args.e, vb_vector)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(ValueError):
        law.calculate_momentum(test_args.e, vb_scalar)

    with raises(TypeError):
        law.calculate_momentum(test_args.e, 100)
    with raises(TypeError):
        law.calculate_momentum(test_args.e, [100])
