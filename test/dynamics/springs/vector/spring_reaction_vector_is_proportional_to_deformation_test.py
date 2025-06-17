from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.dynamics.springs.vector import spring_reaction_is_proportional_to_deformation as spring_law

from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector, CARTESIAN

Args = namedtuple("Args", ["k", "d", "f"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k = Quantity(0.1 * units.newton / units.meter)
    d_x = Quantity(3 * units.meter)
    d_y = Quantity(1 * units.meter)
    d = QuantityCoordinateVector([d_x, d_y, 0], CARTESIAN)
    f_x = Quantity(-0.3 * units.newton)
    f_y = Quantity(-0.1 * units.newton)
    f = QuantityCoordinateVector([f_x, f_y, 0], CARTESIAN)
    return Args(k=k, d=d, f=f)


def test_basic_force(test_args: Args) -> None:
    result = spring_law.calculate_force(test_args.k, test_args.d)
    expected = QuantityCoordinateVector([-0.3 * units.newton, -0.1 * units.newton, 0], CARTESIAN)
    assert_equal_vectors(result, expected)


def test_basic_deformation(test_args: Args) -> None:
    result = spring_law.calculate_deformation(test_args.k, test_args.f)
    expected = QuantityCoordinateVector([3 * units.meter, 1 * units.meter, 0], CARTESIAN)
    assert_equal_vectors(result, expected)


def test_bad_elastic_coefficient(test_args: Args) -> None:
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        spring_law.calculate_force(eb, test_args.d)
    with raises(TypeError):
        spring_law.calculate_force(100, test_args.d)
    with raises(errors.UnitsError):
        spring_law.calculate_deformation(eb, test_args.f)
    with raises(TypeError):
        spring_law.calculate_deformation(100, test_args.f)


def test_bad_deformation(test_args: Args) -> None:
    db = Quantity(1 * units.coulomb)
    vb = QuantityCoordinateVector([db, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        spring_law.calculate_force(test_args.k, vb)
    with raises(TypeError):
        spring_law.calculate_force(test_args.k, 100)


def test_bad_force(test_args: Args) -> None:
    db = Quantity(1 * units.coulomb)
    vb = QuantityCoordinateVector([db, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        spring_law.calculate_deformation(test_args.k, vb)
    with raises(TypeError):
        spring_law.calculate_deformation(test_args.k, 100)
