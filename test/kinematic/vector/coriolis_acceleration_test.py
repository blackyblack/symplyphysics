from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    QuantityVector,
    assert_equal_vectors,
)
from symplyphysics.laws.kinematic.vector import coriolis_acceleration as law

Args = namedtuple("Args", "w v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w_unit = units.radian / units.second
    w = QuantityVector([1 * w_unit, -1 * w_unit, 0])

    v_unit = units.meter / units.second
    v = QuantityVector([1 * v_unit, 1 * v_unit, 1 * v_unit])
    return Args(w=w, v=v)


def test_law(test_args: Args) -> None:
    result = law.calculate_coriolis_acceleration(test_args.w, test_args.v)

    a_unit = units.meter / units.second**2
    assert_equal_vectors(result, QuantityVector([2 * a_unit, 2 * a_unit, -4 * a_unit]))


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityVector([1 * units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_coriolis_acceleration(wb_vector, test_args.v)

    wb_scalar = 3 * units.radian / units.second
    with raises(AttributeError):
        law.calculate_coriolis_acceleration(wb_scalar, test_args.v)

    with raises(TypeError):
        law.calculate_coriolis_acceleration(100, test_args.v)
    with raises(TypeError):
        law.calculate_coriolis_acceleration([100], test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityVector([1 * units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_coriolis_acceleration(test_args.w, vb_vector)

    vb_scalar = 3 * units.meter / units.second
    with raises(AttributeError):
        law.calculate_coriolis_acceleration(test_args.w, vb_scalar)

    with raises(TypeError):
        law.calculate_coriolis_acceleration(test_args.w, 100)
    with raises(TypeError):
        law.calculate_coriolis_acceleration(test_args.w, [100])
