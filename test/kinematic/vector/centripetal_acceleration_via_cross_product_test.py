from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    QuantityVector,
)
from symplyphysics.core.approx import assert_equal_vectors
from symplyphysics.laws.kinematic.vector import centripetal_acceleration_via_cross_product as law

Args = namedtuple("Args", "w r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w_unit = units.radian / units.second
    w = QuantityVector([0, 1 * w_unit, -1 * w_unit])
    r = QuantityVector([1 * units.meter, 2 * units.meter, -1 * units.meter])
    return Args(w=w, r=r)


def test_law(test_args: Args) -> None:
    result = law.calculate_centripetal_acceleration(test_args.w, test_args.r)
    a_unit = units.meter / units.second**2
    assert_equal_vectors(
        result,
        QuantityVector([-2 * a_unit, -1 * a_unit, -1 * a_unit]),
    )


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityVector([1 * units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_centripetal_acceleration(wb_vector, test_args.r)

    wb_scalar = Quantity(1 * units.radian / units.second)
    with raises(AttributeError):
        law.calculate_centripetal_acceleration(wb_scalar, test_args.r)

    with raises(TypeError):
        law.calculate_centripetal_acceleration(100, test_args.r)
    with raises(TypeError):
        law.calculate_centripetal_acceleration([100], test_args.r)


def test_bad_position(test_args: Args) -> None:
    rb_vector = QuantityVector([1 * units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_centripetal_acceleration(test_args.w, rb_vector)
    
    rb_scalar = Quantity(1 * units.meter)
    with raises(AttributeError):
        law.calculate_centripetal_acceleration(test_args.w, rb_scalar)
    
    with raises(TypeError):
        law.calculate_centripetal_acceleration(test_args.w, 100)
    with raises(TypeError):
        law.calculate_centripetal_acceleration(test_args.w, [100])
