from collections import namedtuple
from pytest import fixture, raises
from sympy import abc, cos, sin
from symplyphysics import (
    errors,
    units,
    Vector,
    Quantity,
    QuantityVector,
    assert_equal_vectors,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic.vector import acceleration_due_to_non_uniform_rotation as law


Args = namedtuple("Args", "dw dt r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dw = QuantityVector([0, 0.001 * units.radian / units.second, -0.002 * units.radian / units.second])
    dt = Quantity(0.001 * units.second)
    r = QuantityVector([1 * units.meter, -3 * units.meter, 0.5 * units.meter])
    return Args(dw=dw, dt=dt, r=r)


def test_law(test_args: Args) -> None:
    result = law.calculate_non_uniform_rotation_acceleration(test_args.dw, test_args.dt, test_args.r)

    a_unit = units.meter / units.second**2
    assert_equal_vectors(
        result,
        QuantityVector([-5.5 * a_unit, -2.0 * a_unit, -1.0 * a_unit]),
    )


def test_function_law() -> None:
    t = abc.t
    w = Vector([cos(t), 3 * sin(t), 1])
    r = Vector([1, -3, 2])
    a_result = law.non_uniform_rotation_acceleration_law(w, t, r)
    a_correct = Vector([6 * cos(t), 2 * sin(t), 3 * (sin(t) - cos(t))])
    for _result, _correct in zip(a_result.components, a_correct.components):
        expr_equals(_result, _correct)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_non_uniform_rotation_acceleration(wb_vector, test_args.dt, test_args.r)

    wb_scalar = units.radian / units.second
    with raises(AttributeError):
        law.calculate_non_uniform_rotation_acceleration(wb_scalar, test_args.dt, test_args.r)

    with raises(TypeError):
        law.calculate_non_uniform_rotation_acceleration(100, test_args.dt, test_args.r)
    with raises(TypeError):
        law.calculate_non_uniform_rotation_acceleration([100], test_args.dt, test_args.r)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, tb, test_args.r)
    with raises(TypeError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, 100, test_args.r)


def test_bad_position_vector(test_args: Args) -> None:
    rb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, test_args.dt, rb_vector)

    rb_scalar = units.meter
    with raises(AttributeError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, test_args.dt, rb_scalar)

    with raises(AttributeError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, test_args.dt, 100)
    with raises(AttributeError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, test_args.dt, [100])
