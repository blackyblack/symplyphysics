from collections import namedtuple
from pytest import fixture, raises
from sympy import symbols
from symplyphysics import (
    errors,
    units,
    Quantity,
    QuantityVector,
    scale_vector,
    assert_equal_vectors,
)
from symplyphysics.laws.kinematic.vector import velocity_relative_to_reference_frame as law

Args = namedtuple("Args", "r1 r2 dt v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    r1 = QuantityVector([0, 1, 1], dimension=units.length)
    r2 = QuantityVector([0.1, 0.9, 1.2], dimension=units.length)
    dt = Quantity(0.1 * units.second)
    v = QuantityVector([1, -1, 2], dimension=units.velocity)
    return Args(r1=r1, r2=r2, dt=dt, v=v)


def test_velocity_law(test_args: Args) -> None:
    result = law.calculate_relative_velocity(test_args.r1, test_args.r2, test_args.dt)
    assert_equal_vectors(result, test_args.v)


def test_position_law(test_args: Args) -> None:
    time = symbols("time")

    final_position_vector = law.relative_position_law(
        test_args.r1.to_base_vector(),
        test_args.v.to_base_vector(),
        time,
    )
    final_position_ = QuantityVector.from_base_vector(
        final_position_vector,
        subs={time: test_args.dt},
    )
    assert_equal_vectors(final_position_, test_args.r2)

    initial_position_vector = law.relative_position_law(
        test_args.r2.to_base_vector(),
        scale_vector(-1, test_args.v.to_base_vector()),
        time,
    )
    initial_position_ = QuantityVector.from_base_vector(
        initial_position_vector,
        subs={time: test_args.dt},
    )
    assert_equal_vectors(initial_position_, test_args.r1)


def test_law(test_args: Args) -> None:
    time = symbols("time")

    velocity_ = law.calculate_relative_velocity(test_args.r1, test_args.r2, test_args.dt)

    final_position_vector = law.relative_position_law(
        test_args.r1.to_base_vector(),
        velocity_.to_base_vector(),
        time,
    )
    final_position_ = QuantityVector.from_base_vector(
        final_position_vector,
        subs={time: test_args.dt},
    )
    assert_equal_vectors(final_position_, test_args.r2)

    initial_position_vector = law.relative_position_law(
        test_args.r2.to_base_vector(),
        scale_vector(-1, velocity_.to_base_vector()),
        time,
    )
    initial_position_ = QuantityVector.from_base_vector(
        initial_position_vector,
        subs={time: test_args.dt},
    )
    assert_equal_vectors(initial_position_, test_args.r1)


def test_bad_position(test_args: Args) -> None:
    rb_vector = QuantityVector([Quantity(1 * units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_relative_velocity(rb_vector, test_args.r2, test_args.dt)
    with raises(errors.UnitsError):
        law.calculate_relative_velocity(test_args.r1, rb_vector, test_args.dt)
    
    rb_scalar = Quantity(1 * units.meter)
    with raises(AttributeError):
        law.calculate_relative_velocity(rb_scalar, test_args.r2, test_args.dt)
    with raises(AttributeError):
        law.calculate_relative_velocity(test_args.r1, rb_scalar, test_args.dt)
    
    with raises(TypeError):
        law.calculate_relative_velocity(100, test_args.r2, test_args.dt)
    with raises(TypeError):
        law.calculate_relative_velocity(test_args.r1, 100, test_args.dt)    
    with raises(TypeError):
        law.calculate_relative_velocity([100], test_args.r2, test_args.dt)
    with raises(TypeError):
        law.calculate_relative_velocity(test_args.r1, [100], test_args.dt)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_relative_velocity(test_args.r1, test_args.r2, tb)
    with raises(TypeError):
        law.calculate_relative_velocity(test_args.r1, test_args.r2, 100)
