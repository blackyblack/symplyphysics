from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    QuantityVector,
    assert_equal_vectors,
)
from symplyphysics.laws.kinematics.vector import absolute_velocity_of_arbitrary_motion as law

Args = namedtuple("Args", "vabs vrel vtr")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    vabs = QuantityVector([-1, 3, -4], dimension=units.velocity)
    vrel = QuantityVector([0, 1, -1], dimension=units.velocity)
    vtr = QuantityVector([-1, 2, -3], dimension=units.velocity)
    return Args(vabs=vabs, vrel=vrel, vtr=vtr)


def test_absolute_law(test_args: Args) -> None:
    result = law.calculate_absolute_velocity(test_args.vrel, test_args.vtr)
    assert_equal_vectors(result, test_args.vabs)


def test_relative_law(test_args: Args) -> None:
    result_vector = law.relative_velocity_law(
        test_args.vabs.to_base_vector(),
        test_args.vtr.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(result_vector)
    assert_equal_vectors(result, test_args.vrel)


def test_transfer_law(test_args: Args) -> None:
    result_vector = law.transfer_velocity_law(
        test_args.vabs.to_base_vector(),
        test_args.vrel.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(result_vector)
    assert_equal_vectors(result, test_args.vtr)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityVector([Quantity(1 * units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_absolute_velocity(vb_vector, test_args.vtr)
    with raises(errors.UnitsError):
        law.calculate_absolute_velocity(test_args.vrel, vb_vector)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(AttributeError):
        law.calculate_absolute_velocity(vb_scalar, test_args.vtr)
    with raises(AttributeError):
        law.calculate_absolute_velocity(test_args.vrel, vb_scalar)

    with raises(TypeError):
        law.calculate_absolute_velocity(100, test_args.vtr)
    with raises(TypeError):
        law.calculate_absolute_velocity([100], test_args.vtr)
    with raises(TypeError):
        law.calculate_absolute_velocity(test_args.vrel, 100)
    with raises(TypeError):
        law.calculate_absolute_velocity(test_args.vrel, [100])
