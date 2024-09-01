from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    QuantityVector,
    assert_equal_vectors,
)
from symplyphysics.laws.kinematics.vector import acceleration_of_transfer_between_relative_frames as law

Args = namedtuple("Args", "a_0 a_centr a_rot a_tr")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a_unit = units.meter / units.second**2

    a_0 = QuantityVector([a_unit, 0, 0])
    a_centr = QuantityVector([0, a_unit, 2 * a_unit])
    a_rot = QuantityVector([-1 * a_unit, 0, 3 * a_unit])
    a_tr = QuantityVector([0, a_unit, 5 * a_unit])

    return Args(a_0=a_0, a_centr=a_centr, a_rot=a_rot, a_tr=a_tr)


def test_transfer_law(test_args: Args) -> None:
    result = law.calculate_transfer_acceleration(test_args.a_0, test_args.a_centr, test_args.a_rot)
    assert_equal_vectors(result, test_args.a_tr)


def test_moving_frame_law(test_args: Args) -> None:
    result_vector = law.moving_frame_acceleration_law(
        test_args.a_tr.to_base_vector(),
        test_args.a_centr.to_base_vector(),
        test_args.a_rot.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(result_vector)
    assert_equal_vectors(result, test_args.a_0)


def test_centripetal_law(test_args: Args) -> None:
    result_vector = law.centripetal_acceleration_law(
        test_args.a_tr.to_base_vector(),
        test_args.a_0.to_base_vector(),
        test_args.a_rot.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(result_vector)
    assert_equal_vectors(result, test_args.a_centr)


def test_non_uniform_rotation_law(test_args: Args) -> None:
    result_vector = law.rotation_acceleration_law(
        test_args.a_tr.to_base_vector(),
        test_args.a_0.to_base_vector(),
        test_args.a_centr.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(result_vector)
    assert_equal_vectors(result, test_args.a_rot)


def test_bad_moving_frame_acceleration(test_args: Args) -> None:
    ab_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_transfer_acceleration(ab_vector, test_args.a_centr, test_args.a_rot)

    ab_scalar = units.meter / units.second**2
    with raises(AttributeError):
        law.calculate_transfer_acceleration(ab_scalar, test_args.a_centr, test_args.a_rot)

    with raises(TypeError):
        law.calculate_transfer_acceleration(100, test_args.a_centr, test_args.a_rot)
    with raises(TypeError):
        law.calculate_transfer_acceleration([100], test_args.a_centr, test_args.a_rot)


def test_bad_centripetal_acceleration(test_args: Args) -> None:
    ab_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_transfer_acceleration(test_args.a_0, ab_vector, test_args.a_rot)

    ab_scalar = units.meter / units.second**2
    with raises(AttributeError):
        law.calculate_transfer_acceleration(test_args.a_0, ab_scalar, test_args.a_rot)

    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a_0, 100, test_args.a_rot)
    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a_0, [100], test_args.a_rot)


def test_bad_rotation_acceleration(test_args: Args) -> None:
    ab_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_transfer_acceleration(test_args.a_0, test_args.a_centr, ab_vector)

    ab_scalar = units.meter / units.second**2
    with raises(AttributeError):
        law.calculate_transfer_acceleration(test_args.a_0, test_args.a_centr, ab_scalar)

    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a_0, test_args.a_centr, 100)
    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a_0, test_args.a_centr, [100])
