from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal_vectors,
    errors,
    units,
    Quantity,
    QuantityVector,
)
from symplyphysics.laws.kinematic.vector import centrifugal_acceleration_via_centripetal_acceleration as law

Args = namedtuple("Args", "acp acf")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a_unit = units.meter / units.second**2

    acp = QuantityVector([2 * a_unit, -1 * a_unit, 3 * a_unit])
    acf = QuantityVector([-2 * a_unit, 1 * a_unit, -3 * a_unit])

    return Args(acp=acp, acf=acf)


def test_centrifugal_law(test_args: Args) -> None:
    result = law.calculate_centrifugal_acceleration(test_args.acp)
    assert_equal_vectors(result, test_args.acf)


def test_centripetal_law(test_args: Args) -> None:
    result_vector = law.centripetal_law(test_args.acf.to_base_vector())
    result = QuantityVector.from_base_vector(result_vector)
    assert_equal_vectors(result, test_args.acp)


def test_bad_acceleration(test_args: Args) -> None:
    ab_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_centrifugal_acceleration(ab_vector)

    ab_scalar = units.meter / units.second**2
    with raises(AttributeError):
        law.calculate_centrifugal_acceleration(ab_scalar)

    with raises(TypeError):
        law.calculate_centrifugal_acceleration(100)
    with raises(TypeError):
        law.calculate_centrifugal_acceleration([100])
