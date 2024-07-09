from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal_vectors,
    errors,
    units,
    Quantity,
    QuantityVector,
)
from symplyphysics.laws.dynamics.vector import torque_vector_of_twisting_force as torque_def

# Description
## A force (1, -2, -1) N is applied at a point (0, 2, -3) m. The torque applied
## should amount to (-8, -3, -2) N*m.

Args = namedtuple("Args", "F r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    F = QuantityVector([
        Quantity(1.0 * units.newton),
        Quantity(-2.0 * units.newton),
        Quantity(-1.0 * units.newton)
    ])
    r = QuantityVector([
        Quantity(0.0 * units.meter),
        Quantity(2.0 * units.meter),
        Quantity(-3.0 * units.meter),
    ])
    return Args(F=F, r=r)


def test_basic_law(test_args: Args) -> None:
    result = torque_def.calculate_torque(test_args.r, test_args.F)
    unit = units.newton * units.meter
    assert_equal_vectors(
        result,
        QuantityVector([-8 * unit, -3 * unit, -2 * unit]),
    )


def test_bad_force(test_args: Args) -> None:
    # Fb is vector, but not in the correct dimension
    with raises(errors.UnitsError):
        Fb = QuantityVector([
            Quantity(1.0 * units.second),
            Quantity(-2.0 * units.second),
            Quantity(-1.0 * units.second)
        ])
        torque_def.calculate_torque(test_args.r, Fb)
    with raises(TypeError):
        torque_def.calculate_torque(test_args.r, 100)


def test_bad_distance(test_args: Args) -> None:
    with raises(errors.UnitsError):
        rb = QuantityVector([
            Quantity(1.0 * units.second),
            Quantity(-2.0 * units.second),
            Quantity(-1.0 * units.second)
        ])
        torque_def.calculate_torque(rb, test_args.F)
    with raises(TypeError):
        torque_def.calculate_torque(100, test_args.F)
