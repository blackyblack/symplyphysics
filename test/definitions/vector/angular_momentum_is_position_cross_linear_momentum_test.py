from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    CoordinateSystem,
    assert_equal,
    units,
    Quantity,
    QuantityVector,
    SI,
    errors,
)
from symplyphysics.definitions.vector import (
    angular_momentum_is_position_cross_linear_momentum as angular_momentum_def,)

# Description
## A particle is located at a point with position vector of (1, 2, -1) m and possesses
## a linear momentum (0, -1, 2) kg*m/s. Its angular momentum amounts to (3, -2, -1) kg*m**2/s.

Args = namedtuple("Args", "r p")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    r = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(2.0 * units.meter),
        Quantity(-1.0 * units.meter),
    ])
    p = QuantityVector([
        Quantity(0.0 * units.kilogram * units.meter / units.second),
        Quantity(-1.0 * units.kilogram * units.meter / units.second),
        Quantity(2.0 * units.kilogram * units.meter / units.second),
    ])
    return Args(r=r, p=p)


def test_definition(test_args: Args) -> None:
    result = angular_momentum_def.calculate_angular_momentum(test_args.r, test_args.p)
    assert len(result.components) == 3
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
        units.length * units.momentum)
    for result_component, correct_value in zip(result.components, [3.0, -2.0, -1.0]):
        assert_equal(result_component,
            correct_value * units.kilogram * units.meter**2 / units.second)


def test_bad_position_vector(test_args: Args) -> None:
    r_bad_vector = QuantityVector([
        Quantity(0.0 * units.meter / units.second),
        Quantity(-1.0 * units.meter / units.second),
        Quantity(2.0 * units.meter / units.second),
    ])
    with raises(errors.UnitsError):
        angular_momentum_def.calculate_angular_momentum(r_bad_vector, test_args.p)

    r_scalar = Quantity(1.0 * units.meter)
    with raises(AttributeError):
        angular_momentum_def.calculate_angular_momentum(r_scalar, test_args.p)

    with raises(TypeError):
        angular_momentum_def.calculate_angular_momentum(100, test_args.p)
    with raises(TypeError):
        angular_momentum_def.calculate_angular_momentum([100], test_args.p)

    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    r_c1 = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.radian),
        Quantity(-1.0 * units.meter),
    ], C1)
    with raises(ValueError):
        angular_momentum_def.calculate_angular_momentum(r_c1, test_args.p)
    C2 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    r_c2 = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.radian),
        Quantity(1.0 * units.radian),
    ], C2)
    with raises(ValueError):
        angular_momentum_def.calculate_angular_momentum(r_c2, test_args.p)


def test_bad_linear_momentum(test_args: Args) -> None:
    p_bad_vector = QuantityVector([
        Quantity(0.0 * units.meter / units.second),
        Quantity(-1.0 * units.meter / units.second),
        Quantity(2.0 * units.meter / units.second),
    ])
    with raises(errors.UnitsError):
        angular_momentum_def.calculate_angular_momentum(test_args.r, p_bad_vector)

    p_scalar = Quantity(1.0 * units.kilogram * units.meter / units.second)
    with raises(AttributeError):
        angular_momentum_def.calculate_angular_momentum(test_args.r, p_scalar)

    with raises(TypeError):
        angular_momentum_def.calculate_angular_momentum(test_args.r, 100)
    with raises(TypeError):
        angular_momentum_def.calculate_angular_momentum(test_args.r, [100])

    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    p_c1 = QuantityVector([
        Quantity(1.0 * units.kilogram * units.meter / units.second),
        Quantity(1.0 * units.radian),
        Quantity(-1.0 * units.kilogram * units.meter / units.second),
    ], C1)
    with raises(ValueError):
        angular_momentum_def.calculate_angular_momentum(p_c1, test_args.p)
    C2 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    p_c2 = QuantityVector([
        Quantity(1.0 * units.kilogram * units.meter / units.second),
        Quantity(1.0 * units.radian),
        Quantity(1.0 * units.radian),
    ], C2)
    with raises(ValueError):
        angular_momentum_def.calculate_angular_momentum(p_c2, test_args.p)
