from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, convert_to, Quantity, SI, QuantityVector)
from symplyphysics.laws.dynamics.vector import torque_vector_of_twisting_force as torque_def

# Description
## A force (1, -2, -1) N is applied at a point (0, 2, -3) m. The torque applied
## should amount to (8, 3, 2) N*m.


@fixture(name="test_args")
def test_args_fixture():
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
    Args = namedtuple("Args", "F r")
    return Args(F=F, r=r)


def test_basic_law(test_args):
    result = torque_def.calculate_torque(test_args.F, test_args.r)
    for (result_component, correct_value) in zip(result.components, [8, 3, 2]):
        assert SI.get_dimension_system().equivalent_dims(result_component.dimension,
            units.force * units.length)
        result_value = convert_to(result_component, units.newton * units.meter).evalf(3)
        assert result_value == approx(correct_value, 1e-3)


def test_bad_force(test_args):
    # Fb is vector, but not in the correct dimension
    with raises(errors.UnitsError):
        Fb = QuantityVector([
            Quantity(1.0 * units.second),
            Quantity(-2.0 * units.second),
            Quantity(-1.0 * units.second)
        ])
        torque_def.calculate_torque(Fb, test_args.r)
    # Fb is in the dimension of a force, but not a vector
    with raises(AttributeError):
        Fb = Quantity(1.0 * units.newton)
        torque_def.calculate_torque(Fb, test_args.r)
    with raises(TypeError):
        torque_def.calculate_torque(100, test_args.r)


def test_bad_distance(test_args):
    with raises(errors.UnitsError):
        rb = QuantityVector([
            Quantity(1.0 * units.second),
            Quantity(-2.0 * units.second),
            Quantity(-1.0 * units.second)
        ])
        torque_def.calculate_torque(test_args.F, rb)
    with raises(AttributeError):
        rb = Quantity(1.0 * units.meter)
        torque_def.calculate_torque(test_args.F, rb)
    with raises(TypeError):
        torque_def.calculate_torque(test_args.F, 100)
