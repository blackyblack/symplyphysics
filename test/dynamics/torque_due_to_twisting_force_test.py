from collections import namedtuple
from sympy import pi
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import torque_due_to_twisting_force as torque_def

# Description
## A turning force of magnitude 3 N is applied at a distance of 5 cm away from the rotation axis.
## The angle between the force and the position vector of the point of application of the force
## is pi/6. The resulting torque due to the force should amount to 0.075 N*m.


@fixture(name="test_args")
def test_args_fixture():
    F = Quantity(3.0 * units.newton)
    r = Quantity(5.0 * units.centimeter)
    phi_float = pi/6
    phi_quantity = Quantity(phi_float * units.radian)
    Args = namedtuple("Args", "F r phi_float phi_quantity")
    return Args(F=F, r=r, phi_float=phi_float, phi_quantity=phi_quantity)


def test_basic_law_float_angle(test_args):
    result = torque_def.calculate_torque(test_args.F, test_args.r, test_args.phi_float)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force * units.length)
    result_torque = convert_to(result, units.newton * units.meter).evalf(3)
    assert result_torque == approx(0.075, 1e-4)


def test_basic_law_quantity_angle(test_args):
    result = torque_def.calculate_torque(test_args.F, test_args.r, test_args.phi_quantity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force * units.length)
    result_torque = convert_to(result, units.newton * units.meter).evalf(3)
    assert result_torque == approx(0.075, 1e-4)


def test_bad_force(test_args):
    Fb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        torque_def.calculate_torque(Fb, test_args.r, test_args.phi_float)
    with raises(TypeError):
        torque_def.calculate_torque(100, test_args.r, test_args.phi_float)


def test_bad_distance(test_args):
    rb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        torque_def.calculate_torque(test_args.F, rb, test_args.phi_float)
    with raises(TypeError):
        torque_def.calculate_torque(test_args.F, 100, test_args.phi_float)


def test_bad_angle(test_args):
    phib = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        torque_def.calculate_torque(test_args.F, test_args.r, phib)
