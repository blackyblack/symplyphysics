from collections import namedtuple
from sympy import pi
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import torque_via_force_and_radial_distance as torque_def

# Description
## A turning force of magnitude 3 N is applied at a distance of 5 cm away from the rotation axis.
## The angle between the force and the position vector of the point of application of the force
## is pi/6. The resulting torque due to the force should amount to 0.075 N*m.

Args = namedtuple("Args", "F r phi_float phi_quantity")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    F = Quantity(3.0 * units.newton)
    r = Quantity(5.0 * units.centimeter)
    phi_float = pi / 6
    phi_quantity = Quantity(phi_float * units.radian)
    return Args(F=F, r=r, phi_float=phi_float, phi_quantity=phi_quantity)


def test_basic_law_float_angle(test_args: Args) -> None:
    result = torque_def.calculate_torque(test_args.F, test_args.r, test_args.phi_float)
    assert_equal(result, 0.075 * units.newton * units.meter)


def test_basic_law_quantity_angle(test_args: Args) -> None:
    result = torque_def.calculate_torque(test_args.F, test_args.r, test_args.phi_quantity)
    assert_equal(result, 0.075 * units.newton * units.meter)


def test_bad_force(test_args: Args) -> None:
    Fb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        torque_def.calculate_torque(Fb, test_args.r, test_args.phi_float)
    with raises(TypeError):
        torque_def.calculate_torque(100, test_args.r, test_args.phi_float)


def test_bad_distance(test_args: Args) -> None:
    rb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        torque_def.calculate_torque(test_args.F, rb, test_args.phi_float)
    with raises(TypeError):
        torque_def.calculate_torque(test_args.F, 100, test_args.phi_float)


def test_bad_angle(test_args: Args) -> None:
    phib = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        torque_def.calculate_torque(test_args.F, test_args.r, phib)
