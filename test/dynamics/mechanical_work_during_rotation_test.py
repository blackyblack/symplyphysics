from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import mechanical_work_during_rotation as work_law

# Description
## A torque of magnitude 3.0 N*m is applied to a rigid body. The work done on the body amounts to
## 9 J when the body rotates 3 rad about the rotational axis.


@fixture(name="test_args")
def test_args_fixture():
    tau = Quantity(3.0 * units.newton * units.meter)
    theta = Quantity(3.0 * units.radian)
    Args = namedtuple("Args", "tau theta")
    return Args(tau=tau, theta=theta)


def test_basic_law(test_args):
    result = work_law.calculate_work(test_args.tau, test_args.theta)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_value = convert_to(result, units.joule).evalf(3)
    assert_approx(result_value, 9)


def test_bad_torque(test_args):
    tau_bad = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        work_law.calculate_work(tau_bad, test_args.theta)
    with raises(TypeError):
        work_law.calculate_work(100, test_args.theta)


def test_bad_angular_displacement(test_args):
    theta_bad = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.tau, theta_bad)
