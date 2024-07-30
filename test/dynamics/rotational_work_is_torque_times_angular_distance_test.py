from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import rotational_work_is_torque_times_angular_distance as work_law

# Description
## A torque of magnitude 3.0 N*m is applied to a rigid body. The work done on the body amounts to
## 9 J when the body rotates 3 rad about the rotational axis.

Args = namedtuple("Args", "tau theta")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    tau = Quantity(3.0 * units.newton * units.meter)
    theta = Quantity(3.0 * units.radian)
    return Args(tau=tau, theta=theta)


def test_basic_law(test_args: Args) -> None:
    result = work_law.calculate_work(test_args.tau, test_args.theta)
    assert_equal(result, 9 * units.joule)


def test_bad_torque(test_args: Args) -> None:
    tau_bad = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        work_law.calculate_work(tau_bad, test_args.theta)
    with raises(TypeError):
        work_law.calculate_work(100, test_args.theta)


def test_bad_angular_displacement(test_args: Args) -> None:
    theta_bad = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.tau, theta_bad)
