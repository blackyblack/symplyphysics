from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import rotational_stiffness_is_torque_over_angle as law

Args = namedtuple("Args", "tau phi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    tau = Quantity(4 * units.newton * units.meter)
    phi = Quantity(1 * units.degree)
    return Args(tau=tau, phi=phi)


def test_law(test_args: Args) -> None:
    result = law.calculate_rotational_stiffness(test_args.tau, test_args.phi)
    assert_equal(result, 230 * units.newton * units.meter / units.radian, relative_tolerance=4e-3)


def test_bad_torque(test_args: Args) -> None:
    taub = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_rotational_stiffness(taub, test_args.phi)
    with raises(TypeError):
        law.calculate_rotational_stiffness(100, test_args.phi)


def test_bad_angle(test_args: Args) -> None:
    phib = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_rotational_stiffness(test_args.tau, phib)
