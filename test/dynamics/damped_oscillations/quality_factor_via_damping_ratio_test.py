from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.damped_oscillations import quality_factor_via_damping_ratio as q_factor_law

# Description
## A damped oscillator is described as having damping ratio zeta = 0.01. Then its Q factor
## is Q = 50

Args = namedtuple("Args", "zeta")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    zeta = 0.01
    return Args(zeta=zeta)


def test_law(test_args: Args) -> None:
    result = q_factor_law.calculate_q_factor(test_args.zeta)
    assert_equal(result, 50)


def test_bad_damping_ratio() -> None:
    zeta_bad = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        q_factor_law.calculate_q_factor(zeta_bad)
