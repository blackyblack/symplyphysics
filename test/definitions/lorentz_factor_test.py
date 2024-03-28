from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import lorentz_factor

# Description
## The speed of an object is 2e5 km/s. Its Lorentz factor is 1.32

Args = namedtuple("Args", "v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v = Quantity(2e5 * units.kilometer / units.second)
    return Args(v=v)


def test_law(test_args: Args) -> None:
    result = lorentz_factor.calculate_lorentz_factor(test_args.v)
    assert_equal(result, 1.34, tolerance=2e-3)


def test_bad_velocity() -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        lorentz_factor.calculate_lorentz_factor(vb)
    with raises(TypeError):
        lorentz_factor.calculate_lorentz_factor(100)
