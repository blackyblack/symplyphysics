from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import (
    angular_wavenumber_is_inverse_wavelength as wavenumber_def
)

# Description
## A wave has a wavelength of 3 m, then its angular wavenumber is 2.09 rad/m.

Args = namedtuple("Args", "lambda_")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    lambda_ = Quantity(3.0 * units.meter)
    return Args(lambda_=lambda_)


def test_definition(test_args: Args) -> None:
    result = wavenumber_def.calculate_wavenumber(test_args.lambda_)
    assert_equal(result, 2.09 * units.radian / units.meter, tolerance=3e-3)


def test_bad_wavelength(test_args: Args) -> None:
    lambda_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wavenumber_def.calculate_wavenumber(lambda_bad)
    with raises(TypeError):
        wavenumber_def.calculate_wavenumber(100)
