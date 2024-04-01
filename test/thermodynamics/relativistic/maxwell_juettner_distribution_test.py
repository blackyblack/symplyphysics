from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.relativistic import (
    maxwell_juettner_distribution as distribution_law,)

# Description
## The value of the Maxwell-Juettner distribution for dimensionless temperature theta = 100 and Lorentz
## factor gamma = 100 is 1.84e-3

Args = namedtuple("Args", "g t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    g = 100
    t = 100
    return Args(g=g, t=t)


def test_law(test_args: Args) -> None:
    result = distribution_law.calculate_distribution_function(test_args.g, test_args.t)
    assert_equal(result, 1.84e-3)


def test_bad_input(test_args: Args) -> None:
    b = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        distribution_law.calculate_distribution_function(b, test_args.t)
    with raises(errors.UnitsError):
        distribution_law.calculate_distribution_function(test_args.g, b)
