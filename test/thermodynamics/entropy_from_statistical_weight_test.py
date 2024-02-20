from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import entropy_from_statistical_weight as entropy_law

# Description
## Let the statistical weight of the state of the system be equal to 1000000. Then the entropy of the system will be equal
## to 1.907 joule per kelvin.

Args = namedtuple("Args", ["statistical_weight"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    statistical_weight = 1000000
    return Args(statistical_weight=statistical_weight)


def test_basic_entropy(test_args: Args) -> None:
    result = entropy_law.calculate_entropy(test_args.statistical_weight)
    assert_equal(result, 1.907e-22 * units.joule / units.kelvin)


def test_bad_statistical_weight() -> None:
    statistical_weight = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        entropy_law.calculate_entropy(statistical_weight)
