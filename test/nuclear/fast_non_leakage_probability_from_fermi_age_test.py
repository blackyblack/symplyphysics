from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.core.symbols.probability import Probability
from symplyphysics.laws.nuclear import fast_non_leakage_probability_from_fermi_age as non_leakage_factor

Args = namedtuple("Args", ["Bg", "th"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # sphere with radius = 1 meter
    geometric_buckling = Quantity(9.869 / units.meter**2)
    # water Fermi age = 27 cm^2
    neutron_fermi_age = Quantity(27 * units.centimeter**2)
    return Args(Bg=geometric_buckling, th=neutron_fermi_age)


def test_basic_non_leakage_factor(test_args: Args) -> None:
    result = non_leakage_factor.calculate_probability(test_args.Bg, test_args.th)
    assert isinstance(result, Probability)
    assert_equal(result, 0.9737)


def test_bad_buckling(test_args: Args) -> None:
    Bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        non_leakage_factor.calculate_probability(Bb, test_args.th)
    with raises(TypeError):
        non_leakage_factor.calculate_probability(100, test_args.th)


def test_bad_fermi_age(test_args: Args) -> None:
    thb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        non_leakage_factor.calculate_probability(test_args.Bg, thb)
    with raises(TypeError):
        non_leakage_factor.calculate_probability(test_args.Bg, 100)
