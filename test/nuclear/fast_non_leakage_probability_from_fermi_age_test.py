from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
)
from symplyphysics.core.symbols.probability import Probability
from symplyphysics.laws.nuclear import fast_non_leakage_probability_from_fermi_age as non_leakage_factor


@fixture
def test_args():
    # sphere with radius = 1 meter
    geometric_buckling = Quantity(9.869 / units.meter**2)
    # water Fermi age = 27 cm^2
    neutron_fermi_age = Quantity(27 * units.centimeter**2)
    Args = namedtuple("Args", ["Bg", "th"])
    return Args(Bg=geometric_buckling, th=neutron_fermi_age)


def test_basic_non_leakage_factor(test_args):
    result = non_leakage_factor.calculate_probability(test_args.Bg, test_args.th)
    assert isinstance(result, Probability)
    assert result.value == approx(0.9737, 0.01)


def test_bad_buckling(test_args):
    Bgb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        non_leakage_factor.calculate_probability(Bgb, test_args.th)
    with raises(TypeError):
        non_leakage_factor.calculate_probability(100, test_args.th)


def test_bad_fermi_age(test_args):
    thb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        non_leakage_factor.calculate_probability(test_args.Bg, thb)
    with raises(TypeError):
        non_leakage_factor.calculate_probability(test_args.Bg, 100)
