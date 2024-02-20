from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal
from symplyphysics.core.symbols.probability import Probability
from symplyphysics.laws.nuclear import effective_multiplication_factor

Args = namedtuple("Args", ["kinf", "Pf", "Pt"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    infinite_multiplication_factor = 1.092
    fast_non_leakage_probability = Probability(0.95)
    thermal_non_leakage_probability = Probability(0.96)
    return Args(kinf=infinite_multiplication_factor,
        Pf=fast_non_leakage_probability,
        Pt=thermal_non_leakage_probability)


# Expected to get k_effective = 1 (critical state)
def test_basic_multiplication_factor(test_args: Args) -> None:
    result = effective_multiplication_factor.calculate_multiplication_factor(
        test_args.kinf, test_args.Pf, test_args.Pt)
    assert_equal(result, 1, tolerance=0.01)


def test_bad_fast_non_leakage(test_args: Args) -> None:
    with raises(AttributeError):
        effective_multiplication_factor.calculate_multiplication_factor(
            test_args.kinf, Probability(100), test_args.Pt)


def test_thermal_non_leakage(test_args: Args) -> None:
    with raises(AttributeError):
        effective_multiplication_factor.calculate_multiplication_factor(
            test_args.kinf, test_args.Pf, Probability(100))
