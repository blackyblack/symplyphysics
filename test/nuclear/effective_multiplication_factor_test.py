from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics.core.probability import Probability
from symplyphysics.laws.nuclear import effective_multiplication_factor

@fixture
def test_args():
    infinite_multiplication_factor = 1.092
    fast_non_leakage_probability = Probability(0.95)
    thermal_non_leakage_probability = Probability(0.96)
    Args = namedtuple("Args", ["kinf", "Pf", "Pt"])
    return Args(kinf=infinite_multiplication_factor, Pf=fast_non_leakage_probability, Pt=thermal_non_leakage_probability)

# Expected to get k_effective = 1 (critical state)
def test_basic_multiplication_factor(test_args):
    result = effective_multiplication_factor.calculate_multiplication_factor(test_args.kinf, test_args.Pf, test_args.Pt)
    assert result == approx(1, 0.01)

def test_bad_fast_non_leakage(test_args):
    with raises(AttributeError):
        effective_multiplication_factor.calculate_multiplication_factor(test_args.kinf, 100, test_args.Pt)

def test_thermal_non_leakage(test_args):
    with raises(AttributeError):
        effective_multiplication_factor.calculate_multiplication_factor(test_args.kinf, test_args.Pf, 100)
