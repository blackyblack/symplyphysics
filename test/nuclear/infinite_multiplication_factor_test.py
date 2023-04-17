from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import Probability
from symplyphysics.laws.nuclear import infinite_multiplication_factor

@fixture
def test_args():
    neutron_reproduction = 2.02
    fast_fission = 1.03
    resonance_escape_probability = Probability(0.75)
    thermal_utilisation = Probability(0.7)
    Args = namedtuple("Args", ["n", "e", "p", "f"])
    return Args(n=neutron_reproduction, e=fast_fission, p=resonance_escape_probability, f=thermal_utilisation)

def test_basic_multiplication_factor(test_args):
    result = infinite_multiplication_factor.calculate_multiplication_factor(test_args.n, test_args.e, test_args.p, test_args.f)
    assert result == approx(1.092315, 0.01)

def test_bad_resonance_escape(test_args):
    with raises(AttributeError):
        infinite_multiplication_factor.calculate_multiplication_factor(test_args.n, test_args.e, 100, test_args.f)

def test_bad_thermal_utilisation(test_args):
    with raises(AttributeError):
        infinite_multiplication_factor.calculate_multiplication_factor(test_args.n, test_args.e, test_args.p, 100)
