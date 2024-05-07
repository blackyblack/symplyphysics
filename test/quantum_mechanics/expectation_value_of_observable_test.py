from collections import namedtuple
from pytest import fixture
from sympy import Expr, pi, root, exp, Rational
from symplyphysics import assert_equal
from symplyphysics.laws.quantum_mechanics import expectation_value_of_observable as mean_law

Args = namedtuple("Args", "obs psi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    def psi(position: Expr, _time: Expr) -> Expr:
        return root(pi, -4) * exp(-position**2 / 2)
    
    def obs(psi_in: mean_law.WaveFunction) -> mean_law.WaveFunction:
        def psi_out(position: Expr, time: Expr) -> Expr:
            return position**2 * psi_in(position, time)
        return psi_out
    
    return Args(obs=obs, psi=psi)


def test_law(test_args: Args) -> None:
    result = mean_law.calculate_mean_observable_value(test_args.obs, test_args.psi)
    assert_equal(float(result), Rational(1, 2))
