from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, Quantity, units, errors
from symplyphysics.core.symbols.probability import Probability
from symplyphysics.laws.thermodynamics import probability_of_ideal_gas_macrostate as macrostate_law

Args = namedtuple("Args", "p1 p2 p3 n1 n2 n3")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p1 = Probability(0.1)
    p2 = Probability(0.4)
    p3 = Probability(0.5)
    n1 = 10
    n2 = 13
    n3 = 5
    return Args(p1=p1, p2=p2, p3=p3, n1=n1, n2=n2, n3=n3)


def test_law(test_args: Args) -> None:
    ps = [test_args.p1, test_args.p2, test_args.p3]
    ns = [test_args.n1, test_args.n2, test_args.n3]
    result = macrostate_law.calculate_macrostate_probability(list(zip(ps, ns)))
    assert_equal(result, 2.36e-6)


def test_bad_argument(test_args: Args) -> None:
    qb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        macrostate_law.calculate_macrostate_probability([(qb, test_args.n1)])
    with raises(errors.UnitsError):
        macrostate_law.calculate_macrostate_probability([(test_args.p1, qb)])
    with raises(TypeError):
        macrostate_law.calculate_macrostate_probability(qb)
    with raises(TypeError):
        macrostate_law.calculate_macrostate_probability(test_args.p1)
    with raises(ValueError):
        macrostate_law.calculate_macrostate_probability(
            list(zip([test_args.p1, test_args.p2], [test_args.n1, test_args.n2]))
        )
