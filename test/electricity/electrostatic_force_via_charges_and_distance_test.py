from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity, prefixes)
from symplyphysics.laws.electricity import electrostatic_force_via_charges_and_distance as couloumb_law_scalar

# Description. With help of online calculator at https://www.wolframalpha.com/widgets/gallery/view.jsp?id=5227fa9a19dce7ba113f50a405dcaf09
## The force between the charges q1 = 10 [nC] and q2 = -50 [nC], which are separated by a distance of r = 0.02 [m], is equivalent to -0.011234 [N].

Args = namedtuple("Args", ["q1", "q2", "r"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q1 = Quantity(10 * prefixes.nano * units.coulomb)
    q2 = Quantity(-50 * prefixes.nano * units.coulomb)
    r = Quantity(0.02 * units.meter)
    return Args(q1=q1, q2=q2, r=r)


def test_basic_force(test_args: Args) -> None:
    result = couloumb_law_scalar.calculate_force(test_args.q1, test_args.q2, test_args.r)
    assert_equal(result, -0.0112344 * units.newton)


def test_bad_charge(test_args: Args) -> None:
    q0 = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        couloumb_law_scalar.calculate_force(q0, test_args.q2, test_args.r)
    with raises(errors.UnitsError):
        couloumb_law_scalar.calculate_force(test_args.q1, q0, test_args.r)
    with raises(TypeError):
        couloumb_law_scalar.calculate_force(1, test_args.q2, test_args.r)
    with raises(TypeError):
        couloumb_law_scalar.calculate_force(test_args.q1, 1, test_args.r)


def test_bad_distance(test_args: Args) -> None:
    r0 = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        couloumb_law_scalar.calculate_force(test_args.q1, test_args.q2, r0)
    with raises(TypeError):
        couloumb_law_scalar.calculate_force(test_args.q1, test_args.q2, 1)
