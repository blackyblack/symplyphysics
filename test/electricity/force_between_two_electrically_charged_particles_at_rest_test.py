from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
    prefixes
)
from symplyphysics.laws.electricity import force_between_two_electrically_charged_particles_at_rest as couloumb_law_scalar

# Description. With help of online calculator at https://www.wolframalpha.com/widgets/gallery/view.jsp?id=5227fa9a19dce7ba113f50a405dcaf09
## The force between the charges q1 = 10 [nC] and q2 = -50 [nC], which are separated by a distance of r = 2 [m], is equivalent to -0.011234 [N].

@fixture(name="test_args")
def test_args_fixture():
    q1 = Quantity(10 * prefixes.nano * units.coulomb)
    q2 = Quantity(-50 * prefixes.nano * units.coulomb)
    r = Quantity(0.02 * units.meter)
    Args = namedtuple("Args", ["q1", "q2", "r"])
    return Args(q1=q1, q2=q2, r=r)

def test_basic_force(test_args):
    result = couloumb_law_scalar.calculate_force_btw_two_charged_particles(test_args.q1, test_args.q2, test_args.r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force = convert_to(result, units.newton).evalf(7)
    assert result_force == approx(-0.0112344, 0.00001)

def test_bad_charge(test_args):
    q0 = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        couloumb_law_scalar.calculate_force_btw_two_charged_particles(q0, test_args.q2, test_args.r)
    with raises(errors.UnitsError):
        couloumb_law_scalar.calculate_force_btw_two_charged_particles(test_args.q1, q0, test_args.r)
    with raises(TypeError):
        couloumb_law_scalar.calculate_force_btw_two_charged_particles(1, test_args.q2, test_args.r)
    with raises(TypeError):
        couloumb_law_scalar.calculate_force_btw_two_charged_particles(test_args.q1, 1, test_args.r)

def test_bad_distance(test_args):
    r0 = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        couloumb_law_scalar.calculate_force_btw_two_charged_particles(test_args.q1, test_args.q2, r0)
    with raises(TypeError):
        couloumb_law_scalar.calculate_force_btw_two_charged_particles(test_args.q1, test_args.q2, 1)
