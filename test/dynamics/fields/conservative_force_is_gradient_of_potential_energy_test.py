from collections import namedtuple
from pytest import fixture, raises
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.fields.scalar_field import ScalarField
from symplyphysics.laws.dynamics.fields import conservative_force_is_gradient_of_potential_energy as gradient_law

# Description
## The force associated with a potential of x^2/2 has the vector form of F = -x*e_x, where
## e_x is the unit vector in the direction of the x-axis


@fixture(name="test_args")
def test_args_fixture():
    potential = ScalarField(lambda point: point.x**2 / 2)
    Args = namedtuple("Args", "potential")
    return Args(potential=potential)


def test_basic_law(test_args):
    result_gradient = gradient_law.law(test_args.potential)
    x = result_gradient.coordinate_system.coord_system.base_scalars()[0]
    for i, expr in enumerate([-x, 0, 0]):
        assert expr_equals(result_gradient.components[i], expr)


def test_bad_law():
    with raises(AttributeError):
        gradient_law.law(100)
