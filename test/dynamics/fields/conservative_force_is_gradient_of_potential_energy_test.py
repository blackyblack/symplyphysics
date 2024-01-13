from collections import namedtuple
from pytest import fixture, raises
from sympy import solve
from sympy.vector import CoordSys3D
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.fields.scalar_field import ScalarField
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.laws.dynamics.fields import conservative_force_is_gradient_of_potential_energy as gradient_law

@fixture(name="test_args")
def test_args_fixture():
    potential = ScalarField(lambda point: point.x**2 / 2)
    Args = namedtuple("Args", "potential")
    return Args(potential=potential)


def test_basic_law(test_args):
    result_field = gradient_law.law(test_args.potential)
    result_vector = result_field.apply_to_basis()

    x = result_field.coordinate_system.coord_system.base_scalars()[0]
    for i, expr in enumerate([-x, 0, 0]):
        assert expr_equals(result_vector.components[i], expr)


def test_bad_law(test_args):
    bad_field = VectorField(lambda point: [point.x, point.y, point.z])
    with raises(Exception):
        gradient_law.law(bad_field)
    with raises(AttributeError):
        gradient_law.law(100)
