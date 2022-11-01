from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors, S, pi
)

from sympy.physics.units.definitions.dimension_definitions import angle

from symplyphysics.laws.kinematic import vector_projection_to_axis as projection_law

#We are having force vector of 3 Newtons with angle 60 degrees to axis. Projection of this vector to this axis should be 2 times less than the vector.
@fixture
def test_args():
    input_vector = units.Quantity('input_vector')
    SI.set_quantity_dimension(input_vector, units.force)
    SI.set_quantity_scale_factor(input_vector, 3 * units.newton)

    angle_between_input_vector_and_axis = units.Quantity('angle_between_input_vector_and_axis')
    SI.set_quantity_dimension(angle_between_input_vector_and_axis, angle)
    SI.set_quantity_scale_factor(angle_between_input_vector_and_axis, 60 * units.degree)

    Args = namedtuple('Args', ['input_vector', 'angle_between_input_vector_and_axis'])
    return Args(input_vector = input_vector, angle_between_input_vector_and_axis = angle_between_input_vector_and_axis)

def test_basic_projection(test_args):
    result = projection_law.get_projection(test_args.input_vector, test_args.angle_between_input_vector_and_axis)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_vector = convert_to(result, units.force).subs(units.newton, 1).evalf(2)
    assert result_vector == approx(1.5, 0.01)

