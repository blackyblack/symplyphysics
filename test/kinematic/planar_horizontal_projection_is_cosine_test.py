from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI
)

from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

from symplyphysics.laws.kinematic import planar_horizontal_projection_is_cosine as projection_law

#We are having force vector of 3 Newtons with angle 60 degrees to horizontal axis. Projection of this vector to horizontal axis should be 2 times less than the vector.
@fixture
def test_args():
    force_vector_amplitude = units.Quantity('force_vector_amplitude')
    SI.set_quantity_dimension(force_vector_amplitude, units.force)
    SI.set_quantity_scale_factor(force_vector_amplitude, 3 * units.newton)
    angle_between_vector_and_horizontal_axis = units.Quantity('angle_between_vector_and_horizontal_axis')
    SI.set_quantity_dimension(angle_between_vector_and_horizontal_axis, angle_type)
    SI.set_quantity_scale_factor(angle_between_vector_and_horizontal_axis, 60 * units.degree)

    Args = namedtuple('Args', ['force_vector_amplitude', 'angle_between_vector_and_horizontal_axis'])
    return Args(force_vector_amplitude = force_vector_amplitude, angle_between_vector_and_horizontal_axis = angle_between_vector_and_horizontal_axis)

def test_basic_projection(test_args):
    result = projection_law.calculate_projection(test_args.force_vector_amplitude, test_args.angle_between_vector_and_horizontal_axis)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_vector = convert_to(result, units.newton).subs(units.newton, 1).evalf(2)
    assert result_vector == approx(1.5, 0.01)

