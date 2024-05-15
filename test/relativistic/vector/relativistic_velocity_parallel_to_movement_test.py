from collections import namedtuple
from pytest import fixture, raises
from sympy.physics.units import speed_of_light
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    QuantityVector,
)
from symplyphysics.laws.relativistic.vector import relativistic_velocity_parallel_to_movement as law

# Description
## A body is moving parallel to a frame of reference that is moving relative to the laboratory frame of reference.
## The velocity vector of the body in its proper frame of reference is (0.1*c, 0.2*c, 0). The velocity vector of the
## proper frame of reference relative to the laboratory one is (-0.2*c, -0.4*c, 0). Then the velocity of the body relative
## to the laboratory frame of reference is (-0.111*c, -0.222*c, 0).

Args = namedtuple("Args", "u v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    u = QuantityVector([Quantity(0.1 * speed_of_light), Quantity(0.2 * units.speed_of_light), 0])
    v = QuantityVector([Quantity(-0.2 * speed_of_light), Quantity(-0.4 * units.speed_of_light), 0])
    return Args(u=u, v=v)


def test_law(test_args: Args) -> None:
    result = law.calculate_parallel_velocity_component_in_lab_frame(test_args.u, test_args.v)
    assert len(result.components) == 3
    for result_component, correct_value in zip(result.components, [-0.111, -0.222, 0]):
        assert_equal(result_component, correct_value * speed_of_light, tolerance=1e-2)
