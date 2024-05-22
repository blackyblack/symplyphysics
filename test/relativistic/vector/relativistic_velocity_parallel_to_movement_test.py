from collections import namedtuple
from pytest import fixture, raises
from sympy.physics.units import speed_of_light
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    QuantityVector,
    errors,
)
from symplyphysics.laws.relativistic.vector import relativistic_velocity_parallel_to_movement as law

# Description
## A body is moving parallel to the velocity vector of the proper reference frame relative to the lab frame.
## The velocity vector of the body in its proper frame of reference is (0.1*c, 0.2*c, 0). The velocity vector of the
## proper frame of reference relative to the laboratory one is (-0.2*c, -0.4*c, 0). Then the velocity of the body relative
## to the laboratory frame of reference is (-0.111*c, -0.222*c, 0).

Args = namedtuple("Args", "up ul v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    up = QuantityVector([Quantity(0.1 * speed_of_light), Quantity(0.2 * units.speed_of_light), 0])
    ul = QuantityVector([Quantity(-0.111 * speed_of_light), Quantity(-0.222 * units.speed_of_light), 0])
    v = QuantityVector([Quantity(-0.2 * speed_of_light), Quantity(-0.4 * units.speed_of_light), 0])
    return Args(up=up, ul=ul, v=v)


def test_lab_law(test_args: Args) -> None:
    result = law.calculate_parallel_velocity_component_in_lab_frame(test_args.up, test_args.v)
    assert len(result.components) == 3
    for result_component, correct_component in zip(result.components, test_args.ul.components):
        assert_equal(result_component, correct_component, tolerance=2e-3)


def test_proper_law(test_args: Args) -> None:
    result_vector = law.parallel_velocity_component_in_proper_frame_law(
        test_args.ul.to_base_vector(),
        test_args.v.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(result_vector)
    assert len(result.components) == 3
    for result_component, correct_component in zip(result.components, test_args.up.components):
        assert_equal(result_component, correct_component, tolerance=2e-3)


def test_relative_law(test_args: Args) -> None:
    result_vector = law.proper_frame_velocity_in_lab_frame_law(
        test_args.up.to_base_vector(),
        test_args.ul.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(result_vector)
    assert len(result.components) == 3
    for result_component, correct_component in zip(result.components, test_args.v.components):
        assert_equal(result_component, correct_component, tolerance=2e-3)


def test_non_collinear_velocities() -> None:
    ub = QuantityVector([units.speed_of_light, 0, 0])
    vb_vector = QuantityVector([0, units.speed_of_light, 0])
    with raises(AssertionError):
        law.calculate_parallel_velocity_component_in_lab_frame(ub, vb_vector)


def test_bad_object_velocity(test_args: Args) -> None:
    qb = Quantity(1 * units.coulomb)
    vb_vector = QuantityVector([qb, qb, qb])
    with raises(errors.UnitsError):
        law.calculate_parallel_velocity_component_in_lab_frame(vb_vector, test_args.v)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(AttributeError):
        law.calculate_parallel_velocity_component_in_lab_frame(vb_scalar, test_args.v)

    with raises(TypeError):
        law.calculate_parallel_velocity_component_in_lab_frame(100, test_args.v)
    with raises(TypeError):
        law.calculate_parallel_velocity_component_in_lab_frame([100], test_args.v)


def test_bad_frame_velocity(test_args: Args) -> None:
    qb = Quantity(1 * units.coulomb)
    vb_vector = QuantityVector([qb, qb, qb])
    with raises(errors.UnitsError):
        law.calculate_parallel_velocity_component_in_lab_frame(test_args.up, vb_vector)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(AttributeError):
        law.calculate_parallel_velocity_component_in_lab_frame(test_args.up, vb_scalar)

    with raises(TypeError):
        law.calculate_parallel_velocity_component_in_lab_frame(test_args.up, 100)
    with raises(TypeError):
        law.calculate_parallel_velocity_component_in_lab_frame(test_args.up, [100])
