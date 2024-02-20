from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
    QuantityVector,
)
from symplyphysics.definitions.vector import (
    damping_force_is_proportional_to_velocity as damping_def
)

# Description
## A damping force acts on a body. The body's velocity is (4, -2, 0) m/s, the damping constant
## of the system is 0.05 kg/s. Then the damping force is (-0.2, 0.1, 0.0) N.

Args = namedtuple("Args", "b v f")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    b = Quantity(0.05 * units.kilogram / units.second)
    v = QuantityVector([
        Quantity(4.0 * units.meter / units.second),
        Quantity(-2.0 * units.meter / units.second),
        Quantity(0.0 * units.meter / units.second),
    ])
    f = QuantityVector([
        Quantity(-0.2 * units.newton),
        Quantity(0.1 * units.newton),
        Quantity(0.0 * units.newton)
    ])
    return Args(b=b, v=v, f=f)


def test_damping_force_definition(test_args: Args) -> None:
    result = damping_def.calculate_damping_force(test_args.b, test_args.v)
    assert len(result.components) == 3
    for result_component, correct_component in zip(result.components, test_args.f.components):
        assert_equal(result_component, correct_component)


def test_velocity_law(test_args: Args) -> None:
    result = damping_def.calculate_velocity(test_args.b, test_args.f)
    assert len(result.components) == 3
    for result_component, correct_component in zip(result.components, test_args.v.components):
        assert_equal(result_component, correct_component)


def test_bad_damping_constant(test_args: Args) -> None:
    bb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        damping_def.calculate_damping_force(bb, test_args.v)
    with raises(TypeError):
        damping_def.calculate_damping_force(100, test_args.v)
    with raises(errors.UnitsError):
        damping_def.calculate_velocity(bb, test_args.f)
    with raises(TypeError):
        damping_def.calculate_velocity(100, test_args.f)


def test_bad_velocity(test_args: Args) -> None:
    v_bad_vector = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ])
    with raises(errors.UnitsError):
        damping_def.calculate_damping_force(test_args.b, v_bad_vector)

    v_scalar = Quantity(1.0 * units.meter / units.second)
    with raises(AttributeError):
        damping_def.calculate_damping_force(test_args.b, v_scalar)

    with raises(TypeError):
        damping_def.calculate_damping_force(test_args.b, 100)
    with raises(TypeError):
        damping_def.calculate_damping_force(test_args.b, [100])


def test_bad_damping_force(test_args: Args) -> None:
    f_bad_vector = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ])
    with raises(errors.UnitsError):
        damping_def.calculate_velocity(test_args.b, f_bad_vector)

    f_scalar = Quantity(1.0 * units.newton)
    with raises(AttributeError):
        damping_def.calculate_velocity(test_args.b, f_scalar)

    with raises(TypeError):
        damping_def.calculate_velocity(test_args.b, 100)
    with raises(TypeError):
        damping_def.calculate_velocity(test_args.b, [100])
