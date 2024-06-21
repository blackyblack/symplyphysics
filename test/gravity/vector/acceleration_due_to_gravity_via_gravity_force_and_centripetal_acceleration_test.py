from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    QuantityVector,
)
from symplyphysics.core.approx import assert_equal, assert_equal_vectors
from symplyphysics.laws.gravity.vector import acceleration_due_to_gravity_via_gravity_force_and_centripetal_acceleration as law

Args = namedtuple("Args", "f w r m g")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = QuantityVector([-243 * units.newton, -421 * units.newton, -486 * units.newton])
    w = QuantityVector([0, 0, 7.29e-5 * units.radian / units.second])
    r = QuantityVector([2.25e6 * units.meter, 3.9e6 * units.meter, 4.5e6 * units.meter])
    m = Quantity(70 * units.kilogram)

    g_unit = units.meter / units.second**2
    g = QuantityVector([-3.46 * g_unit, -5.99 * g_unit, -6.94 * g_unit])

    return Args(f=f, w=w, r=r, m=m, g=g)


def test_law(test_args: Args) -> None:
    result = law.calculate_acceleraton_due_to_gravity(test_args.f, test_args.w, test_args.r, test_args.m)
    assert_equal_vectors(result, test_args.g)


def test_gravity_law(test_args: Args) -> None:
    vector = law.gravity_force_law(
        test_args.g.to_base_vector(),
        test_args.w.to_base_vector(),
        test_args.r.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(vector, subs={law.mass: test_args.m})
    assert_equal_vectors(result, test_args.f)


def test_mass_law(test_args: Args) -> None:
    result = law.mass_law(
        test_args.g.to_base_vector(),
        test_args.f.to_base_vector(),
        test_args.w.to_base_vector(),
        test_args.r.to_base_vector(),
    )
    assert_equal(result, test_args.m)


def test_bad_force(test_args: Args) -> None:
    fb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_acceleraton_due_to_gravity(fb_vector, test_args.w, test_args.r, test_args.m)

    fb_scalar = Quantity(units.newton)
    with raises(AttributeError):
        law.calculate_acceleraton_due_to_gravity(fb_scalar, test_args.w, test_args.r, test_args.m)

    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity(100, test_args.w, test_args.r, test_args.m)
    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity([100], test_args.w, test_args.r, test_args.m)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, wb_vector, test_args.r, test_args.m)

    wb_scalar = Quantity(units.radian / units.second)
    with raises(AttributeError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, wb_scalar, test_args.r, test_args.m)

    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, 100, test_args.r, test_args.m)
    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, [100], test_args.r, test_args.m)


def test_bad_position(test_args: Args) -> None:
    rb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, test_args.w, rb_vector, test_args.m)

    rb_scalar = Quantity(units.meter)
    with raises(AttributeError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, test_args.w, rb_scalar, test_args.m)

    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, test_args.w, 100, test_args.m)
    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, test_args.w, [100], test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, test_args.w, test_args.r, mb)
    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, test_args.w, test_args.r, 100)
