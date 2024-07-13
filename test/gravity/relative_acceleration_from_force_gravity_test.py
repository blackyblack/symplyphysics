from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    QuantityVector,
)
from symplyphysics.laws.gravity import relative_acceleration_from_force as law
from symplyphysics.core.approx import assert_equal_vectors

Args = namedtuple("Args", "g ac f m ar")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a_unit = units.meter / units.second**2
    g = QuantityVector([0, 0, -9.81 * a_unit])
    ac = QuantityVector([0, -0.1 * a_unit, 0.25 * a_unit])
    f = QuantityVector([100 * units.newton, 100 * units.newton, 0])
    m = Quantity(100 * units.kilogram)
    ar = QuantityVector([1.0 * a_unit, 1.1 * a_unit, -10.06 * a_unit])
    return Args(g=g, ac=ac, f=f, m=m, ar=ar)


def test_relative_acceleration_law(test_args: Args) -> None:
    result = law.calculate_acceleration(
        test_args.g,
        test_args.ac,
        test_args.f,
        test_args.m,
    )
    assert_equal_vectors(result, test_args.ar)


def test_gravity_acceleration_law(test_args: Args) -> None:
    result_vector = law.acceleration_due_to_gravity_law(
        test_args.ar.to_base_vector(),
        test_args.ac.to_base_vector(),
        test_args.f.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(
        result_vector,
        subs={law.mass: test_args.m},
    )
    assert_equal_vectors(result, test_args.g, absolute_tolerance=1e-10)


def test_coriolis_acceleration_law(test_args: Args) -> None:
    result_vector = law.coriolis_acceleration_law(
        test_args.ar.to_base_vector(),
        test_args.g.to_base_vector(),
        test_args.f.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(
        result_vector,
        subs={law.mass: test_args.m},
    )
    assert_equal_vectors(result, test_args.ac)


def test_force_law(test_args: Args) -> None:
    result_vector = law.force_law(
        test_args.ar.to_base_vector(),
        test_args.g.to_base_vector(),
        test_args.ac.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(
        result_vector,
        subs={law.mass: test_args.m},
    )
    assert_equal_vectors(result, test_args.f)


def test_mass_law(test_args: Args) -> None:
    result = law.mass_law(
        test_args.ar.to_base_vector(),
        test_args.g.to_base_vector(),
        test_args.ac.to_base_vector(),
        test_args.f.to_base_vector(),
    )
    assert_equal(result, test_args.m)


def test_bad_acceleration(test_args: Args) -> None:
    ab_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_acceleration(ab_vector, test_args.ac, test_args.f, test_args.m)
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.g, ab_vector, test_args.f, test_args.m)

    ab_scalar = units.meter / units.second**2
    with raises(AttributeError):
        law.calculate_acceleration(ab_scalar, test_args.ac, test_args.f, test_args.m)
    with raises(AttributeError):
        law.calculate_acceleration(test_args.g, ab_scalar, test_args.f, test_args.m)

    with raises(TypeError):
        law.calculate_acceleration(100, test_args.ac, test_args.f, test_args.m)
    with raises(TypeError):
        law.calculate_acceleration([100], test_args.ac, test_args.f, test_args.m)
    with raises(TypeError):
        law.calculate_acceleration(test_args.g, 100, test_args.f, test_args.m)
    with raises(TypeError):
        law.calculate_acceleration(test_args.g, [100], test_args.f, test_args.m)


def test_bad_force(test_args: Args) -> None:
    fb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.g, test_args.ac, fb_vector, test_args.m)

    fb_scalar = units.newton
    with raises(AttributeError):
        law.calculate_acceleration(test_args.g, test_args.ac, fb_scalar, test_args.m)

    with raises(TypeError):
        law.calculate_acceleration(test_args.g, test_args.ac, 100, test_args.m)
    with raises(TypeError):
        law.calculate_acceleration(test_args.g, test_args.ac, [100], test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.g, test_args.ac, test_args.f, mb)
    with raises(TypeError):
        law.calculate_acceleration(test_args.g, test_args.ac, test_args.f, 100)
