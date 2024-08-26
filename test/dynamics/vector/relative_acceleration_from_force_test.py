from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    assert_equal_vectors,
    errors,
    units,
    Quantity,
    QuantityVector,
)
from symplyphysics.laws.dynamics.vector import relative_acceleration_from_force as law

Args = namedtuple("Args", "m a_rel f a_cor a_tr")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(5 * units.kilogram)

    a_unit = units.meter / units.second**2
    a_rel = QuantityVector([6 * a_unit, 1 * a_unit, -5 * a_unit])
    a_cor = QuantityVector([0, 4 * a_unit, -2 * a_unit])
    a_tr = QuantityVector([-2 * a_unit, 3 * a_unit, 0])

    f = QuantityVector([20 * units.newton, 0, -15 * units.newton])

    return Args(m=m, a_rel=a_rel, f=f, a_cor=a_cor, a_tr=a_tr)


def test_relative_law(test_args: Args) -> None:
    result = law.calculate_relative_acceleration(test_args.m, test_args.f, test_args.a_cor,
        test_args.a_tr)
    assert_equal_vectors(result, test_args.a_rel)


def test_force_law(test_args: Args) -> None:
    result_vector = law.force_law(
        test_args.a_rel.to_base_vector(),
        test_args.a_cor.to_base_vector(),
        test_args.a_tr.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(
        result_vector,
        subs={law.mass: test_args.m},
    )
    assert_equal_vectors(result, test_args.f)


def test_mass_law(test_args: Args) -> None:
    result = law.mass_law(
        test_args.a_rel.to_base_vector(),
        test_args.f.to_base_vector(),
        test_args.a_cor.to_base_vector(),
        test_args.a_tr.to_base_vector(),
    )
    assert_equal(result, test_args.m)


def test_coriolis_law(test_args: Args) -> None:
    result_vector = law.coriolis_acceleration_law(
        test_args.a_rel.to_base_vector(),
        test_args.f.to_base_vector(),
        test_args.a_tr.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(
        result_vector,
        subs={law.mass: test_args.m},
    )
    assert_equal_vectors(result, test_args.a_cor)


def test_translation_law(test_args: Args) -> None:
    result_vector = law.translation_acceleration_law(
        test_args.a_rel.to_base_vector(),
        test_args.f.to_base_vector(),
        test_args.a_cor.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(
        result_vector,
        subs={law.mass: test_args.m},
    )
    assert_equal_vectors(result, test_args.a_tr)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(mb, test_args.f, test_args.a_cor, test_args.a_tr)
    with raises(TypeError):
        law.calculate_relative_acceleration(100, test_args.f, test_args.a_cor, test_args.a_tr)


def test_bad_force(test_args: Args) -> None:
    fb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(test_args.m, fb_vector, test_args.a_cor, test_args.a_tr)

    fb_scalar = Quantity(units.newton)
    with raises(AttributeError):
        law.calculate_relative_acceleration(test_args.m, fb_scalar, test_args.a_cor, test_args.a_tr)

    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, 100, test_args.a_cor, test_args.a_tr)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, [100], test_args.a_cor, test_args.a_tr)


def test_bad_acceleration(test_args: Args) -> None:
    ab_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, ab_vector, test_args.a_tr)
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, test_args.a_cor, ab_vector)

    ab_scalar = Quantity(units.meter / units.second**2)
    with raises(AttributeError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, ab_scalar, test_args.a_tr)
    with raises(AttributeError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, test_args.a_cor, ab_scalar)

    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, 100, test_args.a_tr)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, [100], test_args.a_tr)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, test_args.a_cor, 100)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, test_args.a_cor, [100])
