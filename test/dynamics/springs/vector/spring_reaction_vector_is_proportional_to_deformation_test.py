from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.core.vectors.vectors import QuantityVector
from symplyphysics.laws.dynamics.springs.vector import spring_reaction_is_proportional_to_deformation as spring_law

Args = namedtuple("Args", ["k", "d", "f"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k = Quantity(0.1 * units.newton / units.meter)
    d_x = Quantity(3 * units.meter)
    d_y = Quantity(1 * units.meter)
    d = QuantityVector([d_x, d_y])
    f_x = Quantity(-0.3 * units.newton)
    f_y = Quantity(-0.1 * units.newton)
    f = QuantityVector([f_x, f_y])
    return Args(k=k, d=d, f=f)


def test_basic_force(test_args: Args) -> None:
    result = spring_law.calculate_force(test_args.k, test_args.d)
    assert_equal(result.components[0], -0.3 * units.newton)
    assert_equal(result.components[1], -0.1 * units.newton)


def test_basic_deformation(test_args: Args) -> None:
    result = spring_law.calculate_deformation(test_args.k, test_args.f)
    assert_equal(result.components[0], 3 * units.meter)
    assert_equal(result.components[1], 1 * units.meter)


def test_bad_elastic_coefficient(test_args: Args) -> None:
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        spring_law.calculate_force(eb, test_args.d)
    with raises(TypeError):
        spring_law.calculate_force(100, test_args.d)
    with raises(errors.UnitsError):
        spring_law.calculate_deformation(eb, test_args.f)
    with raises(TypeError):
        spring_law.calculate_deformation(100, test_args.f)


def test_bad_deformation(test_args: Args) -> None:
    db = Quantity(1 * units.coulomb)
    vb = QuantityVector([db])
    with raises(errors.UnitsError):
        spring_law.calculate_force(test_args.k, vb)
    with raises(TypeError):
        spring_law.calculate_force(test_args.k, 100)


def test_bad_force(test_args: Args) -> None:
    db = Quantity(1 * units.coulomb)
    vb = QuantityVector([db])
    with raises(errors.UnitsError):
        spring_law.calculate_deformation(test_args.k, vb)
    with raises(TypeError):
        spring_law.calculate_deformation(test_args.k, 100)
