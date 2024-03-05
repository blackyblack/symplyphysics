from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.springs import (
    spring_reaction_is_proportional_to_deformation as hookes_law,
)

# Description
## A spring of stiffness k = 150 N/m is being compressed by 10 cm. The force exerted by the
## spring to counteract the deformation is thus 15 N.

Args = namedtuple("Args", "k x")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k = Quantity(150.0 * units.newton / units.meter)
    x = Quantity(-10.0 * units.centimeter)
    return Args(k=k, x=x)


def test_law(test_args: Args) -> None:
    result = hookes_law.calculate_spring_reaction(test_args.k, test_args.x)
    assert_equal(result, 15.0 * units.newton)


def test_bad_stiffness(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        hookes_law.calculate_spring_reaction(kb, test_args.x)
    with raises(TypeError):
        hookes_law.calculate_spring_reaction(100, test_args.x)


def test_bad_deformation(test_args: Args) -> None:
    xb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        hookes_law.calculate_spring_reaction(test_args.k, xb)
    with raises(TypeError):
        hookes_law.calculate_spring_reaction(test_args.k, 100)
