from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import potential_energy_from_deformation as hookes_law

# Description
## If we have spring with elastic koefficient 100 N/m deformated for 2cm, it handles potential energy of 0.02 joules.
## Result is independently calculated with https://www.center-pss.ru/math/raschet-potencialnoi-energii-pruzhini.htm

Args = namedtuple("Args", ["k", "x"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k = Quantity(100 * units.newton / units.meter)
    x = Quantity(2 * units.centimeter)
    return Args(k=k, x=x)


def test_basic_energy(test_args: Args) -> None:
    result = hookes_law.calculate_energy(test_args.k, test_args.x)
    assert_equal(result, 0.02 * units.joule)


def test_bad_elastic_koefficient(test_args: Args) -> None:
    kb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hookes_law.calculate_energy(kb, test_args.x)
    with raises(TypeError):
        hookes_law.calculate_energy(100, test_args.x)


def test_bad_deformation(test_args: Args) -> None:
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hookes_law.calculate_energy(test_args.k, xb)
    with raises(TypeError):
        hookes_law.calculate_energy(test_args.k, 100)
