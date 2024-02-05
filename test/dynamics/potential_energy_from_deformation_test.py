from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import potential_energy_from_deformation as hookes_law

# Description
## If we have spring with elastic koefficient 100 N/m deformated for 2cm, it handles potential energy of 0.02 joules.
## Result is independently calculated with https://www.center-pss.ru/math/raschet-potencialnoi-energii-pruzhini.htm


@fixture(name="test_args")
def test_args_fixture():
    k = Quantity(100 * units.newton / units.meter)
    x = Quantity(2 * units.centimeter)
    Args = namedtuple("Args", ["k", "x"])
    return Args(k=k, x=x)


def test_basic_energy(test_args):
    result = hookes_law.calculate_energy(test_args.k, test_args.x)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).evalf(5)
    assert_approx(result_energy, 0.02)


def test_bad_elastic_koefficient(test_args):
    kb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hookes_law.calculate_energy(kb, test_args.x)
    with raises(TypeError):
        hookes_law.calculate_energy(100, test_args.x)


def test_bad_deformation(test_args):
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hookes_law.calculate_energy(test_args.k, xb)
    with raises(TypeError):
        hookes_law.calculate_energy(test_args.k, 100)
