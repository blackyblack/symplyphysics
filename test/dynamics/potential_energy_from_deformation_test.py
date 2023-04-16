from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.dynamics import potential_energy_from_deformation as hookes_law

# Description
## If we have spring with elastic koefficient 100 N/m deformated for 2cm, it handles potential energy of 0.02 joules.
## Result is independently calculated with https://www.center-pss.ru/math/raschet-potencialnoi-energii-pruzhini.htm

@fixture
def test_args():
    k = Quantity(units.force / units.length, 100 * units.newton / units.meter)
    x = Quantity(units.length, 2 * units.centimeter)
    Args = namedtuple("Args", ["k", "x"])
    return Args(k=k, x=x)

def test_basic_energy(test_args):
    result = hookes_law.calculate_energy(test_args.k, test_args.x)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).subs(units.joule, 1).evalf(5)
    assert result_energy == approx(0.02, 0.0001)

def test_bad_elastic_koefficient(test_args):
    kb = Quantity(units.length)
    with raises(errors.UnitsError):
        hookes_law.calculate_energy(kb, test_args.x)
    with raises(TypeError):
        hookes_law.calculate_energy(100, test_args.x)

def test_bad_deformation(test_args):
    xb = Quantity(units.charge)
    with raises(errors.UnitsError):
        hookes_law.calculate_energy(test_args.k, xb)
    with raises(TypeError):
        hookes_law.calculate_energy(test_args.k, 100)
