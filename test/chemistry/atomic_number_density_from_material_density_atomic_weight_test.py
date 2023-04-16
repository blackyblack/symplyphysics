from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.chemistry import atomic_number_density_from_material_density_atomic_weight as atomic_number_density

@fixture
def test_args():
    # boron carbide density is 2.52 g/cm^3
    p = Quantity(units.mass / units.volume, 2.52 * units.gram / units.centimeter**3)
    # boron carbide atomic mass
    M = Quantity(units.mass / units.amount_of_substance, 55.2 * units.gram / units.mole)
    Args = namedtuple("Args", ["p", "M"])
    return Args(p = p, M = M)

def test_basic_number_density(test_args):
    result = atomic_number_density.calculate_atomic_number_density(test_args.p, test_args.M)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.volume)
    result_number_density = convert_to(result, 1 / units.centimeter**3).subs(units.centimeter, 1).evalf(4)
    assert result_number_density == approx(2.75E+22, 0.01)

def test_bad_density(test_args):
    pb = Quantity(units.length)
    with raises(errors.UnitsError):
        atomic_number_density.calculate_atomic_number_density(pb, test_args.M)
    with raises(TypeError):
        atomic_number_density.calculate_atomic_number_density(100, test_args.M)

def test_bad_atomic_mass(test_args):
    Mb = Quantity(units.length)
    with raises(errors.UnitsError):
        atomic_number_density.calculate_atomic_number_density(test_args.p, Mb)
    with raises(TypeError):
        atomic_number_density.calculate_atomic_number_density(test_args.p, 100)
