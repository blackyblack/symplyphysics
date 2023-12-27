from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.hydro.buoyant_force_from_volume import calculate_archimedes_force

# Description
# A body submerged in water (1000 kg/m³) with a volume of 0.01 m³ under 9.81 m/s² gravity should experience a buoyant force of 98.1 N.

@fixture(name="test_args")
def test_args_fixture():
    rho = Quantity(1000 * units.kilogram / units.meter**3)
    V = Quantity(0.01 * units.meter**3)
    g = Quantity(9.81 * units.meter / units.second**2)
    Args = namedtuple("Args", ["rho", "V", "g"])
    return Args(rho=rho, V=V, g=g)

def test_archimedes_force(test_args):
    result = calculate_archimedes_force(test_args.rho, test_args.V, test_args.g)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force = convert_to(result, units.newton).evalf(5)
    assert result_force == approx(98.1, 0.001)

def test_bad_density(test_args):
    bad_density = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        calculate_archimedes_force(bad_density, test_args.V, test_args.g)
    with raises(errors.UnitsError):
        calculate_archimedes_force("not a quantity", test_args.V, test_args.g)

def test_bad_volume(test_args):
    bad_volume = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        calculate_archimedes_force(test_args.rho, bad_volume, test_args.g)
    with raises(errors.UnitsError):
        calculate_archimedes_force(test_args.rho, "not a quantity", test_args.g)
