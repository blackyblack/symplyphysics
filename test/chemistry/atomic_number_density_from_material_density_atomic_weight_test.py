from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
    assert_approx,
)
from symplyphysics.laws.chemistry import atomic_number_density_from_material_density_atomic_weight as atomic_number_density


@fixture(name="test_args")
def test_args_fixture():
    # boron carbide density is 2.52 g/cm^3
    p = Quantity(2.52 * units.gram / units.centimeter**3)
    # boron carbide atomic mass
    M = Quantity(55.2 * units.gram / units.mole)
    Args = namedtuple("Args", ["p", "M"])
    return Args(p=p, M=M)


def test_basic_number_density(test_args):
    result = atomic_number_density.calculate_atomic_number_density(test_args.p, test_args.M)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.volume)
    result_number_density = convert_to(result, 1 / units.centimeter**3).evalf(6)
    assert_approx(result_number_density, 2.75e22)


def test_bad_density(test_args):
    pb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        atomic_number_density.calculate_atomic_number_density(pb, test_args.M)
    with raises(TypeError):
        atomic_number_density.calculate_atomic_number_density(100, test_args.M)


def test_bad_atomic_mass(test_args):
    Mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        atomic_number_density.calculate_atomic_number_density(test_args.p, Mb)
    with raises(TypeError):
        atomic_number_density.calculate_atomic_number_density(test_args.p, 100)
