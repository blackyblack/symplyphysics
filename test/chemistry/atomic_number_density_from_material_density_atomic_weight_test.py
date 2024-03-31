from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)
from symplyphysics.laws.chemistry import atomic_number_density_from_material_density_atomic_weight as atomic_number_density

Args = namedtuple("Args", ["p", "M"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # boron carbide density is 2.52 g/cm^3
    p = Quantity(2.52 * units.gram / units.centimeter**3)
    # boron carbide atomic mass
    M = Quantity(55.2 * units.gram / units.mole)
    return Args(p=p, M=M)


def test_basic_number_density(test_args: Args) -> None:
    result = atomic_number_density.calculate_atomic_number_density(test_args.p, test_args.M)
    assert_equal(result, 2.75e22 / units.centimeter**3)


def test_bad_density(test_args: Args) -> None:
    pb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        atomic_number_density.calculate_atomic_number_density(pb, test_args.M)
    with raises(TypeError):
        atomic_number_density.calculate_atomic_number_density(100, test_args.M)


def test_bad_atomic_mass(test_args: Args) -> None:
    Mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        atomic_number_density.calculate_atomic_number_density(test_args.p, Mb)
    with raises(TypeError):
        atomic_number_density.calculate_atomic_number_density(test_args.p, 100)
