from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import density_from_mass_volume


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(1 * units.kilogram)
    V = Quantity(3 * units.meter**3)
    Args = namedtuple("Args", ["m", "V"])
    return Args(m=m, V=V)


def test_basic_density(test_args):
    result = density_from_mass_volume.calculate_density(test_args.m, test_args.V)
    assert_equal(result, 0.3333 * units.kilogram / units.meter**3)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        density_from_mass_volume.calculate_density(mb, test_args.V)
    with raises(TypeError):
        density_from_mass_volume.calculate_density(100, test_args.V)


def test_bad_volume(test_args):
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        density_from_mass_volume.calculate_density(test_args.m, Vb)
    with raises(TypeError):
        density_from_mass_volume.calculate_density(test_args.m, 100)
