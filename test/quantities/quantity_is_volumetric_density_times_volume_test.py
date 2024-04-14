from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.quantities import quantity_is_volumetric_density_times_volume as density_law

Args = namedtuple("Args", "u q v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    u = Quantity(1 * units.joule / units.liter)
    q = Quantity(1 * units.coulomb / units.liter)
    v = Quantity(1 * units.liter)
    return Args(u=u, q=q, v=v)


def test_volume_energy_density(test_args: Args) -> None:
    result = density_law.calculate_extensive_quantity(test_args.u, test_args.v)
    assert_equal(result, 1 * units.joule)


def test_volume_charge_density(test_args: Args) -> None:
    result = density_law.calculate_extensive_quantity(test_args.q, test_args.v)
    assert_equal(result, 1 * units.coulomb)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        density_law.calculate_extensive_quantity(test_args.u, vb)
    with raises(TypeError):
        density_law.calculate_extensive_quantity(test_args.u, 100)
