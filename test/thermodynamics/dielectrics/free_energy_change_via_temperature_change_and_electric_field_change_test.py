from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.dielectrics import (
    free_energy_change_via_temperature_change_and_electric_field_change as law,)

Args = namedtuple("Args", "s dt e dd")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    s = Quantity(10 * units.joule / units.kelvin / units.meter**3)
    dt = Quantity(-1e-3 * units.kelvin)
    e = Quantity(15 * units.volt / units.meter)
    dd = Quantity(3e-3 * units.coulomb / units.meter**2)
    return Args(s=s, dt=dt, e=e, dd=dd)


def test_law(test_args: Args) -> None:
    result = law.calculate_free_energy_density_change(test_args.s, test_args.dt, test_args.e,
        test_args.dd)
    assert_equal(result, 5.5e-2 * units.joule / units.meter**3)


def test_bad_entropy_density(test_args: Args) -> None:
    sb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_free_energy_density_change(sb, test_args.dt, test_args.e, test_args.dd)
    with raises(TypeError):
        law.calculate_free_energy_density_change(100, test_args.dt, test_args.e, test_args.dd)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_free_energy_density_change(test_args.s, tb, test_args.e, test_args.dd)
    with raises(TypeError):
        law.calculate_free_energy_density_change(test_args.s, 100, test_args.e, test_args.dd)


def test_bad_electric_field(test_args: Args) -> None:
    eb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_free_energy_density_change(test_args.s, test_args.dt, eb, test_args.dd)
    with raises(TypeError):
        law.calculate_free_energy_density_change(test_args.s, test_args.dt, 100, test_args.dd)


def test_bad_electric_displacement(test_args: Args) -> None:
    db = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_free_energy_density_change(test_args.s, test_args.dt, test_args.e, db)
    with raises(TypeError):
        law.calculate_free_energy_density_change(test_args.s, test_args.dt, test_args.e, 100)
