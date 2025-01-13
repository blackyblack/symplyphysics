from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.dielectrics import (
    gibbs_energy_change_via_temperature_change_and_electic_displacement_change as law,)

Args = namedtuple("Args", "s dt d de")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    s = Quantity(10 * units.joule / units.kelvin / units.meter**3)
    dt = Quantity(-1e-3 * units.kelvin)
    d = Quantity(15 * units.coulomb / units.meter**2)
    de = Quantity(3e-3 * units.volt / units.meter)
    return Args(s=s, dt=dt, d=d, de=de)


def test_law(test_args: Args) -> None:
    result = law.calculate_gibbs_energy_density_change(test_args.s, test_args.dt, test_args.d,
        test_args.de)
    assert_equal(result, -3.5e-2 * units.joule / units.meter**3)


def test_bad_entropy_density(test_args: Args) -> None:
    sb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_gibbs_energy_density_change(sb, test_args.dt, test_args.d, test_args.de)
    with raises(TypeError):
        law.calculate_gibbs_energy_density_change(100, test_args.dt, test_args.d, test_args.de)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_gibbs_energy_density_change(test_args.s, tb, test_args.d, test_args.de)
    with raises(TypeError):
        law.calculate_gibbs_energy_density_change(test_args.s, 100, test_args.d, test_args.de)


def test_bad_electric_displacement(test_args: Args) -> None:
    db = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_gibbs_energy_density_change(test_args.s, test_args.dt, db, test_args.de)
    with raises(TypeError):
        law.calculate_gibbs_energy_density_change(test_args.s, test_args.dt, 100, test_args.de)


def test_bad_electric_field_strength(test_args: Args) -> None:
    eb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_gibbs_energy_density_change(test_args.s, test_args.dt, test_args.d, eb)
    with raises(TypeError):
        law.calculate_gibbs_energy_density_change(test_args.s, test_args.dt, test_args.d, 100)
