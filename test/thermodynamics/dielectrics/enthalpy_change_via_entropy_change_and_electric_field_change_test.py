from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.dielectrics import (
    enthalpy_change_via_entropy_change_and_electric_field_change as law,
)

Args = namedtuple("Args", "t ds d de")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(300 * units.kelvin)
    ds = Quantity(1e-9 * units.joule / units.kelvin / units.meter**3)
    d = Quantity(1.5e-5 * units.coulomb / units.meter**2)
    de = Quantity(-5e-3 * units.volt / units.meter)
    return Args(t=t, ds=ds, d=d, de=de)


def test_law(test_args: Args) -> None:
    result = law.calculate_enthalpy_density_change(test_args.t, test_args.ds, test_args.d, test_args.de)
    assert_equal(result, 3.75e-7 * units.joule / units.meter**3)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_enthalpy_density_change(tb, test_args.ds, test_args.d, test_args.de)
    with raises(TypeError):
        law.calculate_enthalpy_density_change(100, test_args.ds, test_args.d, test_args.de)


def test_bad_entropy_density(test_args: Args) -> None:
    sb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_enthalpy_density_change(test_args.t, sb, test_args.d, test_args.de)
    with raises(TypeError):
        law.calculate_enthalpy_density_change(test_args.t, 100, test_args.d, test_args.de)


def test_bad_electric_displacement(test_args: Args) -> None:
    db = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_enthalpy_density_change(test_args.t, test_args.ds, db, test_args.de)
    with raises(TypeError):
        law.calculate_enthalpy_density_change(test_args.t, test_args.ds, 100, test_args.de)


def test_bad_electric_field(test_args: Args) -> None:
    eb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_enthalpy_density_change(test_args.t, test_args.ds, test_args.d, eb)
    with raises(TypeError):
        law.calculate_enthalpy_density_change(test_args.t, test_args.ds, test_args.d, 100)
