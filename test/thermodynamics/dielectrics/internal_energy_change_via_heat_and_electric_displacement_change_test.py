from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.dielectrics import (
    internal_energy_change_via_heat_and_electric_displacement_change as law,)

Args = namedtuple("Args", "dq e dd")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dq = Quantity(1e-3 * units.joule / units.meter**3)
    e = Quantity(10 * units.volt / units.meter)
    dd = Quantity(4e-3 * units.coulomb / units.meter**2)
    return Args(dq=dq, e=e, dd=dd)


def test_law(test_args: Args) -> None:
    result = law.calculate_internal_energy_density_change(test_args.dq, test_args.e, test_args.dd)
    assert_equal(result, 4.1e-2 * units.joule / units.meter**3, tolerance=5e-3)


def test_bad_heat(test_args: Args) -> None:
    qb = Quantity(units.second)
    with raises(errors.UnitsError):
        law.calculate_internal_energy_density_change(qb, test_args.e, test_args.dd)
    with raises(TypeError):
        law.calculate_internal_energy_density_change(100, test_args.e, test_args.dd)


def test_bad_electric_field(test_args: Args) -> None:
    eb = Quantity(units.second)
    with raises(errors.UnitsError):
        law.calculate_internal_energy_density_change(test_args.dq, eb, test_args.dd)
    with raises(TypeError):
        law.calculate_internal_energy_density_change(test_args.dq, 100, test_args.dd)


def test_bad_electric_displacement(test_args: Args) -> None:
    db = Quantity(units.second)
    with raises(errors.UnitsError):
        law.calculate_internal_energy_density_change(test_args.dq, test_args.e, db)
    with raises(TypeError):
        law.calculate_internal_energy_density_change(test_args.dq, test_args.e, 100)
