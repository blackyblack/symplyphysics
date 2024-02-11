from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal
)

from symplyphysics.laws.thermodynamics import van_der_waals_equation

## Source of numbers: https://alexandr4784.narod.ru/izerp/izerp_61_626.pdf

Args = namedtuple("Args", ["amount_of_substance", "mutual_attraction_constant", "volume", "mutual_repulsion_constant", "pressure"])

@fixture(name="test_args")
def test_args_fixture():
    amount_of_substance = Quantity(0.1093 * units.mole)
    mutual_attraction_constant = Quantity(0.136 * units.pascal * units.meter**6 / units.mole**2)
    volume = Quantity(9e-5 * units.meter**3)
    mutual_repulsion_constant = Quantity(3.16e-5 * units.meter**3 / units.mole)
    pressure = Quantity(2.8e6 * units.pascal)
    return Args(
        amount_of_substance=amount_of_substance,
        mutual_attraction_constant=mutual_attraction_constant,
        volume=volume,
        mutual_repulsion_constant=mutual_repulsion_constant,
        pressure=pressure
    )


def test_basic_law(test_args: Args) -> None:
    result = van_der_waals_equation.calculate_temperature(test_args.amount_of_substance, test_args.mutual_attraction_constant, test_args.volume, test_args.mutual_repulsion_constant, test_args.pressure)
    assert_equal(result, 285.7 * units.kelvin)


def test_bad_amount_of_substance(test_args: Args) -> None:
    bad_amount_of_substance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        van_der_waals_equation.calculate_temperature(bad_amount_of_substance, test_args.mutual_attraction_constant, test_args.volume, test_args.mutual_repulsion_constant, test_args.pressure)
    with raises(TypeError):
        van_der_waals_equation.calculate_temperature(100, test_args.mutual_attraction_constant, test_args.volume, test_args.mutual_repulsion_constant, test_args.pressure)


def test_bad_mutual_attraction_constant(test_args: Args) -> None:
    bad_mutual_attraction_constant = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        van_der_waals_equation.calculate_temperature(test_args.amount_of_substance, test_args.mutual_attraction_constant, test_args.volume, bad_mutual_attraction_constant, test_args.pressure)
    with raises(TypeError):
        van_der_waals_equation.calculate_temperature(test_args.amount_of_substance, 100, test_args.volume, test_args.mutual_repulsion_constant, test_args.pressure)


def test_bad_volume(test_args: Args) -> None:
    bad_volume = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        van_der_waals_equation.calculate_temperature(test_args.amount_of_substance, test_args.mutual_attraction_constant, bad_volume, test_args.mutual_repulsion_constant, test_args.pressure)
    with raises(TypeError):
        van_der_waals_equation.calculate_temperature(test_args.amount_of_substance, test_args.mutual_attraction_constant, 100, test_args.mutual_repulsion_constant, test_args.pressure)


def test_bad_mutual_repulsion_constant(test_args: Args) -> None:
    bad_mutual_repulsion_constant = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        van_der_waals_equation.calculate_temperature(test_args.amount_of_substance, test_args.mutual_attraction_constant, test_args.volume, bad_mutual_repulsion_constant, test_args.pressure)
    with raises(TypeError):
        van_der_waals_equation.calculate_temperature(test_args.amount_of_substance, test_args.mutual_attraction_constant, test_args.volume, 100, test_args.pressure)


def test_bad_pressure(test_args: Args) -> None:
    bad_pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        van_der_waals_equation.calculate_temperature(test_args.amount_of_substance, test_args.mutual_attraction_constant, test_args.volume, test_args.mutual_repulsion_constant, bad_pressure)
    with raises(TypeError):
        van_der_waals_equation.calculate_temperature(test_args.amount_of_substance, test_args.mutual_attraction_constant, test_args.volume, test_args.mutual_repulsion_constant, 100)
