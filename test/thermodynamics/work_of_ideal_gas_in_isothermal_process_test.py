from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal
)

from symplyphysics.laws.thermodynamics import work_of_ideal_gas_in_isothermal_process as work_of_ideal_gas

## Source of numbers: https://easy-physic.ru/category/physic/ege1/molekulyarno-kineticheskaya-teoriya/termodinamika/rabota-gaza/rabota-gaza-slozhnye-zadachi

@fixture(name="test_args")
def test_args_fixture():
    mass = Quantity(0.005 * units.kilogram)
    molar_mass = Quantity(0.002 * units.kilogram / units.mole)
    start_volume = Quantity(0.22 * units.meter**3)
    final_volume = Quantity(0.66 * units.meter**3)
    temperature = Quantity(284 * units.kelvin)
    Args = namedtuple("Args", ["mass", "molar_mass", "start_volume", "final_volume", "temperature"])
    return Args(
        mass=mass,
        molar_mass=molar_mass,
        start_volume=start_volume,
        final_volume=final_volume,
        temperature=temperature,
    )


def test_basic_law(test_args):
    result = work_of_ideal_gas.calculate_work(test_args.mass, test_args.molar_mass, test_args.start_volume, test_args.final_volume, test_args.temperature)
    assert_equal(result, 6481.8 * units.joule)


def test_bad_mass(test_args):
    bad_mass = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        work_of_ideal_gas.calculate_work(bad_mass, test_args.molar_mass, test_args.start_volume, test_args.final_volume, test_args.temperature)
    with raises(TypeError):
        work_of_ideal_gas.calculate_work(100, test_args.molar_mass, test_args.start_volume, test_args.final_volume, test_args.temperature)

def test_bad_volume(test_args):
    bad_volume = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        work_of_ideal_gas.calculate_work(test_args.mass, test_args.molar_mass, bad_volume, test_args.final_volume, test_args.temperature)
    with raises(TypeError):
        work_of_ideal_gas.calculate_work(test_args.mass, test_args.molar_mass, 100, test_args.final_volume, test_args.temperature)
    with raises(errors.UnitsError):
        work_of_ideal_gas.calculate_work(test_args.mass, test_args.molar_mass, test_args.start_volume, bad_volume, test_args.temperature)
    with raises(TypeError):
        work_of_ideal_gas.calculate_work(test_args.mass, test_args.molar_mass, test_args.start_volume, 100, test_args.temperature)


def test_bad_temperature(test_args):
    bad_temperature = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        work_of_ideal_gas.calculate_work(test_args.mass, test_args.molar_mass, test_args.start_volume, test_args.final_volume,  bad_temperature)
    with raises(TypeError):
        work_of_ideal_gas.calculate_work(100, test_args.molar_mass, test_args.start_volume, test_args.final_volume, 100)
