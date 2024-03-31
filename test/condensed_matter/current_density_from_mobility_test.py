from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.condensed_matter import current_density_from_mobility as current_law

# Description
## The current density values for a silicon conductor are known.
## It is equal to 120 [A / сm^2].
## https://homework.study.com/explanation/a-p-type-silicon-sample-is-required-to-have-a-drift-current-density-of-120-a-cm-2-under-an-applied-electric-field-of-20-v-cm-determine-the-doping-density-to-meet-this-specification.html
## https://www.ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html

Args = namedtuple("Args", [
    "electrons_concentration", "holes_concentration", "electrons_mobility", "holes_mobility",
    "electric_intensity"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    electrons_concentration = Quantity(1e5 * (1 / units.centimeter**3))
    holes_concentration = Quantity(7.81e16 * (1 / units.centimeter**3))
    electrons_mobility = Quantity(0.17 * (units.meter**2 / units.volt / units.second))
    holes_mobility = Quantity(0.048 * (units.meter**2 / units.volt / units.second))
    electric_intensity = Quantity(2000 * (units.volt / units.meter))

    return Args(electrons_concentration=electrons_concentration,
        holes_concentration=holes_concentration,
        electrons_mobility=electrons_mobility,
        holes_mobility=holes_mobility,
        electric_intensity=electric_intensity)


def test_basic_current(test_args: Args) -> None:
    result = current_law.calculate_current_density(test_args.electrons_concentration,
        test_args.holes_concentration, test_args.electrons_mobility, test_args.holes_mobility,
        test_args.electric_intensity)
    assert_equal(result, 120 * units.ampere / units.centimeter**2)


def test_bad_electrons_concentration(test_args: Args) -> None:
    electrons_concentration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current_density(electrons_concentration,
            test_args.holes_concentration, test_args.electrons_mobility, test_args.holes_mobility,
            test_args.electric_intensity)
    with raises(TypeError):
        current_law.calculate_current_density(100, test_args.holes_concentration,
            test_args.electrons_mobility, test_args.holes_mobility, test_args.electric_intensity)


def test_bad_holes_concentration(test_args: Args) -> None:
    holes_concentration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current_density(test_args.electrons_concentration,
            holes_concentration, test_args.electrons_mobility, test_args.holes_mobility,
            test_args.electric_intensity)
    with raises(TypeError):
        current_law.calculate_current_density(test_args.electrons_concentration, 100,
            test_args.electrons_mobility, test_args.holes_mobility, test_args.electric_intensity)


def test_bad_electrons_mobility(test_args: Args) -> None:
    electrons_mobility = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current_density(test_args.electrons_concentration,
            test_args.holes_concentration, electrons_mobility, test_args.holes_mobility,
            test_args.electric_intensity)
    with raises(TypeError):
        current_law.calculate_current_density(test_args.electrons_concentration,
            test_args.holes_concentration, 100, test_args.holes_mobility,
            test_args.electric_intensity)


def test_bad_holes_mobility(test_args: Args) -> None:
    holes_mobility = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current_density(test_args.electrons_concentration,
            test_args.holes_concentration, test_args.electrons_mobility, holes_mobility,
            test_args.electric_intensity)
    with raises(TypeError):
        current_law.calculate_current_density(test_args.electrons_concentration,
            test_args.holes_concentration, test_args.electrons_mobility, 100,
            test_args.electric_intensity)


def test_bad_electric_intensity(test_args: Args) -> None:
    electric_intensity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current_density(test_args.electrons_concentration,
            test_args.holes_concentration, test_args.electrons_mobility, test_args.holes_mobility,
            electric_intensity)
    with raises(TypeError):
        current_law.calculate_current_density(test_args.electrons_concentration,
            test_args.holes_concentration, test_args.electrons_mobility, test_args.holes_mobility,
            100)
