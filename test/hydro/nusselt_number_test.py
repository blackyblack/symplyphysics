from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, Quantity, units, errors
from symplyphysics.laws.hydro import nusselt_number

# Example from: https://www.omnicalculator.com/physics/nusselt-number
# Conditions:
# Forced convection, moderate speed cross-flow of air over a cylinder.
# Heat transfer coefficient equals 200 W / m*m*K (source https://www.engineersedge.com),
# thermal conductivity for air equals 0.022 W / m*K, characteristic length equals 1 m.

Args = namedtuple("Args", ["heat_transfer_coefficient",
                           "characteristic_length", "thermal_conductivity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    heat_transfer_coefficient = Quantity(
        200 * units.watt / (units.meter**2 * units.kelvin))
    characteristic_length = Quantity(1 * units.meter)
    thermal_conductivity = Quantity(
        0.022 * units.watt / (units.meter * units.kelvin))
    return Args(
        heat_transfer_coefficient=heat_transfer_coefficient,
        characteristic_length=characteristic_length,
        thermal_conductivity=thermal_conductivity,
    )


def test_basic_nusselt_number(test_args: Args) -> None:
    result = nusselt_number.calculate_nusselt_number(
        heat_transfer_coefficient_=test_args.heat_transfer_coefficient,
        characteristic_length_=test_args.characteristic_length,
        thermal_conductivity_=test_args.thermal_conductivity
    )
    assert_equal(result, 9090.91)


def test_bad_heat_transfer_coefficient(test_args: Args) -> None:
    htc = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        nusselt_number.calculate_nusselt_number(
            htc,
            test_args.characteristic_length,
            test_args.thermal_conductivity,
        )
    with raises(TypeError):
        nusselt_number.calculate_nusselt_number(
            1,
            test_args.characteristic_length,
            test_args.thermal_conductivity,
        )


def test_bad_length(test_args: Args) -> None:
    bl = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        nusselt_number.calculate_nusselt_number(
            test_args.heat_transfer_coefficient,
            bl,
            test_args.thermal_conductivity,
        )
    with raises(TypeError):
        nusselt_number.calculate_nusselt_number(
            test_args.heat_transfer_coefficient,
            1,
            test_args.thermal_conductivity,
        )


def test_bad_conductivity(test_args: Args) -> None:
    bc = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        nusselt_number.calculate_nusselt_number(
            test_args.heat_transfer_coefficient,
            test_args.characteristic_length,
            bc,
        )
    with raises(TypeError):
        nusselt_number.calculate_nusselt_number(
            test_args.heat_transfer_coefficient,
            test_args.characteristic_length,
            1,
        )
