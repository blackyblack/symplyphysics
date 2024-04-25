from collections import namedtuple
from pytest import fixture, raises
from sympy import S
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    isobaric_molar_heat_capacity_of_ideal_gas_via_adiabatic_index as heat_capacity_law,
)

# Description
## The isobaric molar heat capacity of a monatomic ideal gas (gamma = 5/3) is C_p = 5 cal/(K*mol).

Args = namedtuple("Args", "gamma")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    gamma = S(5)/3
    return Args(gamma=gamma)


def test_law(test_args: Args) -> None:
    result = heat_capacity_law.calculate_isobaric_molar_heat_capacity(test_args.gamma)
    assert_equal(result, 20.8 * units.joule / (units.kelvin * units.mole))


def test_bad_index() -> None:
    gamma_bad = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heat_capacity_law.calculate_isobaric_molar_heat_capacity(gamma_bad)
