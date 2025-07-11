from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, Quantity, units, errors
from symplyphysics.laws.hydro import reynolds_number_via_fluid_parameters_and_characteristic_length

# Description
# There is pipe with diameter 0.1 m, density 1000 kg/m3, velocity 1 m/s
# and dynamic viscosity 0.000894 Pa*s. The reynolds number should be 111856.823

Args = namedtuple("Args", ["d", "rho", "v", "mu"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    rho = Quantity(1000 * units.kilogram / units.meter**3)
    d = Quantity(0.1 * units.meter)
    v = Quantity(1 * units.meter / units.second)
    mu = Quantity(0.000894 * units.pascal * units.second)
    return Args(d=d, rho=rho, v=v, mu=mu)


def test_reynolds_number(test_args: Args) -> None:
    result = reynolds_number_via_fluid_parameters_and_characteristic_length.calculate_reynolds_number(test_args.d, test_args.rho,
        test_args.v, test_args.mu)
    assert_equal(result, 111856.823)


def test_bad_velocity(test_args: Args) -> None:
    bv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reynolds_number_via_fluid_parameters_and_characteristic_length.calculate_reynolds_number(test_args.d, test_args.rho, bv,
            test_args.mu)
    with raises(TypeError):
        reynolds_number_via_fluid_parameters_and_characteristic_length.calculate_reynolds_number(test_args.d, test_args.rho, 100,
            test_args.mu)


def test_bad_diameter(test_args: Args) -> None:
    bd = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reynolds_number_via_fluid_parameters_and_characteristic_length.calculate_reynolds_number(bd, test_args.rho, test_args.v,
            test_args.mu)
    with raises(TypeError):
        reynolds_number_via_fluid_parameters_and_characteristic_length.calculate_reynolds_number(100, test_args.rho, test_args.v,
            test_args.mu)


def test_bad_density(test_args: Args) -> None:
    bd = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reynolds_number_via_fluid_parameters_and_characteristic_length.calculate_reynolds_number(test_args.d, bd, test_args.v,
            test_args.mu)
    with raises(TypeError):
        reynolds_number_via_fluid_parameters_and_characteristic_length.calculate_reynolds_number(test_args.d, 10, test_args.v,
            test_args.mu)


def test_bad_dynamic_viscosity(test_args: Args) -> None:
    bd = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reynolds_number_via_fluid_parameters_and_characteristic_length.calculate_reynolds_number(test_args.d, test_args.rho, test_args.v,
            bd)
    with raises(TypeError):
        reynolds_number_via_fluid_parameters_and_characteristic_length.calculate_reynolds_number(test_args.d, test_args.rho, test_args.v,
            0.05)
