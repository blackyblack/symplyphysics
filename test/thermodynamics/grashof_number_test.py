from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.core.symbols.celsius import Celsius, to_kelvin_quantity
from symplyphysics.laws.thermodynamics import grashof_number


# Example from: https://www.azcalculator.com/calc/grashof-number.php

Args = namedtuple("Args", ["beta", "t_surf", "t_bulk", "l", "mu"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    beta = Quantity(0.000458 * units.kelvin ** -1)
    t_surf = to_kelvin_quantity(Celsius(20))
    t_bulk = to_kelvin_quantity(Celsius(60))
    l = Quantity(1 * units.meter)
    mu = Quantity(0.0000010034 * units.meter**2 / units.second)
    return Args(beta=beta, t_surf=t_surf, t_bulk=t_bulk, l=l, mu=mu)


def test_basic_grashof_number(test_args: Args) -> None:
    result = grashof_number.calculate_grashof_number(
        test_args.beta,
        test_args.t_surf,
        test_args.t_bulk,
        test_args.l,
        test_args.mu,
    )
    assert_equal(result, -178503313966.7169)


def test_bad_coefficient_of_volume_expansion(test_args: Args) -> None:
    beta = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        grashof_number.calculate_grashof_number(
            beta,
            test_args.t_surf,
            test_args.t_bulk,
            test_args.l,
            test_args.mu
        )
    with raises(TypeError):
        grashof_number.calculate_grashof_number(
            100,
            test_args.t_surf,
            test_args.t_bulk,
            test_args.l,
            test_args.mu
        )


def test_bad_temperature(test_args: Args) -> None:
    ts = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        grashof_number.calculate_grashof_number(
            test_args.beta,
            ts,
            test_args.t_bulk,
            test_args.l,
            test_args.mu
        )
    with raises(errors.UnitsError):
        grashof_number.calculate_grashof_number(
            test_args.beta,
            test_args.t_surf,
            ts,
            test_args.l,
            test_args.mu
        )
    with raises(TypeError):
        grashof_number.calculate_grashof_number(
            test_args.beta,
            100,
            test_args.t_bulk,
            test_args.l,
            test_args.mu
        )
    with raises(TypeError):
        grashof_number.calculate_grashof_number(
            test_args.beta,
            test_args.t_surf,
            100,
            test_args.l,
            test_args.mu
        )


def test_bad_length(test_args: Args) -> None:
    l = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        grashof_number.calculate_grashof_number(
            test_args.beta,
            test_args.t_surf,
            test_args.t_bulk,
            l,
            test_args.mu
        )
    with raises(TypeError):
        grashof_number.calculate_grashof_number(
            test_args.beta,
            test_args.t_surf,
            test_args.t_bulk,
            100,
            test_args.mu
        )


def test_bad_viscosity(test_args: Args) -> None:
    mu = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        grashof_number.calculate_grashof_number(
            test_args.beta,
            test_args.t_surf,
            test_args.t_bulk,
            test_args.l,
            mu
        )
    with raises(TypeError):
        grashof_number.calculate_grashof_number(
            test_args.beta,
            test_args.t_surf,
            test_args.t_bulk,
            test_args.l,
            100
        )
