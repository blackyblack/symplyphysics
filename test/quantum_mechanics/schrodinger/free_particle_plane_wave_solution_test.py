from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.quantum_mechanics.schrodinger import (
    free_particle_plane_wave_solution as solution_law,
)

Args = namedtuple("Args", "p e x t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(0.1 * units.planck_momentum)
    e = Quantity(1 * units.planck_energy)
    x = Quantity(2 * units.planck_length)
    t = Quantity(1 * units.planck_time)
    return Args(p=p, e=e, x=x, t=t)


def test_law(test_args: Args) -> None:
    result = solution_law.calculate_wave_function_value(test_args.p, test_args.e, test_args.x, test_args.t)
    assert_equal(result.real, 0.70, tolerance=5e-3)
    assert_equal(result.imag, -0.72, tolerance=5e-3)


def test_bad_momentum(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        solution_law.calculate_wave_function_value(pb, test_args.e, test_args.x, test_args.t)
    with raises(TypeError):
        solution_law.calculate_wave_function_value(100, test_args.e, test_args.x, test_args.t)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        solution_law.calculate_wave_function_value(test_args.p, eb, test_args.x, test_args.t)
    with raises(TypeError):
        solution_law.calculate_wave_function_value(test_args.p, 100, test_args.x, test_args.t)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        solution_law.calculate_wave_function_value(test_args.p, test_args.e, xb, test_args.t)
    with raises(TypeError):
        solution_law.calculate_wave_function_value(test_args.p, test_args.e, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        solution_law.calculate_wave_function_value(test_args.p, test_args.e, test_args.x, tb)
    with raises(TypeError):
        solution_law.calculate_wave_function_value(test_args.p, test_args.e, test_args.x, 100)
