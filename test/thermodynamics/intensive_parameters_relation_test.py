from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics import intensive_parameters_relation as gibbs_duhem_law

# Description
## The change in chemical potential of the system with entropy S = 10 J/K, volume V = 1 m**3 and 1000 particles
## when the pressure change is dp = 0.1 Pa and the temperature change dT = 0.1 K, is d(mu) = -0.9 mJ.

Args = namedtuple("Args", "s dt v dp n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    s = Quantity(10 * units.joule / units.kelvin)
    dt = Quantity(0.1 * units.kelvin)
    v = Quantity(1 * units.meter**3)
    dp = Quantity(0.1 * units.pascal)
    n = 1000
    return Args(s=s, dt=dt, v=v, dp=dp, n=n)


def test_law(test_args: Args) -> None:
    result = gibbs_duhem_law.calculate_chemical_potential_change(test_args.s, test_args.dt,
        test_args.v, test_args.dp, test_args.n)
    assert_equal(result, -0.9 * prefixes.milli * units.joule)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_duhem_law.calculate_chemical_potential_change(sb, test_args.dt, test_args.v,
            test_args.dp, test_args.n)
    with raises(TypeError):
        gibbs_duhem_law.calculate_chemical_potential_change(100, test_args.dt, test_args.v,
            test_args.dp, test_args.n)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_duhem_law.calculate_chemical_potential_change(test_args.s, tb, test_args.v,
            test_args.dp, test_args.n)
    with raises(TypeError):
        gibbs_duhem_law.calculate_chemical_potential_change(test_args.s, 100, test_args.v,
            test_args.dp, test_args.n)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_duhem_law.calculate_chemical_potential_change(test_args.s, test_args.dt, vb,
            test_args.dp, test_args.n)
    with raises(TypeError):
        gibbs_duhem_law.calculate_chemical_potential_change(test_args.s, test_args.dt, 100,
            test_args.dp, test_args.n)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_duhem_law.calculate_chemical_potential_change(test_args.s, test_args.dt, test_args.v,
            pb, test_args.n)
    with raises(TypeError):
        gibbs_duhem_law.calculate_chemical_potential_change(test_args.s, test_args.dt, test_args.v,
            100, test_args.n)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_duhem_law.calculate_chemical_potential_change(test_args.s, test_args.dt, test_args.v,
            test_args.dp, nb)
