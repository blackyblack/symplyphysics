from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.chemistry import lennard_jones_potential

# Description
## The value of the Lennard-Jones potential for two interacting Xenon atoms at the distance of
## 1 Å is 2.65e-13 J (e = 2.94e-21 J, sigma = 4.1 Å)


@fixture(name="test_args")
def test_args_fixture():
    e = Quantity(2.94e-21 * units.joule)
    sigma = Quantity(4.10e-10 * units.angstrom)
    r = Quantity(1.0 * units.angstrom)
    Args = namedtuple("Args", "e sigma r")
    return Args(e=e, sigma=sigma, r=r)


def test_basic_law(test_args):
    result = lennard_jones_potential.calculate_potential(test_args.e, test_args.sigma, test_args.r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_value = convert_to(result, units.joule).evalf(3)
    assert result_value == approx(2.65e-13, 1e-15)


def test_bad_energy(test_args):
    eb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        lennard_jones_potential.calculate_potential(eb, test_args.sigma, test_args.r)
    with raises(TypeError):
        lennard_jones_potential.calculate_potential(100, test_args.sigma, test_args.r)


def test_bad_distances(test_args):
    rb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        lennard_jones_potential.calculate_potential(test_args.e, rb, test_args.r)
    with raises(errors.UnitsError):
        lennard_jones_potential.calculate_potential(test_args.e, test_args.sigma, rb)
    with raises(TypeError):
        lennard_jones_potential.calculate_potential(test_args.e, 100, test_args.r)
    with raises(TypeError):
        lennard_jones_potential.calculate_potential(test_args.e, test_args.sigma, 100)
