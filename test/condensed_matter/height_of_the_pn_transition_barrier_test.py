from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    units,
    errors,
    SI,
    convert_to,
    Quantity,
)
from symplyphysics.laws.condensed_matter import height_of_the_pn_transition_barrier as barrier_law

# Description
## For concentration of donors equal to 3e16 [1 / cm^3] and concentration of acceptors equal to 2e12 [1 / cm^3]
## the height of the potential barrier is known. It is equal to 0.529 volts. You can check it on the website below:
## https://cleanroom.byu.edu/pn_junction
## It is known that the concentration of intrinsic charge carriers in silicon is equal to 1e10 [1 / cm^3].
## https://www.ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html


@fixture(name="test_args")
def test_args_fixture():
    concentration_donors = Quantity(3e16 * (1 / units.centimeter**3))
    concentration_acceptors = Quantity(2e12 * (1 / units.centimeter**3))
    concentration_intrinsic = Quantity(1e10 * (1 / units.centimeter**3))

    temperature = Quantity(300 * units.kelvin)

    Args = namedtuple("Args", ["concentration_donors", "concentration_acceptors", "concentration_intrinsic", "temperature"])
    return Args(
        concentration_donors=concentration_donors,
        concentration_acceptors=concentration_acceptors,
        concentration_intrinsic=concentration_intrinsic,
        temperature=temperature)


def test_basic_height_barrier(test_args):
    result = barrier_law.calculate_height_barrier(test_args.concentration_donors,
        test_args.concentration_acceptors, test_args.concentration_intrinsic, test_args.temperature)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage)
    result = convert_to(result, units.V).evalf(5)
    assert result == approx(0.529, rel=0.1)


def test_bad_concentration_donors(test_args):
    concentration_donors = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        barrier_law.calculate_height_barrier(concentration_donors, test_args.concentration_acceptors, test_args.concentration_intrinsic, test_args.temperature)
    with raises(TypeError):
        barrier_law.calculate_height_barrier(100, test_args.concentration_acceptors, test_args.concentration_intrinsic, test_args.temperature)


def test_bad_concentration_acceptors(test_args):
    concentration_acceptors = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        barrier_law.calculate_height_barrier(test_args.concentration_donors, concentration_acceptors, test_args.concentration_intrinsic, test_args.temperature)
    with raises(TypeError):
        barrier_law.calculate_height_barrier(test_args.concentration_donors, 100, test_args.concentration_intrinsic, test_args.temperature)


def test_bad_concentration_intrinsic(test_args):
    concentration_intrinsic = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        barrier_law.calculate_height_barrier(test_args.concentration_donors, test_args.concentration_acceptors, concentration_intrinsic, test_args.temperature)
    with raises(TypeError):
        barrier_law.calculate_height_barrier(test_args.concentration_donors, test_args.concentration_acceptors, 100, test_args.temperature)


def test_bad_temperature(test_args):
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        barrier_law.calculate_height_barrier(test_args.concentration_donors, test_args.concentration_acceptors, test_args.concentration_intrinsic, temperature)
    with raises(TypeError):
        barrier_law.calculate_height_barrier(test_args.concentration_donors, test_args.concentration_acceptors, test_args.concentration_intrinsic, 100)
