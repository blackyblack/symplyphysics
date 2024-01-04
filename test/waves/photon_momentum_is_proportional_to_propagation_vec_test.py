from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import pi
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.waves import photon_momentum_is_proportional_to_propagation_vec as photon_momentum_law

# Description
## Assert we have ultraviolet radiation with frequency of 3e16 Hz.
## The value of the modulus of the wave propagation vector for a given frequency in a vacuum will be 2pi*1e8.
## With online calculator
## https://www.fxyz.ru/%D1%84%D0%BE%D1%80%D0%BC%D1%83%D0%BB%D1%8B_%D0%BF%D0%BE_%D1%84%D0%B8%D0%B7%D0%B8%D0%BA%D0%B5/%D0%B0%D1%82%D0%BE%D0%BC%D0%BD%D0%B0%D1%8F_%D1%84%D0%B8%D0%B7%D0%B8%D0%BA%D0%B0/%D0%BA%D0%B2%D0%B0%D0%BD%D1%82%D1%8B/%D1%84%D0%BE%D1%82%D0%BE%D0%BD/%D0%B8%D0%BC%D0%BF%D1%83%D0%BB%D1%8C%D1%81_%D1%84%D0%BE%D1%82%D0%BE%D0%BD%D0%B0/
## we obtain momentum of single photone equal to 6.6307632061911e-26 Newton*second.


@fixture(name="test_args")
def test_args_fixture():
    module_propagation_vector = Quantity(2e8 * pi * (1/units.meter))
    Args = namedtuple("Args", ["module_propagation_vector"])
    return Args(module_propagation_vector=module_propagation_vector)


def test_basic_momentum(test_args):
    result = photon_momentum_law.calculate_momentum(test_args.module_propagation_vector)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.momentum)
    result_current = convert_to(result, units.newton * units.second).evalf(6)
    assert result_current == approx(6.6307632061911e-26, 0.00001)


def test_bad_propogation_vec():
    module_propagation_vector = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        photon_momentum_law.calculate_momentum(module_propagation_vector)
    with raises(TypeError):
        photon_momentum_law.calculate_momentum(100)
