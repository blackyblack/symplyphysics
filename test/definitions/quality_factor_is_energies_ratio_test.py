from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_approx, errors, units, convert_to, Quantity, SI)
from symplyphysics.definitions import quality_factor_is_energies_ratio as quality_factor_def

# Description
## If the oscillating system has resonant frequency 2 rad/sec, stores total 13 Joules of energy and dissipates 0.4 Watt power, it's quality factor should be 65.
## No external calculators were used for such computation.


@fixture(name="test_args")
def test_args_fixture():
    w = Quantity(2 * units.radian / units.second)
    W = Quantity(13 * units.joule)
    P = Quantity(0.4 * units.watt)
    Args = namedtuple("Args", ["w", "W", "P"])
    return Args(w=w, W=W, P=P)


def test_basic_quality_factor(test_args):
    result = quality_factor_def.calculate_quality_factor(test_args.w, test_args.W, test_args.P)
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
        quality_factor_def.definition_units_SI)
    result_factor = convert_to(result, quality_factor_def.definition_units_SI).evalf(2)
    assert_approx(result_factor, 65)


def test_bad_frequency(test_args):
    fb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        quality_factor_def.calculate_quality_factor(fb, test_args.W, test_args.P)
    with raises(TypeError):
        quality_factor_def.calculate_quality_factor(100, test_args.W, test_args.P)


def test_bad_energy(test_args):
    eb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        quality_factor_def.calculate_quality_factor(test_args.w, eb, test_args.P)
    with raises(TypeError):
        quality_factor_def.calculate_quality_factor(test_args.w, 100, test_args.P)


def test_bad_power(test_args):
    pb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        quality_factor_def.calculate_quality_factor(test_args.w, test_args.W, pb)
    with raises(TypeError):
        quality_factor_def.calculate_quality_factor(test_args.w, test_args.W, 100)
