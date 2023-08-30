from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import impedance_is_resistance_and_reactance as impedance_def

# Description
## If the object has 10 Ohm resistance in series with capacitor with reactance of -1.59 Ohm,
## the impedance magnitude of the circuit is 10.12 Ohm.
## Taken from https://wiraelectrical.com/impedance-and-admittance/


@fixture(name="test_args")
def test_args_fixture():
    R = Quantity(10 * units.ohm)
    X = Quantity(-1.59 * units.ohm)
    Args = namedtuple("Args", ["R", "X"])
    return Args(R=R, X=X)


def test_basic_impedance(test_args):
    result = impedance_def.calculate_impedance_magnitude(test_args.R, test_args.X)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_impedance = convert_to(result, impedance_def.definition_units_SI).evalf(5)
    assert result_impedance == approx(10.12, 0.001)


def test_bad_resistance(test_args):
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        impedance_def.calculate_impedance_magnitude(Rb, test_args.X)
    with raises(TypeError):
        impedance_def.calculate_impedance_magnitude(100, test_args.X)


def test_bad_reactance(test_args):
    Xb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        impedance_def.calculate_impedance_magnitude(test_args.R, Xb)
    with raises(TypeError):
        impedance_def.calculate_impedance_magnitude(test_args.R, 100)
