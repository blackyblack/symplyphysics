from collections import namedtuple
from sympy import I
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import admittance_is_inversed_impedance as admittance_def

# Description
## If the dipole has 2i Ohm impedance, it should have 0.5i Siemens conductivity. No external calculators were used for such computation.


@fixture(name="test_args")
def test_args_fixture():
    Z = Quantity(2 * I * units.ohm)
    Args = namedtuple("Args", ["Z"])
    return Args(Z=Z)


def test_basic_admittance(test_args):
    result = admittance_def.calculate_admittance(test_args.Z)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.conductance)
    result_admittance = abs(convert_to(result, admittance_def.definition_units_SI).evalf(2))
    assert result_admittance == approx(0.5, 0.001)


def test_bad_impedance():
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        admittance_def.calculate_admittance(Rb)
    with raises(TypeError):
        admittance_def.calculate_admittance(100)
