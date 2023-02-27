from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)

from symplyphysics.laws.waves import frequency_from_velocity as dopler_law

# Man stands near railroad and hears horn of incoming train. Train comes with velocity of 60km/h and horns with 2000Hz frequency. What frequency hears the man? 
## We have online calc for Dopler effect here: https://planetcalc.ru/2351/. With our parameters we should obtain 1907Hz observed frequency.

@fixture
def test_args():
    sound_velocity = units.Quantity('sound_velocity')
    SI.set_quantity_dimension(sound_velocity, units.velocity)
    SI.set_quantity_scale_factor(sound_velocity, 340 * units.meter/units.second)

    train_velocity = units.Quantity('train_velocity')
    SI.set_quantity_dimension(train_velocity, units.velocity)
    SI.set_quantity_scale_factor(train_velocity, 60 * units.kilometer/units.hour)

    man_velocity = units.Quantity('man_velocity')
    SI.set_quantity_dimension(man_velocity, units.velocity)
    SI.set_quantity_scale_factor(man_velocity, 0)

    horn_frequency = units.Quantity('horn_frequency')
    SI.set_quantity_dimension(horn_frequency, units.frequency)
    SI.set_quantity_scale_factor(horn_frequency, 2000 * units.hertz)
    
    Args = namedtuple('Args', ['sound_velocity', 'man_velocity', 'train_velocity', 'horn_frequency'])
    return Args(sound_velocity = sound_velocity, man_velocity = man_velocity, train_velocity = train_velocity, horn_frequency = horn_frequency)

def test_basic_frequency(test_args):
    result = dopler_law.calculate_frequency(test_args.horn_frequency, test_args.sound_velocity, test_args.train_velocity, test_args.man_velocity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.frequency)
    result_focus = convert_to(result, units.hertz).subs(units.hertz, 1).evalf(4)
    assert result_focus == approx(1907, 0.1)


def test_bad_velocity(test_args):
    vb = units.Quantity('vb')
    SI.set_quantity_dimension(vb, units.charge)
    SI.set_quantity_scale_factor(vb, 1 * units.coulomb)

    with raises(errors.UnitsError):
        dopler_law.calculate_frequency(test_args.horn_frequency, vb, test_args.train_velocity, test_args.man_velocity)

    with raises(TypeError):
        dopler_law.calculate_frequency(test_args.horn_frequency, 100, test_args.train_velocity, test_args.man_velocity)

    with raises(errors.UnitsError):
        dopler_law.calculate_frequency(test_args.horn_frequency, test_args.sound_velocity, vb, test_args.man_velocity)

    with raises(TypeError):
        dopler_law.calculate_frequency(test_args.horn_frequency, test_args.sound_velocity, 100, test_args.man_velocity)

    with raises(errors.UnitsError):
        dopler_law.calculate_frequency(test_args.horn_frequency, test_args.sound_velocity, test_args.train_velocity, vb)

    with raises(TypeError):
        dopler_law.calculate_frequency(test_args.horn_frequency, test_args.sound_velocity, test_args.train_velocity, 100)


def test_bad_frequency(test_args):
    fb = units.Quantity('fb')
    SI.set_quantity_dimension(fb, units.charge)
    SI.set_quantity_scale_factor(fb, 1 * units.coulomb)

    with raises(errors.UnitsError):
        dopler_law.calculate_frequency(fb, test_args.sound_velocity, test_args.train_velocity, test_args.man_velocity)

    with raises(TypeError):
        dopler_law.calculate_frequency(100, test_args.sound_velocity, test_args.train_velocity, test_args.man_velocity)
        