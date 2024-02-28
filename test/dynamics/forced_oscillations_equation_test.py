from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import forced_oscillations_equation as forced_eqn

# Description
## An oscillating external force is acting on an oscillator of mass m = 1 kg and natural
## angular frequency w0 = 1 rad/s. 
