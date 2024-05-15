from sympy import Eq, sqrt
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals

# Description
## Let us consider two inertial reference frames `S` (lab frame) and `S'` (proper frame). The proper frame moves
## with speed `v` relative to the lab frame. Both frames contain identical fixed (relative to their own frame)
## clocks that are synchronized according to the [Einstein rule](https://en.wikipedia.org/wiki/Einstein_synchronisation).
## Let `x, y, z, t` be the coordinates and time of some event in frame `S`, and `x', y', z', t'` be the coordinates
## and time of the same event in frame `S'`. Assuming that the space is uniform and isotropic and that the time is
## uniform, there exists a linear dependence between `x, y, z, t` and `x', y', z', t'`, which is called the Lorentz
## transformation of space and time.

# Law: t' = (t - (V * x) / c**2) / sqrt(1 - (V / c)**2)
## t' - time in proper (S') reference frame
## t - time in lab (S) reference frame
## x - position in lab (S) reference fram
## V - speed of proper (S') frame relative to lab (S') frame
## c - speed of light

# Conditions
## - Space is uniform and isotropic.
## - Time is uniform.
## - The relative frame velocity is parallel to the `x` axis.

# Notes
## - In this law, the Lorentz transformation from the lab frame `S` into the proper frame `S'` is described. In order to get
##   an opposite transformation (from the proper frame `S'` into the lab frame `S`), replace all primed variables with unprimed
##   ones and vice verce, and replace `V` with `-V`. This is consistent with the fact that frame `S` can be viewed as moving
##   with speed `-V` relative to frame `S'`, and hence the same Lorentz transformation can be applied.
## - In the limit of `V/c << 1` the formula reduces to the classical Galilean transformation `t' = t`.

time_in_proper_frame = Symbol("time_in_proper_frame", units.time)
time_in_lab_frame = Symbol("time_in_lab_frame", units.time)
position_in_lab_frame = Symbol("position_in_lab_frame", units.length)
relative_frame_speed = Symbol("relative_frame_speed", units.velocity)  # speed of proper frame relative to lab frame

law = Eq(
    time_in_proper_frame,
    (time_in_lab_frame - relative_frame_speed * position_in_lab_frame / speed_of_light**2)
    / sqrt(1 - (relative_frame_speed / speed_of_light)**2)
)

# In the limit `V/c << 1`, the formula reduces to the classical relation `t' = t`.
_classical_time_expr = law.rhs.series(relative_frame_speed, 0, 1).removeO()
assert expr_equals(_classical_time_expr, time_in_lab_frame)


@validate_input(
    time_in_lab_frame_=time_in_lab_frame,
    position_in_lab_frame_=position_in_lab_frame,
    relative_frame_speed_=relative_frame_speed,
)
@validate_output(time_in_proper_frame)
def calculate_time_in_proper_frame(
    time_in_lab_frame_: Quantity,
    position_in_lab_frame_: Quantity,
    relative_frame_speed_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        time_in_lab_frame: time_in_lab_frame_,
        position_in_lab_frame: position_in_lab_frame_,
        relative_frame_speed: relative_frame_speed_,
    })
    return Quantity(result)
