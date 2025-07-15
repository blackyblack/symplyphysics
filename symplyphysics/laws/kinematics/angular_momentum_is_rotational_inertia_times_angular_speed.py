"""
Angular momentum is rotational inertia times angular speed
==========================================================

For a rigid body rotating around a fixed axis, the component of its angular momentum parallel
to the rotational axis is found as the product of the body's rotational inertia and the magnitude
of its angular velocity.

**Conditions:**

#. The body is rigid.
#. The axis of rotation is fixed.
#. The origin of the coordinate system lies on the axis of rotation.

**Links:**

#. `Wikipedia, vector counterpart of this law <https://en.wikipedia.org/wiki/Angular_momentum#Orbital_angular_momentum_in_three_dimensions>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorNorm, VectorDot

from symplyphysics.definitions import (
    rotational_inertia_is_mass_times_squared_radius as _rotational_inertia_def,)
from symplyphysics.definitions.vector import (
    angular_momentum_is_position_cross_linear_momentum as _angular_momentum_def,
    momentum_is_mass_times_velocity_vector as _linear_momentum_def,
)
from symplyphysics.laws.kinematics.vector import (
    absolute_velocity_of_arbitrary_motion as _absolute_velocity_law,
    velocity_of_transfer_between_reference_frames as _transfer_velocity_law,
)

angular_momentum = symbols.angular_momentum
"""
Component of the vector of :symbols:`angular_momentum` parallel to the rotational axis.
"""

rotational_inertia = symbols.rotational_inertia
"""
:symbols:`rotational_inertia` of the body about the given rotational axis.
"""

angular_speed = symbols.angular_speed
"""
:symbols:`angular_speed` of the body.
"""

law = Eq(angular_momentum, rotational_inertia * angular_speed)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law in case of a material point (particle)

# Angular velocity pseudovector `w`.
_w_vec = clone_as_vector_symbol(symbols.angular_speed)

# Unit (pseudo)vector in the direction of `w`.
_e_w_vec = _w_vec / VectorNorm(_w_vec)

# 1. Decompose the position vector of the particle

# `r = r_t + r_n` i.e. the position vector can be split into a component tangential to the angular
# velocity pseudovector and a component normal to it.
_r_t = clone_as_symbol(symbols.distance_to_origin, subscript="t")
_r_t_vec = _r_t * _e_w_vec  # this expresses the tangentiality condition
_r_n_vec = clone_as_vector_symbol(symbols.distance_to_axis, subscript="n")

_orthogonality_condition = {
    VectorDot(_w_vec, _r_n_vec): 0,
}

# NOTE: the origin of the coordinate system must lie on the axis of rotation in order for us to be
# able to do such a decomposition *and* for `VectorNorm(_r_n_vec)` to be the distance from the
# particle to the rotational axis.
_r_vec = _r_t_vec + _r_n_vec

# 2. Prepare the expression for the angular momentum pseudovector

_p_vec = _linear_momentum_def.law.rhs.subs({
    _linear_momentum_def.mass: _rotational_inertia_def.mass,
})

_l_vec = _angular_momentum_def.law.rhs.subs({
    _angular_momentum_def.position_vector: _r_vec,
    _angular_momentum_def.linear_momentum: _p_vec,
})

# 3. Suppose a reference frame is attached to the particle, rotating with it. Then we can use the
#    law of absolute-relative-transfer velocities to find the expression for the particle's
#    velocity in the original, stationary reference frame.

_v_tr_vec = _transfer_velocity_law.law.rhs.subs({
    _transfer_velocity_law.moving_frame_velocity: 0,
    _transfer_velocity_law.angular_velocity: _w_vec,
    _transfer_velocity_law.position_vector: _r_vec,
})

_v_vec = _absolute_velocity_law.law.rhs.subs({
    _absolute_velocity_law.relative_velocity: 0,
    _absolute_velocity_law.transfer_velocity: _v_tr_vec,
})

# 4. Substitute the expressions for the velocity and the rotational inertia.

_l_vec = _l_vec.subs(_linear_momentum_def.velocity, _v_vec).subs(_orthogonality_condition)

_l_t_expr = VectorDot(_l_vec, _e_w_vec).subs(_orthogonality_condition)
_l_t_expr = _l_t_expr.subs({
    VectorNorm(_w_vec): angular_speed,
    VectorNorm(_r_n_vec): _rotational_inertia_def.radial_distance,  # see note in (1)
})

_l_t_expr = solve(
    (Eq(angular_momentum, _l_t_expr), _rotational_inertia_def.definition),
    (angular_momentum, _rotational_inertia_def.mass),
    dict=True,
)[0][angular_momentum]

assert expr_equals(_l_t_expr, law.rhs)

# In case of a finite rigid body, sum the expression for `L_t` across all points of the body; note
# that the body's angular velocity is constant and could be taken out of the summation, thus
# leaving us with the sum `m_i * r_i^2`, equal to the body's rotational inertia about the given
# axis.


@validate_input(rotational_inertia_=rotational_inertia, angular_velocity_=angular_speed)
@validate_output(angular_momentum)
def calculate_angular_momentum(rotational_inertia_: Quantity,
    angular_velocity_: Quantity) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia: rotational_inertia_,
        angular_speed: angular_velocity_,
    })
    return Quantity(result)
