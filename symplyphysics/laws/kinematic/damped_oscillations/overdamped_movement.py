from sympy import Eq

# Description
## ...

# ---------------------------------------------------------------------- #
# z > 1, overdamped system                                               #
# # the system exponentially decays to steady state without oscillations #
# # larger values of the damping ratio return to equilibrium more slowly #
# ---------------------------------------------------------------------- #

# Law: x(t) = exp(-omega*zeta*t) * (C1 * exp(omega*sqrt(zeta**2 - 1)*t) + C2 * exp(-omega*sqrt(zeta**2 - 1)*t))
## x(t) - position of overdamped oscillator
## t - time
## zeta - damping ratio
## omega - undamped angular frequency
