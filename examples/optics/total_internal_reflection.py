#!/usr/bin/env python3

import math
from sympy import solve, pi
from symplyphysics.laws.optics import refraction_angle_from_environments as refraction_law

# This example show usefulness of refraction law.
## If ray of light tries to transfer from one medium to another, part of it transfers and refracts and another part reflects.
## There is a possibility to have the whole ray to be reflected and nothing refracted to another medium.
## It is possible when ray transfers from more denser medium to less denser.
## This phenomenon is utilized by fiber optics technology.

## As law says, if alpha is accidence angle, beta is refraction angle, n1*sin(alpha) = n2*sin(beta), or sin(beta) / sin(alpha) = n1/n2. If n1 > n2, beta will be always more then alpha.
## If we increase alpha, beta increases as well. There is such value of alpha when beta reaches it's max possible value - 90 degrees, it means refracted ray goes along mediums border.
## If we continue increasing alpha, there will be no more refraction anymore, the whole ray will be only reflected.

solutions = solve(refraction_law.law, refraction_law.incedence_angle, dict=True)
result_expr = solutions[0][refraction_law.incedence_angle]

## Optical fiber usually consists of core with refractive index
core_refractive_index = 1.479
## and coating with refractive index
coat_refractive_index = 1.474

angle_applied = result_expr.subs({
    refraction_law.incedence_refractive_index: core_refractive_index,
    refraction_law.resulting_refractive_index: coat_refractive_index,
    refraction_law.refraction_angle: pi / 2
})

if (angle_applied > pi / 2 or angle_applied < -pi / 2):
    result_expr = solutions[1][refraction_law.incedence_angle]
    angle_applied = result_expr.subs({
        refraction_law.incedence_refractive_index: core_refractive_index,
        refraction_law.resulting_refractive_index: coat_refractive_index,
        refraction_law.refraction_angle: pi / 2
    })

assert angle_applied <= pi / 2
assert angle_applied >= -pi / 2

extreme_angle_degrees = angle_applied.evalf(3) / math.pi * 180
print(
    f"For full internal reflection angle should be above {extreme_angle_degrees.round(3)} degrees")
