#!/usr/bin/env python3
"""
Calculate the angle of total internal reflection for two given media.

If a ray of light tries to transfer from one medium to another, part of it transfers and refracts,
and another part reflects. There is a possibility for the whole ray to be reflected with nothing
refracted into the other medium. It is possible when the ray transfers from a denser medium to a
less dense one. This phenomenon is utilized by fiber-optics technology.

According to the law, if a1 is the incidence angle and a2 is the refraction angle, then
n1*sin(a1) = n2*sin(a2), or sin(a2) / sin(a1) = n1 / n2. If n1 > n2, a2 will always be greater
than a1. If we increase a1, a2 increases as well. There is such a value of a1 when a2 reaches its
maximum possible value (90 degrees) meaning the refracted ray goes along the border of the media.
If we continue increasing a1, there will be no refraction anymore; the whole ray will be
reflected.
"""

import math
from sympy import solve, pi
from symplyphysics.laws.optics import refraction_angle_from_environments as refraction_law

## Optical fiber usually consists of core with refractive index
CORE_REFRACTIVE_INDEX = 1.479
## and coating with refractive index
COAT_REFRACTIVE_INDEX = 1.474

equation = refraction_law.law.subs({
    refraction_law.incidence_refractive_index: CORE_REFRACTIVE_INDEX,
    refraction_law.resulting_refractive_index: COAT_REFRACTIVE_INDEX,
    refraction_law.refraction_angle: pi / 2
})

solutions = solve(equation, refraction_law.incidence_angle)

angle_expr = solutions[0]
if abs(angle_expr) > pi / 2:
    angle_expr = solutions[1]

assert abs(angle_expr) <= pi / 2

angle_degrees = (angle_expr / math.pi * 180).evalf(3)
print(f"For full internal reflection angle should be above {angle_degrees} degrees")
