#!/usr/bin/env python3

from sympy import symbols, solve
from symplyphysics import units
from symplyphysics.laws.optics import lens_focus_from_object_and_image as thin_lens

# Two thin converging lenses of focal lengths F1 = 10 cm and F2 = 20 cm are separated
# by 20 cm. An object is placed 15 cm to the left of the lens 1. What is the position 
# of the final image?

F1, F2 = symbols("F1 F2")  # focal lengths
d1 = symbols("d1")  # distance from lens 1 to object
delta = symbols("delta")  # distance between the two lenses

values = {
    F1: 10,
    F2: 20,
    d1: 15,
    delta: 20,
}

# Calculate the position of the image formed by lens 1 using the thin-lens equation

thin_lens_1 = thin_lens.law.subs({
    thin_lens.focus_distance: F1,
    thin_lens.distance_to_object: d1,
})
f1 = solve(thin_lens_1, thin_lens.distance_to_image)[0]

# Calculate the distance between lens 2 and the image from lens 1

d2 = delta - f1

# Calculate the position of the image formed by lens 2 using the thin-lens equation

thin_lens_2 = thin_lens.law.subs({
    thin_lens.focus_distance: F2,
    thin_lens.distance_to_object: d2
})
f2 = solve(thin_lens_2, thin_lens.distance_to_image)[0]

result_f2 = f2.subs(values)

side = "right" if result_f2 > 0 else "left"
print(f"The image is located {abs(result_f2.evalf(3))} cm to the {side} of lens 2.")
