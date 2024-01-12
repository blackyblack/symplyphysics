#!/usr/bin/env python3

from sympy import symbols, solve
from symplyphysics.laws.optics import lens_focus_from_object_and_image as thin_lens

# Two thin converging lenses of focal lengths 10 cm and 20 cm are separated by 20 cm. 
# An object is placed 15 cm to the left of lens 1. What is the position of the final image?

focal_length_1, focal_length_2 = symbols("focal_length_1 focal_length_2")
distance_to_object_1 = symbols("distance_to_object_1")
distance_between_lenses = symbols("distance_between_lenses")

values = {
    focal_length_1: 10,
    focal_length_2: 20,
    distance_to_object_1: 15,
    distance_between_lenses: 20,
}

# Calculate the position of the image formed by lens 1 using the thin-lens equation

thin_lens_1 = thin_lens.law.subs({
    thin_lens.focus_distance: focal_length_1,
    thin_lens.distance_to_object: distance_to_object_1,
})
distance_to_image_1 = solve(thin_lens_1, thin_lens.distance_to_image)[0]

# Calculate the distance between lens 2 and the image from lens 1

distance_to_object_2 = distance_between_lenses - distance_to_image_1

# Calculate the position of the image formed by lens 2 using the thin-lens equation

thin_lens_2 = thin_lens.law.subs({
    thin_lens.focus_distance: focal_length_2,
    thin_lens.distance_to_object: distance_to_object_2
})
distance_to_image_2 = solve(thin_lens_2, thin_lens.distance_to_image)[0]
distance_to_image_2_value = distance_to_image_2.subs(values).evalf(3)

side = "right" if distance_to_image_2_value > 0 else "left"
print(f"The image is located {abs(distance_to_image_2_value)} cm to the {side} of lens 2.")
