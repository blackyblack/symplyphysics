from sympy import Symbol, solve, Eq, simplify

from symplyphysics import print_expression, Quantity, prefixes, units, \
    convert_to
from symplyphysics.laws.optics import lens_focus_from_object_and_image as focus_law
from symplyphysics.laws.optics import focal_length_of_a_concave_spherical_mirror as mirror_law

# Example from https://www.test-uz.ru/solutions_view.php?cat=6&num=15.5

# In a concave mirror with a radius of curvature R = 40 cm, they want to get a real image,
# the height of which is half the height of the object itself.
# Where should the object be placed and where will the image be obtained?

heigth_object = Symbol("heigth_object")
heigth_image = Symbol("heigth_image")

# Because it said in example question
heigth_object_eq = Eq(heigth_object, 2 * heigth_image)

focus = Symbol("focus")
radius_curvature = Symbol("radius_curvature")

distance_from_object = Symbol("distance_from_object")
distance_from_image = Symbol("distance_from_image")

heigth_object_value = heigth_object_eq.rhs

# From picture in example:
# heigth_image / heigth_object  = focus / (distance_from_object - focus)
similarity_sides_triangle_eq = Eq(heigth_image / heigth_object, focus / (distance_from_object - focus))

focus_eq = focus_law.law.subs({
    focus_law.focus_distance: focus,
    focus_law.distance_to_image: distance_from_image,
    focus_law.distance_to_object: distance_from_object
})
focus_value = solve(focus_eq, focus, dict=True)[0][focus]

similarity_sides_triangle_simplify_eq = simplify(similarity_sides_triangle_eq.subs({
    focus: focus_value,
    heigth_object: heigth_object_value
}))
distance_lens_to_object_value = solve(similarity_sides_triangle_simplify_eq, distance_from_object,
                                      dict=True)[0][distance_from_object]

updated_focus_value = focus_value.subs({
    distance_from_object: distance_lens_to_object_value
})

focus_lens_eq = mirror_law.law.subs({
    mirror_law.focus_distance: updated_focus_value,
    mirror_law.curvature_radius: radius_curvature
})
distance_from_image_value = solve(focus_lens_eq, distance_from_image, dict=True)[0][distance_from_image]
distance_from_object_value = distance_lens_to_object_value.subs({distance_from_image: distance_from_image_value})

answer1 = Eq(distance_from_object, distance_from_object_value)
answer2 = Eq(distance_from_image, distance_from_image_value)
print(f"Equation for distance from object to lens is:\n{print_expression(distance_from_object_value)}")
print(f"Equation for distance from image to lens is:\n{print_expression(distance_from_image_value)}")

distance_from_object_m = distance_from_object_value.subs({
    radius_curvature: Quantity(20 * prefixes.centi * units.meters)
})
distance_from_image_m = distance_from_image_value.subs({
    radius_curvature: Quantity(20 * prefixes.centi * units.meters)
})
print(f"Distance from object to lens is: {print_expression(convert_to(Quantity(distance_from_object_m), prefixes.centi * units.meters).evalf(3))} cm")
print(f"Distance from image to lens is: {print_expression(convert_to(Quantity(distance_from_image_m), prefixes.centi * units.meters).evalf(3))} cm")
