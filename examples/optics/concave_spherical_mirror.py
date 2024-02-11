from sympy import Symbol, solve, Eq, simplify

from symplyphysics import print_expression, Quantity, prefixes, units, convert_to
from symplyphysics.laws.optics import focal_length_of_a_concave_spherical_mirror as mirror_law
from symplyphysics.laws.optics import linear_magnification_from_distance_to_object_and_distance_to_image as magnification_distance_law
from symplyphysics.laws.optics import linear_magnification_from_object_height_and_image_height as magnification_height_law
from symplyphysics.laws.optics import lens_focus_from_object_and_image as lens_focus_law

# Example from https://www.test-uz.ru/solutions_view.php?cat=6&num=15.5

# In a concave mirror with a radius of curvature R = 40 cm, they want to get a real image,
# the height of which is half the height of the object itself.
# Where should the object be placed and where will the image be obtained?

linear_magnification = Symbol("linear_magnification")
radius_curvature = Symbol("radius_curvature")

focus = Symbol("focus")

distance_from_object = Symbol("distance_from_object")
distance_from_image = Symbol("distance_from_image")
object_height = Symbol("object_height")
image_height = Symbol("image_height")

magnification_height_eq = magnification_height_law.law.subs({
    magnification_height_law.image_height: image_height,
    magnification_height_law.object_height: object_height,
    magnification_height_law.magnification: linear_magnification,
})

magnification_distance_eq = magnification_distance_law.law.subs({
    magnification_distance_law.distance_to_image: distance_from_image,
    magnification_distance_law.distance_to_object: distance_from_object,
    magnification_distance_law.magnification: linear_magnification,
})

distance_from_object_value = solve(
    [magnification_height_eq, magnification_distance_eq],
    (image_height/object_height, distance_from_object), dict=True)[0][distance_from_object]

focus_eq = mirror_law.law.subs({
    mirror_law.curvature_radius: radius_curvature,
    mirror_law.focus_distance: focus
})

lens_focus_eq = lens_focus_law.law.subs({
    lens_focus_law.focus_distance: focus,
    lens_focus_law.distance_to_image: distance_from_image,
    lens_focus_law.distance_to_object: distance_from_object_value,
})

distance_from_image_value = solve(
    [focus_eq, lens_focus_eq],
    (focus, distance_from_image), dict=True)[0][distance_from_image]

magnification_eq = magnification_distance_law.law.subs({
    magnification_distance_law.magnification: linear_magnification,
    magnification_distance_law.distance_to_object: distance_from_object,
    magnification_distance_law.distance_to_image: distance_from_image_value
})
distance_from_object_value = simplify(
    solve(magnification_eq, distance_from_object, dict=True)[0][distance_from_object])

answer1 = Eq(distance_from_object, distance_from_object_value)
answer2 = Eq(distance_from_image, distance_from_image_value)
print(f"Equation for distance from object to mirror is:\n{print_expression(answer1)}")
print(f"Equation for distance from image to mirror is:\n{print_expression(answer2)}")

distance_from_object_m = distance_from_object_value.subs({
    radius_curvature: Quantity(40 * prefixes.centi * units.meters),
    linear_magnification: 1 / 2
})
distance_from_image_m = distance_from_image_value.subs({
    radius_curvature: Quantity(40 * prefixes.centi * units.meters),
    linear_magnification: 1 / 2
})
print(f"Distance from object to mirror is: {print_expression(convert_to(Quantity(distance_from_object_m), prefixes.centi * units.meters).evalf(3))} cm")
print(f"Distance from image to mirror is: {print_expression(convert_to(Quantity(distance_from_image_m), prefixes.centi * units.meters).evalf(3))} cm")
