from sympy.vector import vector_integrate

# Description
## Field circulation along closed curve is curvilinear integral of field by this curve.
## C = CurveIntegral(F * dl, Curve), where
## C is circulation
## F is field vector function
## l if radius-vector of the curve

#TODO: implement and use unevaluated form of vector_integrate

# field_ should be instance of CoordSys3D
# trajectory_ is usually a ParametricRegion, eg ParametricRegion((cos(theta), sin(theta)), (theta, 0, pi/2))
def calculate_circulation(field_, trajectory_):
    return vector_integrate(field_, trajectory_)
