from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
#If there is no external force is applied to system of objects, the summary momentum of this system remains constant during and after any interactions between objects
#Scalar version for two colliding objects in 1-dimentional env. m1, V1, m2, V2 are masses and Velocities befor collision, m3, V3, m4 and V4 - after.
#Also applicable for reactive engine simulation

mass1, mass2, velocity1, velocity2, mass3, velocity3, mass4, velocity4 = symbols('mass1 velocity1 mass2 velocity2 mass3 velocity3 mass4 velocity4')
law = Eq(mass3 * velocity3 + mass4 * velocity4, mass1 * velocity1 + mass2 * velocity2)

def print():
    return pretty(law, use_unicode=False)

@validate_input(mass1_=units.mass, velocity1_=units.velocity, mass2_=units.mass, valocity2_=units.velocity, mass3_=units.mass, velocity3_=units.velocity, mass4_=units.mass, velocity4_=units.velocity)
@validate_output(units.mass)
@validate_output(units.velocity)

def calculate_result_velocity(mass1_: Quantity, velocity1_: Quantity, mass2_ : Quantity, velocity2_ : Quantity, mass3_ : Quantity, velocity3_ : Quantity, mass4_ : Quantity) -> Quantity:
    result_velocity_expr = solve(law, mass1_, velocity1_, mass2_, velocity2_, mass3_, velocity3_, mass4_, dict=True)[0]
    result_expr = result_velocity_expr.subs({mass1: mass1_, velocity1: velocity1_, mass2: mass2_, velocity2: velocity2_, mass3: mass3_, velocity3: velocity3_, mass4: mass4_})
    return expr_to_quantity(result_expr, 'Velocity4')

def calculate_result_mass(mass1_: Quantity, velocity1_: Quantity, mass2_ : Quantity, velocity2_ : Quantity, mass3_ : Quantity, velocity3_ : Quantity, velocity4_ : Quantity) -> Quantity:
    result_velocity_expr = solve(law, mass1_, velocity1_, mass2_, velocity2_, mass3_, velocity3_, velocity4_, dict=True)[0]
    result_expr = result_velocity_expr.subs({mass1: mass1_, velocity1: velocity1_, mass2: mass2_, velocity2: velocity2_, mass3: mass3_, velocity3: velocity3_, velocity4: velocity4_})
    return expr_to_quantity(result_expr, 'mass4')