from sympy import solve, Eq

from symplyphysics import Symbol, Quantity, validate_input, \
    validate_output, units
from symplyphysics.core.symbols.fraction import Fraction

from symplyphysics.laws.optics import linear_magnification_from_object_height_and_image_height as magnification_law


class BaseLaw:
    def __init__(self, law: Eq) -> None:
        self.law = law
        self.symbols = law.atoms(Symbol)

    def calculate_symbol_value(self, variable: Symbol, dict_quantities: dict[Symbol, Quantity | Fraction]
    ):
        if variable in self.symbols:
            pass
        else:
            raise ValueError("The specified variable is not in the law")

        symbols_set = self.symbols.copy()
        symbols_set.remove(variable)

        variables_list = list(symbols_set)

        variables_for_calculate = list(str(symbol) for symbol in variables_list)
        dimensions = list(symbol.dimension for symbol in variables_list)

        dict_validate = dict(
            zip(variables_for_calculate, dimensions)
        )

        @validate_input(**dict_validate)
        @validate_output(variable.dimension)
        def __calculate_variable(variable_: Symbol, *quantities_: Quantity | Fraction) -> Quantity | Fraction:
            solved = solve(self.law, variable_, dict=True)[0][variable]

            result_expr = solved.subs(zip(variables_list, quantities_))

            result = Quantity(result_expr)
            return result

        return __calculate_variable(variable, *dict_quantities.values())


law = BaseLaw(magnification_law.law)
print(magnification_law.object_height.dimension)
print(
    law.calculate_symbol_value(
        variable=magnification_law.magnification,
        dict_quantities={
            magnification_law.object_height: Quantity(5 * units.meters),
            magnification_law.image_height: Quantity(7 * units.meters)
        }
    )
)


try:
    law.calculate_symbol_value(
        variable=magnification_law.magnification,
        dict_quantities={
            magnification_law.object_height: Quantity(5 * units.meters),
            magnification_law.image_height: Quantity(7 * units.coulomb)
        }
    )
except Exception as ex:
    print(ex)


# BUT not exceptions
try:
    print(
        law.calculate_symbol_value(
            variable=magnification_law.magnification,
            dict_quantities={
                magnification_law.object_height: Quantity(5 * units.meters),
                magnification_law.image_height: Quantity(7 * units.length)
            }
        )
    )
except Exception as ex:
    print(ex)
