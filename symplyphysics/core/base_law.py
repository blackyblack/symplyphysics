from sympy import solve, Eq

from symplyphysics import Symbol, Quantity, validate_input, validate_output
from symplyphysics.core.symbols.fraction import Fraction


class BaseLaw:
    def __init__(self, law: Eq) -> None:
        self.law = law
        self.symbols = law.atoms(Symbol)

    def calculate_symbol_value(self, variable: Symbol, dict_quantities: dict[Symbol, Quantity | Fraction]
    ) -> Quantity | Fraction:
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
