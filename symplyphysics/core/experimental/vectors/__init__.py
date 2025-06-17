# pylint: disable=too-many-lines

from __future__ import annotations

from typing import Any, Optional, Sequence, Self, Callable, TypeVar, Iterable
from collections import defaultdict
from itertools import product as py_product

from sympy import (Basic, Expr, S, Mul as SymMul, Add as SymAdd, Derivative as SymDerivative,
    fraction, sympify as sym_sympify, IndexedBase, Symbol as SymSymbol, Idx)
from sympy.core import function as sym_fn
from sympy.core.parameters import global_parameters
from sympy.combinatorics.permutations import Permutation
from sympy.physics.units import Dimension
from sympy.printing.printer import Printer
from sympy.tensor.indexed import Indexed

from symplyphysics.core.operations.sum_indexed import IndexedSum
from symplyphysics.core.symbols.symbols import (DimensionSymbol, next_name, Symbol,
    process_subscript_and_names, global_index)
from symplyphysics.core.symbols.id_generator import last_id, next_id
from ..miscellaneous import cacheit, sympify_expr

_T = TypeVar("_T")


def sort_with_sign(
    it: Iterable[_T],
    key: Optional[Callable[[_T], Any]] = None,
) -> tuple[int, list[_T]]:
    """
    Sorts `it` with an optional `key`. Also returns the sign of the permutation of the elements of
    `it` as the result of sorting, which is `1` if the number of swaps performed between any two
    elements is even, `-1` if the number is odd, and `0` if `it` contains any equal elements.
    """

    old_it = it if isinstance(it, list) else list(it)

    if key is None:
        new_it = sorted(old_it)
        indices = [old_it.index(v) for v in new_it]
    else:
        old_keys = [key(v) for v in old_it]
        new_keys = sorted(old_keys)
        indices = [old_keys.index(k) for k in new_keys]
        new_it = [old_it[i] for i in indices]

    if len(set(indices)) != len(indices):
        sign = 0
    else:
        sign = Permutation(indices).signature()

    return sign, new_it


def _check_vector(value: Any) -> Expr:
    """
    Raises error if `value` is not a vector expression. Converts the input into an `Expr` if
    needed.
    """

    if not is_vector_expr(value):
        raise ValueError(f"Expected '{value}' to be a vector.")

    return sympify_expr(value)


@cacheit
def is_vector_expr(value: Any) -> bool:  # pylint: disable=too-many-return-statements
    """
    Checks if `value` is a vector expression, for which `value` must either be `0` or a linear
    combination of other vector expressions, i.e. it can at most be a sum of vectors multiplied
    by scalars.
    """

    if value == 0:
        return True

    if isinstance(value, VectorExpr):
        return True

    if isinstance(value, SymAdd):
        return all(is_vector_expr(arg) for arg in value.args)

    if isinstance(value, SymMul):
        _, denominator = fraction(value)

        if is_vector_expr(denominator):
            return False

        n_vectors = 0
        has_zero = False

        for arg in value.args:
            if arg == 0:
                has_zero = True
                continue

            if is_vector_expr(arg):
                n_vectors += 1
                continue

        match n_vectors:
            case 0:
                return has_zero
            case 1:
                return True
            case _:
                raise ValueError("A vector can only be multiplied by a scalar.")

    if isinstance(value, SymDerivative):
        func, *args = value.args

        return is_vector_expr(func) and not any(is_vector_expr(arg) for (arg, _) in args)

    if isinstance(value, Indexed):
        return is_vector_expr(value.base)

    if isinstance(value, IndexedSum):
        return is_vector_expr(value.args[0])

    return False


def is_atomic_vector(value: Any) -> bool:
    """
    Checks if `value` represents an atomic vector, in the sense that it is irreducible and does
    not represent any vector operation, such as `VectorCross`.
    """

    return getattr(value, "is_atomic_vector", None) or isinstance(value,
        (VectorSymbol, AppliedVectorFunction))


@cacheit
def split_factor(value: Any) -> tuple[Expr, Expr]:
    """
    If `value` is a product of a vector and a scalar, returns the vector and the scalar. Otherwise
    returns `value` and `1`.

    Raises an error if `value` isn't a vector expression.
    """

    expr = _check_vector(value)

    if not isinstance(expr, SymMul):
        return expr, S.One

    for i, arg in enumerate(expr.args):
        if is_vector_expr(arg):
            factor = SymMul(*(arg for j, arg in enumerate(expr.args) if j != i))
            return arg, factor

    raise ValueError("The execution should not get here.")


@cacheit
def into_terms(value: Any) -> tuple[Expr, ...]:
    """
    Returns a tuple consisting of the terms which, when when added up, would give `value` again.

    Raises an error if `value` isn't a vector expression.
    """

    expr = _check_vector(value).expand()

    if expr == 0:
        return ()

    if isinstance(expr, SymAdd):
        return expr.args

    return (expr,)


def _ordered_mul(
    *values: Expr,
    key: Optional[Callable[[Expr], Any]] = None,
) -> dict[int, dict[tuple[Expr, ...], Expr]]:
    """
    Multiplies the vectors in `values` together ensuring that the order of multiplication is
    preserved. Returns the mapping of the permutation sign to the mapping of the (sorted) vector
    tuple to its factor.
    
    In other words, if given two vectors `v = sum(a_i * v_i, i)` and `w = sum(b_j * w_j, j)`, the
    result of their multiplication is `v * w = sum(a_i * b_j * v_i * w_j, i, j)`. The
    aforementioned tuple is the pair `(v_i, w_j)`, which is sorted and the sign of the permutation
    obtained as a result of sorting is taken into account. If `v_i = w_j`, the sign is `0`. The
    same considerations can be applied for more than two vectors.

    Raises an error if any value in `values` isn't a vector expression.
    """

    exprs = tuple(_check_vector(arg) for arg in values)

    if len(values) < 2:
        raise ValueError(f"Expected to multiply at least two vectors, got {len(values)}")

    key = key or id

    mapping: dict[int, dict[tuple[Expr, ...],
        Expr]] = defaultdict(lambda: defaultdict(lambda: S.Zero))

    for terms in py_product(*(into_terms(expr) for expr in exprs)):
        vectors: list[Expr] = []
        factors: list[Expr] = []
        has_zero = False

        for term in terms:
            if term == 0:
                has_zero = True
                break

            vector, factor = split_factor(term)
            vectors.append(vector)
            factors.append(factor)

        if has_zero:
            continue

        sign, sorted_vectors = sort_with_sign(vectors, key=key)
        factor = SymMul(*factors)

        mapping[sign][tuple(sorted_vectors)] += factor

    # Remove tuples corresponding to a factor of zero
    for sign, tuple_to_factor in mapping.items():
        mapping[sign] = {
            vectors: factor for vectors, factor in tuple_to_factor.items() if factor != 0
        }

    return mapping


class VectorExpr(Expr):  # pylint: disable=too-few-public-methods
    """
    Base class for vector expressions, excluding their sum and scalar multiplication, which are
    covered by `sympy.Add` and `sympy.Mul` respectively.
    """

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _eval_vector_norm(self) -> Optional[Expr]:
        return None

    @classmethod
    def _eval_vector_dot(cls, _lhs: Expr, _rhs: Expr) -> Optional[Expr]:
        return None

    @classmethod
    def _eval_vector_cross(cls, _lhs: Expr, _rhs: Expr) -> Optional[Expr]:
        return None


def _process_vector_names(
    code: Optional[str] = None,
    latex: Optional[str] = None,
    *,
    base: str = "VEC",
    i: Optional[int] = None,
) -> tuple[str, str]:
    if not code:
        code = f"{base}{i}"
        if not latex:
            base = base[0].lower()
            latex = f"{{\\vec {base}}}_{{{i}}}"
    elif not latex:
        latex = f"{{\\vec {code}}}"

    return code, latex


# TODO: Add support for axial vectors.
class VectorSymbol(Symbol, VectorExpr):  # pylint: disable=too-many-ancestors
    """
    Class representing a symbolic vector.

    The ``display_symbol``, ``display_latex``, and ``dimension`` parameters are used to instantiate
    the ``DimensionSymbol`` class.
    """

    def __new__(
        cls,
        display_symbol: Optional[str] = None,
        dimension: Dimension = Dimension(1),
        *,
        display_latex: Optional[str] = None,
    ) -> VectorSymbol:
        return super().__new__(cls)  # pylint: disable=no-value-for-parameter

    def __init__(
            self,
            display_symbol: Optional[str] = None,
            dimension: Dimension = Dimension(1),
            *,
            display_latex: Optional[str] = None,
    ) -> None:
        id_ = next_id("VEC")

        display_symbol, display_latex = _process_vector_names(display_symbol, display_latex, i=id_)

        Symbol.__init__(
            self,
            display_symbol=display_symbol,
            dimension=dimension,
            display_latex=display_latex,
        )
        VectorExpr.__init__(self)

    def _hashable_content(self) -> tuple[Any, ...]:
        return (id(self),)

    def _eval_derivative(self, s: Expr) -> Expr:
        # Needed for `sympy.solve` to work
        if is_vector_expr(s):
            return SymDerivative(self, s, evaluate=False)

        return S.Zero


class IndexedVectorSymbol(DimensionSymbol, VectorExpr, IndexedBase):
    index: Idx

    def __new__(
        cls,
        name_or_symbol: Optional[str | SymSymbol] = None,
        index: Optional[Idx] = None,
        dimension: Dimension = Dimension(S.One),
        *,
        display_latex: Optional[str] = None,
    ) -> IndexedVectorSymbol:
        # SymPy subs() and solve() creates dummy symbols. Allow create new indexed symbols
        # without renaming
        if isinstance(name_or_symbol, SymSymbol):
            obj = IndexedBase.__new__(cls, name_or_symbol)
        else:
            obj = IndexedBase.__new__(cls, next_name("VEC"))

        return obj

    def __init__(
            self,
            name_or_symbol: Optional[str | SymSymbol] = None,
            index: Optional[Idx] = None,
            dimension: Dimension = Dimension(S.One),
            *,
            display_latex: Optional[str] = None,
    ) -> None:
        if name_or_symbol is None:
            display_name = str(self.name)
        else:
            display_name = str(name_or_symbol)

        self.index = index or global_index

        DimensionSymbol.__init__(
            self,
            display_name=display_name,
            dimension=dimension,
            display_latex=display_latex,
        )
        VectorExpr.__init__(self)
        IndexedBase.__init__(self)


class VectorNorm(Expr):  # pylint: disable=too-few-public-methods
    """
    Class representing the Euclidean norm (see link 1, *Euclidean norm*) of a vector expression.

    The vector argument is stored in position `0` of `self.args`.

    The vector norm has the following **properties**:

    1. **Subadditivity**: for all vectors `a` and `b`, `norm(a + b) <= norm(a) + norm(b)`.

    2. **Absolute homogeneity**: for all scalars `k` and vectors `a`, `norm(k * a) = abs(k) * norm(a)`.

    3. **Positive definiteness**: for all vectors `a`, `norm(a) = 0` if and only if `a = 0`.

    **Links:**

    1. `Wikipedia <https://en.wikipedia.org/wiki/Norm_(mathematics)>`__.
    """

    is_nonnegative = True
    is_negative = False
    is_real = True
    is_commutative = True

    @cacheit
    def __new__(cls, value: Any, *, evaluate: Optional[bool] = None) -> Expr:
        vector = _check_vector(value)

        if evaluate is None:
            evaluate = global_parameters.evaluate

        if not evaluate:
            obj = super().__new__(cls)  # pylint: disable=no-value-for-parameter
            obj._args = (vector,)
            return obj

        if vector == 0:
            return S.Zero

        vector = vector.doit()

        if isinstance(vector, VectorExpr):
            norm = vector._eval_vector_norm()

            if norm is not None:
                return norm

        # Assuming `a` is a scalar and `v, w` are vectors, converts
        # `a * v + a * w` to `a * (v + w)` if possible
        if isinstance(vector, SymAdd):
            vector = vector.together()

        # The above simplification isn't possible
        if not isinstance(vector, SymMul):
            return cls(vector, evaluate=False)

        # Also handles the case when `vector` is simply another vector multiplied by a scalar
        vector, factor = split_factor(vector)

        return cls(vector, evaluate=False) * abs(factor)

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        return f"norm({p.doprint(self.args[0])})"

    def _eval_derivative(self, symbol: Expr) -> Expr:
        # Needed for `sympy.solve` to work
        if is_vector_expr(symbol):
            return SymDerivative(self, symbol, evaluate=False)

        done = self.doit()

        if not isinstance(done, VectorNorm):
            return done.diff(symbol)

        # (d/dx)[norm(v(x))]
        # = (d/dx)[sqrt(dot(v(x), v(x)))]
        # = (d/dx)[dot(v(x), v(x))] / (2 * sqrt(v(x), v(x)))
        # = (2 * dot(v(x), (d/dx)[v(x)])) / (2 * norm(v(x)))
        # = dot(v(x), (d/dx)[v(x)]) / norm(v(x))

        vector: Expr = done.args[0]
        return VectorDot(vector, vector.diff(symbol)) / done


class VectorDot(Expr):
    """
    The **dot product**, or **scalar product**, is a binary operation that takes two vectors and
    returns a single number.

    Geometrically, the dot product can be expressed using the length (`norm`) of the vectors and
    the (non-directional) `angle` between them: `dot(a, b) = norm(a) * norm(b) * cos(angle(a, b))`
    where `a` and `b` are vectors.

    In particular,

    1. `dot(a, b) = 0` if and only if `a` and `b` are orthogonal.

    2. `dot(a, b) = norm(a) * norm(b)` if and only if `a` and `b` are codirectional.

    3. `dot(a, a) = norm(a)^2` as a result of (2).

    The dot product has the following **properties**:

    1. **Commutativity**: for all vectors `a, b`, `dot(a, b) = dot(b, a)`.

    2. **Linearity** in the **first** argument: for all vectors `a, b, c` and scalars `k, l`,
       `dot(a * k + b * l, c) = k * dot(a, c) + l * dot(b, c)`.

    3. **Linearity** in the **second** argument, which follows from (1) and (2).

    4. **Absense of cancellation**: for all vectors `a, b, c` s.t. `a ≠ 0`, `dot(a, b) = dot(a, c)`
       does not imply `b = c`.

    5. Applicability of the **product rule**: for all vector-valued differentiable functions
       `a, b`, `d[dot(a, b)] = dot(d[a], b) + dot(a, d[b])` where `d[v]` denotes the derivative of
       vector `v`.

    Note that the properties and relations mentioned only apply to **real-valued vectors**.

    The dot product is a *true scalar* in a sense that it is unchanged if the orientation of the
    frame is reversed.
    """

    @property
    def lhs(self) -> VectorExpr:
        return self.args[0]

    @property
    def rhs(self) -> VectorExpr:
        return self.args[1]

    @cacheit
    # pylint: disable-next=too-many-branches
    def __new__(cls, *values: Any, evaluate: Optional[bool] = None) -> Expr:
        if len(values) != 2:
            raise ValueError(f"{cls.__name__} accepts two positional arguments.")

        lhs, rhs = map(_check_vector, values)

        if evaluate is None:
            evaluate = global_parameters.evaluate

        if not evaluate:
            obj = super().__new__(cls)  # pylint: disable=no-value-for-parameter
            obj._args = (lhs, rhs)
            return obj

        lhs, rhs = lhs.doit(), rhs.doit()

        if isinstance(lhs, VectorExpr):
            result = lhs._eval_vector_dot(lhs, rhs)

            if result is not None:
                return result

        if isinstance(rhs, VectorExpr):
            result = rhs._eval_vector_dot(lhs, rhs)

            if result is not None:
                return result

        sign_to_mapping = _ordered_mul(lhs, rhs)
        result = S.Zero

        for sign, tuple_to_factor in sign_to_mapping.items():
            # Vector is multiplied by itself, reduce product to norm squared
            if sign == 0:
                for (v, _), factor in tuple_to_factor.items():
                    dot = VectorNorm(v)**2
                    result += dot * factor

                continue

            # The dot product is commutative, i.e. the sign of the permutation of the symbols doesn't
            # matter
            for (v, w), factor in tuple_to_factor.items():
                if is_atomic_vector(v) and is_atomic_vector(w):
                    dot = cls(v, w, evaluate=False)
                else:
                    dot = cls(v, w)

                result += dot * factor

        return result

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        return f"dot({p.doprint(self.lhs)}, {p.doprint(self.rhs)})"

    def _eval_derivative(self, symbol: Expr) -> Expr:
        # Needed for `sympy.solve` to work
        if is_vector_expr(symbol):
            return SymDerivative(self, symbol, evaluate=False)

        lhs = self.lhs
        rhs = self.rhs

        derived_lhs = VectorDot(lhs.diff(symbol), rhs)
        derived_rhs = VectorDot(lhs, rhs.diff(symbol))

        return derived_lhs + derived_rhs


class VectorCross(VectorExpr):  # pylint: disable=too-few-public-methods
    """
    The **cross product**, or **vector product**, is a binary operation that takes two vectors and
    returns another vector. The cross product is only defined in a *3-dimensional space* (although
    its construction is also possible in a 7-dimensional space, the following properties do not
    hold there).

    Geometrically, the cross product between vectors `a` and `b` can be defined as `cross(a, b) =
    norm(a) * norm(b) * sin(angle(a, b)) * n` where `norm` is the length operator and `n` is a unit
    vector orthogonal to the plane containing `a` and `b` with such direction that the ordered set
    `(a, b, n)` is positively oriented.

    The cross product has the following **properties**:

    1. **Anticommutativity**: for all vectors `a, b`, `cross(a, b) = -cross(b, a)`.

    2. For any vector `a`, `cross(a, a) = 0`, which follows from (1).

    3. **Distributivity over addition**: for all vectors `a, b, c`, `cross(a, b + c) =
       cross(a, b) + cross(a, c)`.

    4. **Absense of associativity**: for all vectors `a, b, c`, `cross(a, cross(b, c)) ≠
       cross(cross(a, b), c)`.

    5. However, the **Jacobi identity** is satisfied: for all vectors `a, b, c`,
       `cross(a, cross(b, c)) + cross(b, cross(c, a)) + cross(c, cross(b, a)) = 0`.

    6. **Absense of cancellation**: for all vectors `a, b, c` s.t. `a ≠ 0`, `cross(a, b) =
       cross(a, c)` does not imply `b = c`. This only happens if `dot(a, b) = dot(a, c)` holds.

    7. Applicability of the **product rule**: for all vector-valued differentiable functions
       `a, b`, `d[cross(a, b)] = cross(d[a], b) + cross(a, d[b])` where `d[v]` denotes the
       derivative of vector `v`.

    It is related to the dot product by the following relation: for all vectors `a, b`,
    `norm(cross(a, b))^2 + dot(a, b)^2 = norm(a)^2 * norm(b)^2`.

    The cross product is a *pseudovector*, i.e. it is negated if the orientation of the frame is
    reversed.
    """

    @cacheit
    def __new__(cls, *values: Any, evaluate: Optional[bool] = None) -> VectorExpr:
        if len(values) != 2:
            raise ValueError(f"{cls.__name__} accepts two positional arguments.")

        lhs, rhs = map(_check_vector, values)

        if evaluate is None:
            evaluate = global_parameters.evaluate

        if not evaluate:
            obj = super().__new__(cls)  # pylint: disable=no-value-for-parameter
            obj._args = (lhs, rhs)
            return obj

        lhs, rhs = lhs.doit(), rhs.doit()

        if isinstance(lhs, VectorExpr):
            result = lhs._eval_vector_cross(lhs, rhs)

            if result is not None:
                return result

        if isinstance(rhs, VectorExpr):
            result = rhs._eval_vector_cross(lhs, rhs)

            if result is not None:
                return result

        result = S.Zero
        sign_to_mapping = _ordered_mul(lhs, rhs)

        for sign, tuple_to_factor in sign_to_mapping.items():
            # Cross product is zero when arguments are equal
            if sign == 0:
                continue

            # Note that `vectors` is already sorted
            for vectors, factor in tuple_to_factor.items():
                v, w = vectors

                if is_atomic_vector(v) and is_atomic_vector(w):
                    cross = cls(v, w, evaluate=False)
                else:
                    cross = cls(v, w)

                result += cross * factor * sign

        return result

    def _sympystr(self, p: Printer) -> str:
        lhs, rhs = self.args
        return f"cross({p.doprint(lhs)}, {p.doprint(rhs)})"

    @classmethod
    def _eval_vector_dot(cls, lhs: Expr, rhs: Expr) -> Optional[Expr]:
        lhs_is_cross = isinstance(lhs, VectorCross)
        rhs_is_cross = isinstance(rhs, VectorCross)

        if lhs_is_cross and rhs_is_cross:
            a, b = lhs.args
            c, d = rhs.args

            return VectorDot(a, b) * VectorDot(c, d) - VectorDot(b, c) * VectorDot(a, d)

        if lhs_is_cross and not rhs_is_cross:
            a, b = lhs.args

            return VectorMixedProduct(rhs, a, b)

        if not lhs_is_cross and rhs_is_cross:
            c, d = rhs.args

            return VectorMixedProduct(lhs, c, d)

        return None

    @classmethod
    def _eval_vector_cross(cls, lhs: Expr, rhs: Expr) -> Optional[Expr]:
        lhs_is_cross = isinstance(lhs, VectorCross)
        rhs_is_cross = isinstance(rhs, VectorCross)

        if lhs_is_cross and rhs_is_cross:
            a, b = lhs.args
            c, d = rhs.args

            return c * VectorMixedProduct(d, a, b) - d * VectorMixedProduct(c, a, b)

        if lhs_is_cross and not rhs_is_cross:
            a, b = lhs.args

            return b * VectorDot(rhs, a) - a * VectorDot(rhs, b)

        if not lhs_is_cross and rhs_is_cross:
            c, d = rhs.args

            return c * VectorDot(lhs, d) - d * VectorDot(lhs, c)

        return None

    def _eval_derivative(self, symbol: Expr) -> Expr:
        # Needed for `sympy.solve` to work
        if is_vector_expr(symbol):
            return SymDerivative(self, symbol, evaluate=False)

        lhs, rhs = self.args

        derived_lhs = VectorCross(lhs.diff(symbol), rhs)
        derived_rhs = VectorCross(lhs, rhs.diff(symbol))

        return derived_lhs + derived_rhs


class VectorMixedProduct(Expr):  # pylint: disable=too-few-public-methods
    """
    The **scalar triple product**, or **mixed product**, is defined as the dot product of one
    vector with the cross product with the other two.

    Geometrically, given vectors `a, b, c` the scalar triple product `dot(a, cross(b, c))` can be
    interpreted as the signed volume of the parallelepiped defined by these vectors.

    The mixed product has the following **properties**:

    1. It does not change under a positive permutation of the arguments: for all vectors `a, b, c`,
       `mixed(a, b, c) = mixed(b, c, a) = mixed(c, a, b)`.

    2. It is negated under a negative permutation of the arguments: for all vectors `a, b, c`,
       `mixed(a, c, b) = mixed(b, a, c) = mixed(c, b, a) = -mixed(a, b, c)`.

    3. The scalar triple product is zero if and only if the three vectors `a, b, c` are coplanar.

    4. It is linear in all arguments.

    The scalar triple product is a *pseudoscalar*, i.e. it is negated if the orientation of the
    frame is reversed.
    """

    is_real = True

    @cacheit
    def __new__(cls, *values: Expr, evaluate: Optional[bool] = None) -> Expr:
        if len(values) != 3:
            raise ValueError(f"{cls.__name__} requires three positional arguments.")

        a, b, c = map(_check_vector, values)

        if evaluate is None:
            evaluate = global_parameters.evaluate

        if not evaluate:
            obj = super().__new__(cls)  # pylint: disable=no-value-for-parameter
            obj._args = (a, b, c)
            return obj

        a, b, c = a.doit(), b.doit(), c.doit()

        sign_to_mapping = _ordered_mul(a, b, c)
        result = S.Zero

        for sign, tuple_to_factor in sign_to_mapping.items():
            # Mixed product is zero if some arguments are equal to one another
            if sign == 0:
                continue

            for vectors, factor in tuple_to_factor.items():
                if all(is_atomic_vector(vector) for vector in vectors):
                    mixed = cls(*vectors, evaluate=False)
                else:
                    u, v, w = vectors

                    mixed = VectorDot(u, VectorCross(v, w))

                result += mixed * factor * sign

        return result

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        a, b, c = self.args
        return f"mixed({p.doprint(a)}, {p.doprint(b)}, {p.doprint(c)})"

    def _eval_derivative(self, symbol: Expr) -> Expr:
        # Needed for `sympy.solve` to work
        if is_vector_expr(symbol):
            return SymDerivative(self, symbol, evaluate=False)

        a, b, c = self.args
        return VectorDot(a, VectorCross(b, c)).diff(symbol)


class VectorJacobian(VectorExpr):

    def __new__(
        cls,
        function: VectorFunction,
        index: int,
        argument: Expr,
        /,
    ) -> VectorJacobian:
        return super().__new__(cls, function, sympify_expr(index), sympify_expr(argument))


class AppliedVectorFunction(sym_fn.Application, VectorExpr):
    """This class represents the result of applying a vector-valued function to some arguments."""

    @cacheit
    def __new__(cls, *args: Any, **kwargs: Any) -> Self:
        if cls is AppliedVectorFunction:
            raise TypeError("Call `VectorFunction` instead to instantiate a new function.")

        n = len(args)

        if not cls._valid_nargs(n):
            template = "{name} takes {qual} {args} argument{plural} ({given} given)"
            arguments = getattr(cls, "arguments", None)
            nargs = min(cls.nargs) if arguments is None else len(arguments)
            message = template.format(
                name=cls,
                qual="exactly" if len(cls.nargs) == 1 else "at least",
                args=nargs,
                plural="s" * (nargs != 1),
                given=n,
            )
            raise TypeError(message)

        args = tuple(sym_sympify(arg, strict=True) for arg in args)
        undefineds = [arg.name for arg in args if isinstance(arg, sym_fn.FunctionClass)]
        if undefineds:
            template = "Invalid argument: expecting an expression, not undefined function{plural}: {names}"
            message = template.format(
                plural="s" * (len(undefineds) > 1),
                names=", ".join(undefineds),
            )
            raise TypeError(message)

        result = super().__new__(cls, *args, **kwargs)
        return result

    def _eval_derivative(self, symbol: Expr) -> Expr:
        # For now, the derivative of a function application is only implemented when no more than
        # one argument of the function is equal to `symbol`. All other cases require a vector
        # alternative of `sympy.Subs` or a Jacobian matrix.

        # Needed for `sympy.solve` to work
        if is_vector_expr(symbol):
            return SymDerivative(self, symbol, evaluate=False)

        derived_args = []

        found = False

        for i, arg in enumerate(self.args):
            if is_vector_expr(arg) and arg != 0:
                # Assuming `v` and `w` are vector-valued functions:
                #   (d/dx)[v(w(x))] = jacobian[v](w(x)) * (d/dx)[w(x)]
                # Here, `jacobian[v](w(x))` denotes the Jacobian matrix of `v` evaluated at `w(x)`.

                # TODO: finish implementing the `VectorJacobian` class logic
                return VectorJacobian(self, i, arg)

            scale = arg.diff(symbol)

            if scale == 0:
                continue

            if scale != 1 or found:
                raise NotImplementedError("VectorSubs has not been implemented yet.")

            found = True

            partial = VectorDerivative(self, arg)
            derived_args.append(partial * scale)

        if not found:
            return S.Zero

        return SymAdd(*derived_args)


class UndefinedVectorFunction(sym_fn.FunctionClass):
    """The (meta)class of undefined vector functions."""

    def __new__(
        mcs,
        name: str,
        bases: Optional[Sequence[type]] = None,
        __dict__: Optional[dict[str, Any]] = None,
        **kwargs: Any,
    ) -> Self:
        bases = bases or (AppliedVectorFunction,)

        if __dict__ is None:
            __dict__ = {}
        __dict__ |= kwargs
        __dict__["_kwargs"] = kwargs
        __dict__["__module__"] = None

        obj = super().__new__(mcs, name, bases, __dict__)
        obj.name = name
        return obj


class VectorFunction(DimensionSymbol, UndefinedVectorFunction):
    _arguments: Optional[tuple[Basic]]

    @property
    def arguments(cls) -> Optional[tuple[Basic]]:
        return cls._arguments

    def __new__(  # pylint: disable=signature-differs
        mcs,
        display_name: Optional[str] = None,
        arguments: Optional[Sequence[Basic]] = None,
        *,
        dimension: Dimension = Dimension(1),
        display_latex: Optional[str] = None,
        **kwargs: Any,
    ) -> Self:
        name = next_name("FUN")
        obj = UndefinedVectorFunction.__new__(mcs, name, **kwargs)

        return obj

    def __init__(
        cls,
        display_name: Optional[str] = None,
        arguments: Optional[Sequence[Basic]] = None,
        *,
        dimension: Dimension = Dimension(1),
        display_latex: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        if arguments is not None and not isinstance(arguments, tuple):
            arguments = tuple(arguments)
        cls._arguments = arguments

        display_name, display_latex = _process_vector_names(
            display_name,
            display_latex,
            base="FUN",
            i=last_id("FUN"),
        )
        DimensionSymbol.__init__(
            cls,
            display_name=display_name,
            dimension=dimension,
            display_latex=display_latex,
        )

        if arguments is not None:
            kwargs["nargs"] = len(arguments)
        UndefinedVectorFunction.__init__(cls, **kwargs)

    def __repr__(cls) -> str:
        return str(cls.display_name)


class VectorDerivative(SymDerivative, VectorExpr):  # pylint: disable=too-few-public-methods

    def __new__(cls, expr: Expr, *variables: Expr, **kwargs: Any) -> Expr:
        expr = _check_vector(expr)

        kwargs.setdefault("evaluate", False)
        derivative = super().__new__(cls, expr, *variables, **kwargs)

        if not isinstance(derivative, SymDerivative):
            return derivative

        for symbol, _ in derivative.args[1:]:
            if is_vector_expr(symbol):
                raise ValueError(f"Unable to differentiate with respect to vector: {symbol}")

        return derivative


def vector_diff(expr: Expr, *variables: Expr, **kwargs: Any) -> Expr:
    """Evaluate the derivative of `expr` with respect to `variables`."""

    return VectorDerivative(
        expr,
        *(sympify_expr(variable) for variable in variables),
        **kwargs,
    ).doit()


def clone_as_vector_symbol(
    source: DimensionSymbol,
    *,
    display_symbol: Optional[str] = None,
    display_latex: Optional[str] = None,
    subscript: Optional[str] = None,
) -> VectorSymbol:
    name = display_symbol or source.display_name

    # NOTE: assuming the source latex code allows to be placed after "\vec"
    latex = display_latex or f"{{\\vec {source.display_latex}}}"

    name, latex = process_subscript_and_names(name, latex, subscript)

    return VectorSymbol(
        display_symbol=name,
        dimension=source.dimension,
        display_latex=latex,
    )


def clone_as_indexed_vector(
    source: DimensionSymbol,
    index: Optional[Idx] = None,
    *,
    display_symbol: Optional[str] = None,
    display_latex: Optional[str] = None,
) -> IndexedVectorSymbol:
    name = display_symbol or source.display_name

    # NOTE: assuming the source latex code allows to be placed after "\vec"
    latex = display_latex or f"{{\\vec {source.display_latex}}}"

    return IndexedVectorSymbol(
        name_or_symbol=name,
        index=index,
        dimension=source.dimension,
        display_latex=latex,
    )


def clone_as_vector_function(
    source: DimensionSymbol,
    arguments: Optional[Sequence[Any]] = None,
    *,
    display_symbol: Optional[str] = None,
    display_latex: Optional[str] = None,
    subscript: Optional[str] = None,
    **kwargs: Any,
) -> VectorFunction:
    name = display_symbol or source.display_name
    latex = display_latex or source.display_latex

    name, latex = process_subscript_and_names(name, latex, subscript)

    return VectorFunction(
        display_name=name,
        arguments=arguments,
        dimension=source.dimension,
        display_latex=latex,
        **kwargs,
    )


def split_into_tangential_and_normal_components(original: Expr, target: Expr) -> tuple[Expr, Expr]:
    """
    Splits `original` into the component tangential to `target` and the component normal to
    `target`. Note that when added together, they give back `original`.
    """

    tangential_component = VectorDot(original, target) / VectorDot(target, target) * target

    normal_component = original - tangential_component

    return tangential_component, normal_component


__all__ = [
    "AppliedVectorFunction",
    "VectorCross",
    "VectorDot",
    "VectorExpr",
    "VectorFunction",
    "VectorMixedProduct",
    "VectorNorm",
    "VectorSymbol",
    "VectorDerivative",
    "clone_as_vector_symbol",
    "clone_as_vector_function",
    "split_into_tangential_and_normal_components",
]
