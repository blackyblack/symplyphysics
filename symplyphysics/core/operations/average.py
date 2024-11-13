from .single_argument_operation import SingleArgumentOperation

class Average(SingleArgumentOperation):  # pylint: disable=too-many-ancestors
    r"""
    The average value :math:`\langle A \rangle` of a quantity :math:`A` which
    is a function of a random variable :math:`t` (*time*) over the interval
    :math:`(t_0, t_1)` is defined as follows:

    :code:`avg(A) = Integral(A(t), (t, t_0, t_1)) / (t_1 - t_0)`

    Latex:
        .. math::
            \langle A \rangle = \frac{1}{t_1 - t_0} \int \limits_{t_0}^{t_1} A(t) dt

    The interval is usually assumed to be the biggest possible domain of the
    variable :math:`t`

    Another possible definition might be derived from a *statistical ensemble*
    featuring the probability distribution :math:`\rho` over the microstates:

    :code:`avg(A) = Integral(rho(p) * A(p), p)`

    Latex:
        .. math::
            \langle A \rangle = \int \rho A(p) dp

    where `p` is the canonical momentum and one degree of freedom is assumed
    to exist in the system. The full formula can be found `here
    <https://en.wikipedia.org/wiki/Ensemble_(mathematical_physics)#Classical_mechanical>`__.

    Under the `ergodic hypothesis <https://en.wikipedia.org/wiki/Ergodic_hypothesis>`__,
    both the time average and the ensemble average are the same.
    """
