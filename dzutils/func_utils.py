import functools as _functools


def compose(*functions):
    """
    Takes any number of functions, returns the composed function
    """

    def compose2(f, g):
        return lambda x: f(g(x))

    return _functools.reduce(compose2, functions, lambda x: x)
