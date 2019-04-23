import functools as _functools


def compose(*functions):
    """
    Takes any number of functions, returns the composed function
    """

    def compose2(f, g):
        return lambda x: f(g(x))

    return _functools.reduce(compose2, functions, lambda x: x)

def boil_fuction (func,*args,**kwargs):
    """
    Returns the output of func(*args,**kwargs), on exception returns None

    Jenky utility for working with ambiguous inputs
    """
    try:
        return func(*args, **kwargs)
    except:
        return None
def boil_generator (func,*args,**kwargs):
    """
    Returns the output of func(*args,**kwargs), on exception returns empty iterable

    Jenky utility for working with ambiguous inputs
    """
    try:
        return func(*args,**kwargs)
    except:
        return iter(())
