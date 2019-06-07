import functools as _functools
import numpy as _np
import os as _os
import sys as _sys


def compose(*functions):
    """
    Takes any number of functions, returns the composed function
    """

    def compose2(f, g):
        return lambda x: f(g(x))

    return _functools.reduce(compose2, functions, lambda x: x)


def boil_fuction(func, *args, **kwargs):
    """
    Returns the output of func(*args,**kwargs), on exception returns None

    Jenky utility for working with ambiguous inputs
    """
    try:
        return func(*args, **kwargs)
    except:
        return None


def boil_generator(func, *args, **kwargs):
    """
    Returns the output of func(*args,**kwargs), on exception returns empty iterable

    Jenky utility for working with ambiguous inputs
    """
    try:
        return func(*args, **kwargs)
    except:
        return iter(())


def index_to_parameters(index, *args):
    """
    takes tuples where (minimum,maximum,samples) and index, returns all params

    Samples should be an integer value!
    Don't take half a sample and put the rest back.
    It's okay if the type is not int.
    """
    output = [
        arg[0]
        + ((arg[1] - arg[0] + 1) // (arg[2]))
        * ((index // _np.prod([a[2] for a in args[:i]])) % arg[2])
        for i, arg in enumerate(args)
    ]
    return output


def index_to_pair(index, sample, radius):
    """
    Takes the index and the distance within that index not to pair with the index in the output,
    returns a pair of values in the ranges allowed by sample and radius

    Sample should be a tuple  (minimum,maximum,samples)  radius should be an int limiting the sampling of the paired value with that in sample

    Radius of 0 means sample all pairs, radius at 1 means only use i +/- (1+) etc.
    This docstring is jenky and confusing, and if you're using this function, you're probably doing something jenky and confusing

    """
    relative = index_to_parameters(
        index,
        (sample[0], sample[1], sample[2]),
        (
            sample[0] + radius,
            sample[1] - radius,
            sample[2] - (radius != 0) - 2 * (radius),
        ),
    )
    absolute = [
        relative[0],
        (relative[0] + relative[1]) % num_residues
        + ((relative[0] + relative[1]) % num_residues == 0) * num_residues,
    ]
    return absolute


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = _sys.stdout
        _sys.stdout = open(_os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        _sys.stdout.close()
        _sys.stdout = self._original_stdout


######### The following stuff is yanked from Guido Van Rossum's blog ##########
_registry = {}


class MultiMethod(object):
    def __init__(self, name):
        self.name = name
        self.typemap = {}

    def __call__(self, *args):
        types = tuple(arg.__class__ for arg in args)  # a generator expression!
        function = self.typemap.get(types)
        if function is None:
            raise TypeError("no match")
        return function(*args)

    def register(self, types, function):
        if types in self.typemap:
            raise TypeError("duplicate registration")
        self.typemap[types] = function


# Guido Notes that this is not thread safe, and also the risk of namespace
# collision of __lastreg__
def multimethod(*types):
    def register(function):
        function = getattr(function, "__lastreg__", function)
        name = function.__name__
        mm = _registry.get(name)
        if mm is None:
            mm = _registry[name] = MultiMethod(name)
        mm.register(types, function)
        mm.__lastreg__ = function
        return mm

    return register


######### End Guido's code ####################################################
