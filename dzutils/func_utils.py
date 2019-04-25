import functools as _functools


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
