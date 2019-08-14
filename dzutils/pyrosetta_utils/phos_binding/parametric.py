from pyrosetta.rosetta.protocols.helical_bundle import MakeBundle
from numpy import linspace
from itertools import product


class BundleGridRangeParam(object):
    """
    Param range class
    """

    def __init__(self, name, start, stop=0, steps=1):
        self.name_string = name
        self.lower_bound = start
        self.upper_bound = stop
        self.steps = steps

    def get_range(self):
        return linspace(self.lower_bound, self.upper_bound, num=self.steps)

    def to_dict(self):
        return {
            "name": self.name,
            "type": "range",
            "start": self.start,
            "stop": self.stop,
            "steps": self.steps,
        }


class BundleGridBinaryParam(BundleGridRangeParam):
    """
    Derived from param range except the range is just two or one value
    """

    def __init__(self, name, base_value, allow_other=False):
        self.name_string = name
        self.first_value = base_value
        self.second_value = bool(not base_value) if allow_other else base_value

    def get_range(self):
        return list(set((self.first_value, self.second_value)))

    def to_dict(self):
        return {
            "name": self.name,
            "type": "binary",
            "base_value": self.first_value,
            "allow_others": True
            if self.first_value != self.second_value
            else False,
        }


class HelixGridParam(object):
    def __init__(self, num, length):
        self.num = num
        self.length = length
        self._params = {}

    def set_range_param(self, name, start, stop=0, steps=1):
        """
        takes a BundleGridParam value and registers it with the params dict
        just giving start will yield a range containing a single value
        """
        self._params[name] = BundleGridRangeParam(name, start, stop, steps)

    def set_binary_param(self, name, base_value, allow_other=False):
        """
        takes a BundleGridParam value and registers it with the params dict

        Can either have a single value or two options
        """
        self._params[name] = BundleGridBinaryParam(
            name, base_value, allow_other
        )

    def set_r0(self, start, stop=0, steps=1):
        """
        Set the r0 param range

        Start required, default is no sampling (steps=1)
        """
        self.set_range_param("r0", start, stop, steps)

    def set_omega0(self, start, stop=0, steps=1):
        """
        Set the omega0 param range

        Start required, default is no sampling (steps=1)
        """
        self.set_range_param("omega0", start, stop, steps)

    def set_omega1(self, start, stop=0, steps=1):
        """
        Set the omega1 param range

        Start required, default is no sampling (steps=1)
        """
        self.set_range_param("omega1", start, stop, steps)

    def set_delta_omega0(self, start, stop=0, steps=1):
        """
        Set the delta_omega0 param range

        Start required, default is no sampling (steps=1)
        """
        self.set_range_param("delta_omega0", start, stop, steps)

    def set_delta_omega1(self, start, stop=0, steps=1):
        """
        Set the delta_omega1 param range

        Start required, default is no sampling (steps=1)
        """
        self.set_range_param("delta_omega1", start, stop, steps)

    def set_delta_t(self, start, stop=0, steps=1):
        """
        Set the delta_t param range

        Start required, default is no sampling (steps=1)
        """
        self.set_range_param("delta_t", start, stop, steps)

    def set_z0_offset(self, start, stop=0, steps=1):
        """
        Set the z0_offset param range

        Start required, default is no sampling (steps=1)
        """
        self.set_range_param("z0_offset", start, stop, steps)

    def set_z1_offset(self, start, stop=0, steps=1):
        """
        Set the z1_offset param range

        Start required, default is no sampling (steps=1)
        """
        self.set_range_param("z1_offset", start, stop, steps)

    def set_epsilon(self, start, stop=0, steps=1):
        """
        Set the epsilon param range

        Start required, default is no sampling (steps=1)
        """
        self.set_range_param("epsilon", start, stop, steps)

    def set_invert(self, base_value, allow_other=False):
        """
        Set the invert param range

        Start required, default is no sampling (steps=1)
        """
        self.set_binary_param("invert", base_value, allow_other=allow_other)

    def get_param_range(self, param_name):
        return self._params[param_name].get_range()

    def get_param_names(self):
        return self._params.keys()

    def get_params_dicts(self):
        """
        Converts the params to a list of dicts with name, type, values etc
        """
        return [param.to_dict() for param in self._params]

    def get_grid_dicts(self):
        return [
            dict(zip(self._params.keys(), prod))
            for prod in product(
                *[param.get_range() for param in self._params.values()]
            )
        ]


class PyBundleGridSampler(object):
    """
    Object to generate helical parameters for wide sampling on the bundle grid

    this object exists to patch that bundle grid sampler uses JD2 and
    therefore cannot be used in pyrosetta trivially.

    potential gotcha: helix numbers are 0-indexed, but the output for the
    grids is one-indexed (because Rosetta MakeBundle helices are one-indexed)
    """

    def __init__(self, num_helices, helix_length):
        self.num_helices = num_helices
        self.helix_length = helix_length
        self._helices = {
            num: HelixGridParam(num, self.helix_length)
            for num in range(1, self.num_helices + 1)
        }

    def set_helix_params_from_dicts(self, num, *params):
        for param in params:
            if param["type"] == "range":
                self._helices[num].set_range_param(
                    **{key: param[key] for key in param if key != "type"}
                )
            else:
                self._helices[num].set_binary_param(
                    **{key: param[key] for key in param if key != "type"}
                )

    def get_helix_grid(self, num):
        return self._helices[num].get_grid_dicts()

    def all_helix_grids(self):
        return [
            {**{"helix_length":self.helix_length,"num_helices": self.num_helices}, **dict(prod)}
            for prod in product(
                *[
                    product([num], helix.get_grid_dicts())
                    for num, helix in self._helices.items()
                ]
            )
        ]


def configure_helix(helix, length, **params):
    helix.set_helix_length(length)
    helix_bundle_params = {
        param_name: helix.calculator_op().parameter(
            helix.calculator_op().parameter_enum_from_name(param_name)
        )
        for param_name in params
    }
    for param in params:
        helix_bundle_params[param].set_value(params[param])


def helix_bundle_maker_wrapper(length, *helix_params, degrees=True):
    """
    Wrapper to make bundle generation a bit
    """
    bundle_maker = MakeBundle()
    bundle_maker.set_use_degrees(degrees)
    for i, param_set in enumerate(helix_params, 1):
        bundle_maker.add_helix()
        configure_helix(bundle_maker.helix(i),length,**param_set)
    return bundle_maker