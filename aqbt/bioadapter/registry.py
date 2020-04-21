import inspect
from typing import Any
from typing import Callable
from typing import List
from typing import Type
from typing import Union

import networkx as nx
from tqdm import tqdm


class ConversionError(Exception):
    """Generic conversion error."""


class RegistryError(Exception):
    """Generic error raised during registration."""


class BioAdapter:
    """Conversion type registry."""

    def __init__(self):
        """Initialize the registry."""
        self.g = nx.DiGraph()

    def _add_function(
        self, frm: str, to: str, func: Callable, weight: int, override: bool
    ):
        """Add a function to the register.

        :param frm: 'from' key
        :param to: 'to' key
        :param func: function
        :return: None
        """
        if not override and frm in self.g and to in self.g[frm]:
            raise RegistryError("Conversion {}->{} already exists".format(frm, to))
        self.g.add_edge(frm, to, function=func, weight=weight)

    def _get_path(self, frm: str, to: str) -> List[Callable]:
        """Find the shortest path from the 'from' to the 'to' node, returning
        the list of functions in that path.

        :param frm: 'from' node
        :param to: 'to' node
        :return: list of functions
        """
        path = self.path(frm, to)
        funcs = []
        for n1, n2 in nx.utils.pairwise(path):
            funcs.append(self.g[n1][n2]["function"])
        return funcs, path

    def path(self, frm, to):
        if frm not in self.g:
            raise ConversionError(
                "Registry has no key '{}'. Keys are {}".format(frm, list(self.g.nodes))
            )
        if to not in self.g:
            raise ConversionError(
                "Registry has no key '{}'. Keys are {}".format(to, list(self.g.nodes))
            )
        if to in self.g[frm]:
            return [frm, to]
        else:
            try:
                return nx.shortest_path(self.g, frm, to, weight="weight")
            except nx.NetworkXNoPath:
                raise ConversionError(
                    "There is no conversion path from '{}' to '{}'".format(frm, to)
                )

    def convert(
        self, data, frm=None, to=None, pbar=None, no_raise=False, **arguments
    ) -> Any:
        """Convert data from one type to another. If 'frm' is not provided, the
        type will be infered from the type of the 'data' provided.

        If functions in the shortest path algorith require arguments,
        these arguments should be provided in in the 'arguments' kwargs.

        Note that conflicts can arise in the arugment names of the shortest path.
        It is good practice to keep unique names for argument when registering
        functions.

        For example, to convert from a Bio.SeqRecord to a Benchling API DNASequence
        object, we must include the benchling_session argument for one of the conversion
        functions. The convert function will only use the argument when necessary.

        .. code-block:: python

            registry.convert(seqrecord, to="DNASequence", benchling_session=session)

        :param data: data to convert
        :param frm: the optional type as a string
        :param to: the require type to convert to
        :param pbar: if provided or True, list dispatches will use a tqdm progress bar.
        :param no_raise: if True, conversion function will return None instead
                            of raising a ConversionError.
        :param arguments: optional arguments conversion functions may require.
        :return: the converted data
        :raises: TypeError when functions are missing arguments.
        :raises: ConversionError if there is no path between the conversion types
        """
        if arguments is None:
            arguments = {}

        if to is None:
            raise ConversionError("Must provide 'to' string")
        if isinstance(to, type):
            to = to.__name__
        if isinstance(frm, type):
            frm = frm.__name__
        if frm is None:
            frm = data.__class__.__name__
            if frm == "list":
                frm = "list({})".format(data[0].__class__.__name__)
                to = "list({})".format(to)

        converted = data
        if frm == to:
            return converted
        func_paths, path = self._get_path(frm, to)
        if not func_paths:
            raise ConversionError("No path from '{}' to '{}'".format(frm, to))
        for i, f in enumerate(func_paths):
            if hasattr(f, "__source__"):
                argspec = inspect.getfullargspec(f.__source__)
            else:
                argspec = inspect.getfullargspec(f)
            args = argspec.args[1:]
            if argspec.defaults:
                required_args = argspec.args[1 : -len(argspec.defaults)]
            else:
                required_args = argspec.args[1:]
            missing = [r for r in required_args if r not in arguments]
            if missing:
                raise TypeError(
                    "{func}() is missing the following arguments: {args}."
                    " Did you mean `{this_func}(data, {key}=<value>)`?".format(
                        func=f.__name__,
                        this_func=self.convert.__name__,
                        args=", ".join(missing),
                        key=missing[0],
                    )
                )
            pulled_kwargs = {k: arguments[k] for k in args if k in arguments}
            if hasattr(f, "__source__"):
                pulled_kwargs["pbar"] = pbar
                pulled_kwargs["no_raise"] = no_raise
            try:
                converted = f(converted, **pulled_kwargs)
            except Exception as e:
                raise e
                if no_raise:
                    return None
                fname = f.__name__
                if hasattr(f, "__source__"):
                    fname = "[{}(x) for x in data]".format(f.__source__.__name__)
                raise ConversionError(
                    "An error occured during conversion:"
                    "\n\tin method: {func}()"
                    "\n\tduring step {step}, {path}"
                    "\n\t{errcls}: {err}".format(
                        step=i,
                        func=fname,
                        path=" > ".join(path),
                        errcls=e.__class__.__name__,
                        err=str(e),
                    )
                ) from e
        return converted

    def has_path(self, frm: str, to: str):
        try:
            path = nx.shortest_path(self.g, frm, to, weight="weight")
            if path:
                return True
        except nx.NetworkXNoPath:
            pass
        return False

    def register(
        self,
        input_name: Union[str, Type],
        output_name: Union[str, Type],
        many: bool = False,
        weight: int = 1,
        override: bool = False,
    ):
        """Decorator to register a function as a conversion function.

        .. code-block::

            from bioadapter import registry

            @registry("type1", "type2")
            def foo(data):
                # conversion code
                pass

            # many=True will also autoregister 'list(type1)' and 'list(type2)'
            # and autodispatch the registered function as a list comprehension
            @registry("type1", "type2", many=True)
            def bar(data):
                # conversion code
                pass


            # you can manually register a 'many' conversion function as well
            @registery('list(type1)', 'list(type2')
            def baz(data):
                # conversion code
                pass


        :param input_name: the name of the input data type
        :param output_name: the name of the ouptut data type
        :param many: if True, will also auto-register a conversion of `list` types.
                        The input and output names for this type will be called
                        `list(<input_name>)` and `list(<output_name>)` respectively.
        :return: decorated function
        """
        if not isinstance(input_name, str):
            input_name = input_name.__name__
        if not isinstance(output_name, str):
            output_name = output_name.__name__

        def registered_function(f):
            self._add_function(
                input_name, output_name, f, weight=weight, override=override
            )
            if many:

                def dispatch_list(*args, **kwargs):
                    pbar = kwargs["pbar"]
                    no_raise = kwargs["no_raise"]
                    del kwargs["pbar"]
                    del kwargs["no_raise"]
                    if pbar is True:
                        iterator = tqdm(args[0])
                    elif pbar:
                        iterator = pbar(args[0])
                    else:
                        iterator = args[0]
                    results = []
                    for a in iterator:
                        try:
                            results.append(f(a, *args[1:], **kwargs))
                        except Exception as e:
                            if no_raise:
                                results.append(None)
                            else:
                                raise e
                    return results

                dispatch_list.__source__ = f
                self._add_function(
                    "list({})".format(input_name),
                    "list({})".format(output_name),
                    dispatch_list,
                    weight=weight,
                    override=override,
                )
            return f

        return registered_function

    def __call__(self, input_name, output_name, many=False, weight=1, override=False):
        return self.register(
            input_name, output_name, many=many, weight=weight, override=override
        )


bioadapter = BioAdapter()

convert = bioadapter.convert
