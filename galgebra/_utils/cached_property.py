class cached_property:
    """ A property that is only computed once """
    def __init__(self, getter):
        self.fget = getter
        self.__name__ = getter.__name__
        self.__doc__ = getter.__doc__

    def __set_name__(self, owner, name):
        self.__name__ = name
        # copy across the annotation for sphinx
        try:
            prop_ann = self.fget.__annotations__['return']
        except KeyError:
            return

        try:
            owner_anns = owner.__annotations__
        except AttributeError:
            owner_anns = owner.__annotations__ = {}
        owner_anns[name] = prop_ann

    def __get__(self, obj, cls):
        if obj is None:
            return self
        val = self.fget(obj)
        # this entry hides the _cached_property
        setattr(obj, self.__name__, val)
        return val
