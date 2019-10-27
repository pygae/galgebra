# -*- coding: utf-8 -*-

from sympy import Symbol

from ga import Ga
from mv import J, Jinv


def create_multivector(GA, name):
    blades = [1] + GA.blades_lst
    mv = GA.mv(0, 'scalar')
    for blade_index, blade in enumerate(blades):
        mv += Symbol('{name}[{i}]'.format(name=name, i=blade_index)) * blade
    return mv


_CLASS_TEMPLATE = '''# -*- coding: utf-8 -*-


class FlatMv(object):
    def __init__(self, coefs):
        assert len(coefs) == {class_blade_count}
        self.coefs = coefs

    def __getitem__(self, index):
        return self.coefs[index]
'''

_BINARY_OPERATOR_TEMPLATE = '''
    def __{op_name}__(self, other):
        x = self.coefs
        y = other.coefs
        return FlatMv([
            {op_list}
        ])
'''

_UNARY_METHOD_TEMPLATE = '''
    def {method_name}(self):
        x = self.coefs
        return FlatMv([
            {method_list}
        ])
'''

_BINARY_METHOD_TEMPLATE = '''
    def {method_name}(self, other):
        x = self.coefs
        y = other.coefs
        return FlatMv([
            {method_list}
        ])
'''


def format_class(class_blade_count):
    return _CLASS_TEMPLATE.format(class_blade_count=class_blade_count)


def format_op_name(name):
    return name


def format_op_list(mv):
    return ',\n            '.join(str(blade_coef) for blade_coef in mv.blade_coefs())


def format_binary_operator(name, mv):
    return _BINARY_OPERATOR_TEMPLATE.format(op_name=format_op_name(name), op_list=format_op_list(mv))


def format_method_name(name):
    return name


def format_method_list(mv):
    return ',\n            '.join(str(blade_coef) for blade_coef in mv.blade_coefs())


def format_unary_method(name, mv):
    return _UNARY_METHOD_TEMPLATE.format(method_name=format_method_name(name), method_list=format_method_list(mv))


def format_binary_method(name, mv):
    return _BINARY_METHOD_TEMPLATE.format(method_name=format_method_name(name), method_list=format_method_list(mv))


def format_geometric_algebra(GA):

    X = create_multivector(GA, 'x')
    Y = create_multivector(GA, 'y')

    flat_geometric_algebra = format_class(len(GA.blades_lst0))
    flat_geometric_algebra += format_binary_operator('add', X + Y)
    flat_geometric_algebra += format_binary_operator('sub', X - Y)
    flat_geometric_algebra += format_binary_operator('mul', X * Y)
    flat_geometric_algebra += format_binary_operator('and', Jinv(J(X) ^ J(Y)))
    flat_geometric_algebra += format_binary_operator('xor', X ^ Y)
    flat_geometric_algebra += format_binary_operator('lshift', X << Y)
    flat_geometric_algebra += format_binary_operator('rshift', X >> Y)
    flat_geometric_algebra += format_binary_method('meet', Jinv(J(X) ^ J(Y)))
    flat_geometric_algebra += format_binary_method('join', X ^ Y)
    flat_geometric_algebra += format_unary_method('rev', X.rev())

    return flat_geometric_algebra


def flatten(flat_ga_module, mv):
    return flat_ga_module.FlatMv([blade_coef for blade_coef in mv.blade_coefs()])


def expand(GA, flat_mv):
    assert len(flat_mv.coefs) == len(GA.blades_lst0)
    mv = GA.mv(0, 'scalar')
    for blade_coef, blade in zip(flat_mv.coefs, GA.blades_lst0):
        mv += blade_coef * blade
    return mv


if __name__ == "__main__":
    GA = Ga('e*1|2|3', g=[1, 1, 1])
    print format_geometric_algebra(GA)
