# -*- coding: utf-8 -*-
"""
IPython (0.11) extension for physical quantity input.
See README.txt for usage examples.

Author: Georg Brandl <georg@python.org>.
This file has been placed in the public domain.
"""

import re
import sys
from math import pi

def parse_operators(line):
    return(line)

class QTransformer(object):
    # XXX: inheriting from PrefilterTransformer as documented gives TypeErrors,
    # but apparently is not needed after all
    priority = 99
    enabled = True
    def transform(self, line, continue_prompt):
        if line[0] == '>':
            line = line[1:]
            line = parse_operators(line)
        return line

q_transformer = QTransformer()

def load_ipython_extension(ip):
    ip.prefilter_manager.register_transformer(q_transformer)

    # active true float division
    exec ip.compile('from __future__ import division', '<input>', 'single') \
        in ip.user_ns

    print '***sympy geometric algebra extensions activated.'

def unload_ipython_extension(ip):
    ip.prefilter_manager.unregister_transformer(q_transformer)
