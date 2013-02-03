#!/usr/bin/env python
#
# pretty - A miniature library that provides a Python print and stdout
# wrapper that makes colored terminal text easier to use (eg. without
# having to mess around with ANSI escape sequences). This code is public
# domain - there is no license except that you must leave this header.
#
# Copyright (C) 2008 Brian Nez <thedude at bri1 dot com>
#

import sys

codeCodes = {
        'black':        '0;30',         'bright gray':  '0;37',
        'blue':         '0;34',         'white':                '1;37',
        'green':        '0;32',         'bright blue':  '1;34',
        'cyan':         '0;36',         'bright green': '1;32',
        'red':          '0;31',         'bright cyan':  '1;36',
        'purple':       '0;35',         'bright red':   '1;31',
        'yellow':       '0;33',         'bright purple':'1;35',
        'dark gray':'1;30',             'bright yellow':'1;33',
        'normal':       '0'
}

def printc(text, color):
        """Print in color."""
        print "\033["+codeCodes[color]+"m"+text+"\033[0m"

def writec(text, color):
        """Write to stdout in color."""
        sys.stdout.write("\033["+codeCodes[color]+"m"+text+"\033[0m")

def switchColor(color):
        """Switch console color."""
        sys.stdout.write("\033["+codeCodes[color]+"m")

if __name__ == '__main__':
        print "Welcome to the test routine!"
        print "I will now try to print a line of text in each color."
        for color in codeCodes.keys():
                writec("Hello, world!", color)
                print "\t", color
