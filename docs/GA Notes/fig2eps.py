#!/usr/bin/python
#fig2eps.py

import os,sys

preamble = '\\documentclass[12pt,letter,fleqn]{report}\n'+\
           '\\pagestyle{empty}\n'+\
           '\\usepackage[latin1]{inputenc}\n'+\
           '\\usepackage[pdflatex]{geometry,graphics,graphicx,color}\n'+\
           '\\usepackage{amsmath}\n'+\
           '\\usepackage{bm}\n'+\
           '\\usepackage{color}\n'+\
           '\\usepackage{amsfonts}\n'+\
           '\\usepackage{amssymb}\n'+\
           '\\newcommand{\\bfrac}[2]{\\displaystyle\\frac{#1}{#2}}\n'+\
           '\\newcommand{\\lp}{\\left (}\n'+\
           '\\newcommand{\\rp}{\\right )}\n'+\
           '\\newcommand{\\half}{\\frac{1}{2}}\n'+\
           '\\newcommand{\\llt}{\\left <}\n'+\
           '\\newcommand{\\rgt}{\\right >}\n'+\
           '\\newcommand{\\abs}[1]{\\left |{#1}\\right | }\n'+\
           '\\newcommand{\\pdiff}[2]{\\bfrac{\\partial {#1}}{\\partial {#2}}}\n'+\
           '\\newcommand{\\lbrc}{\\left \\{}\n'+\
           '\\newcommand{\\rbrc}{\\right \\}}\n'+\
           '\\newcommand{\\W}{\\wedge}\n'+\
           "\\newcommand{\\prm}[1]{{#1}'}\n"+\
           '\\newcommand{\\ddt}[1]{\\bfrac{d{#1}}{dt}}\n'+\
           '\\newcommand{\\R}{\\dagger}\n'+\
           '\\begin{document}\n'
postscript = '\\end{document}\n'

if __name__ == '__main__':

    name = sys.argv[1]
    #name = 'space1d'
    pstex  = name+'.pdf'
    pstex_t = name+'.pdf_t'
    tmp_name = name+'_tmp'
    tmp_tex_name = tmp_name+'.tex'

    body = preamble+'\\input{'+pstex_t+'}\n'+postscript

    texfile = open(tmp_tex_name,'w')
    texfile.write(body)
    texfile.close()
    os.system('pdflatex '+tmp_tex_name)
    os.system('acroread '+tmp_name+'.pdf')
    op = raw_input('Action ([File OK]):')
    if op = '\n':
        print 'File OK'
    else:
        print 'Do More'

    #os.system('rm '+tmp_name+'.*')
    sys.exit(0)

