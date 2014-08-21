#!/usr/bin/python
#fig2eps.py

import os,sys

preamble = '\\documentclass[12pt,letter,fleqn]{report}\n'+\
           '\\pagestyle{empty}\n'+\
           '\\usepackage[latin1]{inputenc}\n'+\
           '\\usepackage[pdftex]{geometry,graphics,graphicx,color}\n'+\
           '\\usepackage{amsmath}\n'+\
           '\\usepackage{bm}\n'+\
           '\\usepackage{color}\n'+\
           '\\usepackage{amsfonts}\n'+\
           '\\usepackage{amssymb}\n'+\
           '\\include{bookmacros}\n'+\
           '\\begin{document}\n'
postscript = '\\end{document}\n'

def make_PDF(name,delx):
    print name,delx
    pstex  = name+'.pdf'
    pstex_t = name+'.pdf_t'
    tmp_name = name+'_tmp'
    tmp_tex_name = tmp_name+'.tex'

    body = preamble+'\\hspace{'+delx+'in}\\input{'+pstex_t+'}\n'+postscript

    texfile = open(tmp_tex_name,'w')
    texfile.write(body)
    texfile.close()
    os.system('pdflatex '+tmp_tex_name)
    os.system('evince '+tmp_name+'.pdf')
    os.system('rm '+tmp_name+'.tex '+tmp_name+'.log '+tmp_name+'.aux')
    return(tmp_name+'.pdf')

if __name__ == '__main__':

    name = sys.argv[1]
    delx = '0'

    while delx !='':
        print delx
        if delx in ('q','Q'):
            sys.exit(0)
        print delx
        tmp_PDF = make_PDF(name,delx)
        delx = raw_input('Action ([PDF OK],delx value,q/Q):')

    os.system('pdfcrop.pl '+tmp_PDF+' '+name+'_crop.pdf')
    os.system('rm '+tmp_PDF)
    sys.exit(0)

