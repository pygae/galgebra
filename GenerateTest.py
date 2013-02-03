#!/usr/bin/python
#GenerateTest.py

import os,sys,GAdir
import re as regrep

if GAdir.GA == 'GA':
    from GAdebug import oprint
else:
    from sympy.galgebra.GAdebug import oprint

def get_function_strings(pstr):
    iold = 0
    funct_dict = {}
    post = ''
    pre  = ''
    while True:
        i = pstr.find('\ndef ',iold)
        if i == -1:
            post = pstr[iold-4:]
            break
        if iold > 0:
            icolon = pstr.find(':',iold-4)
            funct_dict[pstr[iold:icolon]] = pstr[iold-4:i-1]
        else:
            pre = pstr[:i-1]
        iold = i+5
    include_fct = []
    exclude_fct = []
    for key in funct_dict:
        if key in post:
            ikey = post.find(key)
            if post[ikey-1] != '#':
                include_fct.append(key)
            else:
                exclude_fct.append(key)
    return(funct_dict,pre,include_fct,exclude_fct)

def ExecuteProgramString(pstr):
    tmp_file = open('tmp.py','w')
    tmp_file.write(pstr)
    tmp_file.close()
    os.system('python tmp.py > tmp.out')
    tmp_file = open('tmp.out','r')
    out_lst = tmp_file.readlines()
    tmp_file.close()
    os.system('rm tmp.py tmp.out')
    return(out_lst)

def get_print_statments(fstr):
    plst = []
    istart = fstr.find('\n',2)+1
    while True:
        iprint = fstr.find('print ',istart)
        iFmt   = fstr.find('.Fmt(',istart)
        if iprint < 0 and iFmt < 0:
            break
        if (iprint < 0 and iFmt >= 0) or iFmt < iprint:
            ibeg = fstr.rfind('\n',iFmt)
            iend = fstr.find('\n',ibeg+1)
            plst.append(fstr[ibeg+1:iend])
        if (iprint >= 0 and iFmt < 0) or iprint < iFmt:
            iend = fstr.find('\n',iprint)
            plst.append(fstr[iprint:iend])
        istart = iend+1
    return(plst)


program_name = 'manifold_test.py'

program_file = open(program_name,'r')
program_str = program_file.read()
program_file.close()

(fdict,pre,include_fct,exclude_fct) = get_function_strings(program_str)

revised_program = pre+'\n\n'

for key in fdict:
    if key not in include_fct and key not in exclude_fct:
        revised_program += fdict[key]+'\n\n'
for key in fdict:
    if key in include_fct:
        revised_program += fdict[key]+'\n\n'
revised_program += '\n'
for f in include_fct:
    revised_program += '\n'+f

revised_program = revised_program.replace('oprint','#oprint')

print revised_program

out_lst = ExecuteProgramString(revised_program)

oprint('Revised Output',out_lst)

test_program = pre+'\n\n'

iresult = 0

for key in fdict:
    if key not in include_fct and key not in exclude_fct:
        test_program += fdict[key]+'\n\n'
for key in fdict:
    if key in include_fct:
        tmp = fdict[key].replace('def '+key,'def test_'+key)
        tmp = tmp.replace('return','')
        plst = get_print_statments(tmp)
        for p in plst:
            print p

#test_program += '\n'

#test_program = test_program.replace('oprint','#oprint')

#print test_program

