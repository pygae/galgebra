from sympy import *
from sympy.galgebra.GA import *
from sympy import cos, sin
from numpy import linalg
import sys,os
import re as regrep
import subprocess
#op_cntrct = regrep.compile(r'(([A-Za-z0-9\_\#]+)(\||<<|>>)([A-Za-z0-9\_\#]+))')
op_cntrct = regrep.compile(r'(([A-Za-z0-9\_\#]+)(\||<|>)([A-Za-z0-9\_\#]+))')
op_wedge  = regrep.compile(r'(([A-Za-z0-9\_\#]+)[\^]{1}([A-Za-z0-9\_\#]+)([\^]{1}([A-Za-z0-9\_\#]+))*)')

ops = r'[\^\|\<\>]+'

#ops_search = regrep.compile(r'(\^|\||<<|>>)+')
ops_search = regrep.compile(r'(\^|\||<|>)+')

parse_paren_calls = 0


def print_lst(lst):  # This is not used
    ilst = 0
    for x in lst:
        print str(ilst)+': '+str(x)
        ilst += 1
    return

def contains_interval(interval1,interval2): #interval1 inside interval2
    if interval1[0] > interval2[0] and interval1[1] < interval2[1]:
        return(True)
    else:
        return(False)

def parse_paren(line):
    global parse_paren_calls
    parse_paren_calls += 1

    if ('(' not in line) or (')' not in line):
        return([[[line]]])
    level = 0
    max_level = 0
    ich = 0
    paren_lst = []
    for ch in line:
        if ch == '(':
            level += 1
            paren_lst.append([level,ich])
        if ch == ')':
            if level < 1:
                sys.stderr.write('Mismathed Parenthesis in: '+line+'\n')
                sys.exit(1)
            paren_lst.reverse()
            iparen = 0
            for elem in paren_lst:
                if elem[0] == level:
                    paren_lst[iparen].append(ich)
                    break
                iparen += 1
            paren_lst.reverse()
            level -= 1
        max_level = max(max_level,level)
        ich += 1
    if level != 0:
        sys.stderr.write('Mismathed Parenthesis in: '+line+'\n')
        sys.exit(1)
    if max_level > 0:
        level_lst = []
        for x in range(max_level+1):
            level_lst.append([])
        for group in paren_lst:
            level_lst[group[0]].append(group[1:])
        ilevel = max_level
        while ilevel > 1:
            level = level_lst[ilevel]
            level_down = level_lst[ilevel-1]
            igroup = 0
            for group in level:
                igroup_down = 0
                for group_down in level_down:
                    if contains_interval(group,group_down):
                        level_lst[ilevel][igroup].append(igroup_down)
                    igroup_down += 1
                igroup += 1
            ilevel -= 1
        ilevel = 1
        for level in level_lst[1:]:
            igroup = 0
            for group in level:
                token = '#'+str(parse_paren_calls)+'_'+str(ilevel)+'_'+str(igroup)+'#'
                level_lst[ilevel][igroup].append(line[group[0]:group[1]+1])
                level_lst[ilevel][igroup].append(token)
                igroup += 1
            ilevel += 1
        ilevel = 1
        for level in level_lst[1:]:
            igroup = 0
            for group in level:
                group.append(group[-2])
                level_lst[ilevel][igroup] = group
                igroup += 1
            ilevel += 1
        ilevel = max_level
        while ilevel > 1:
            igroup = 0
            for group in level_lst[ilevel]:
                group_down = level_lst[ilevel-1][group[2]]
                replace_text = group_down[-1].replace(group[-3],group[-2])
                level_lst[ilevel-1][group[2]][-1] = replace_text
                igroup += 1
            ilevel -= 1
        for group in level_lst[1]:
            line = line.replace(group[2],group[3])
        ilevel = 1
        level_lst[0] = [[line]]
    return(level_lst)

def unparse_paren(level_lst):
    line = level_lst[0][0][0]
    for level in level_lst[1:]:
        for group in level:
            new_string = group[-1]
            if new_string[:2] == '((' and new_string[-2:] == '))':
                new_string = new_string[1:-1]
            line = line.replace(group[-2],new_string)
    return(line)

def sub_paren(s):
    string = s.group(0)
    return('(%s)' % string)

def add_paren(line,re_exprs):
    paren_flg = False
    if (line[0] == '(') and (line[-1] == ')'):
        paren_flg = True
        line = line[1:-1]
    if ('(' in line) or (')' in line):
        line_levels = parse_paren(line)
        ilevel = 0
        for level in line_levels:
            igroup = 0
            for group in level:
                group[-1] = regrep.sub(re_exprs,sub_paren,group[-1])
                line_levels[ilevel][igroup] = group
                igroup += 1
            ilevel += 1
        line = unparse_paren(line_levels)
    else:
        line = regrep.sub(re_exprs,sub_paren,line)
    if paren_flg:
        line = '('+line+')'
    return(line)

def parse_line(line):
    level_lst = parse_paren(line)
    ilevel = 0
    for level in level_lst:
        igroup = 0
        for group in level:
            string = group[-1]
            string = add_paren(string,op_cntrct)
            string = add_paren(string,op_wedge)
            level_lst[ilevel][igroup][-1] = string
            igroup += 1
        ilevel += 1
    line = unparse_paren(level_lst)
    return(line)

"""
def int_replace(s):
    s1 = s.group(1)
    s2 = s.group(2)
    s3 = s.group(3)
    s4 = s.group(4)
    if len(s1) > 0:
        return('%s' % s.group(0) )
    else:
        return('S(%s)%sS(%s)' % (s2,s3,s4))

int_search = regrep.compile(r'([A-Za-z\_]*)([0-9]+)(\||<<|>>|<|>|\^)([0-9]+)')
"""

def int_replace(s):
    s0 = s.group(0)
    s1 = s.group(1)
    #s2 = s.group(2)
    #sys.stderr.write('('+s0+','+s1+','+s2+')\n')
    if len(s1) == 0:
        return('S(%s)' % s0)
    else:
        return('%s' % s0)

int_search = regrep.compile(r'([A-Za-z\_]*)([0-9]+)')

def parse_integers(line):
    #sys.stderr.write('parse_integers (input): '+line+'\n')
    newline = regrep.sub(int_search,int_replace,line)
    #sys.stderr.write('parse_integers (output): '+newline+'\n')
    return(newline)

def parse_all(line):
    if "'''" in line:
        return("")
    line_split = line.split('@')
    new_line = line_split[0]
    for i in range(1,len(line_split)):         # Process each segment starting after an '@'.
        level  = 0
        for j in range(len(line_split[i])):    # Scan segment
            ch = line_split[i][j]
            if ch == ',' or (ch == ')' and level == 0):
                break
            if ch == '(':
                level += 1
            if ch == ')':
                level -= 1
        new_line += parse_integers(parse_line(line_split[i][:j])) + line_split[i][j:]
    return(new_line)

############################ Routines for LAGA ###########################

def preprocess(input_file_name):
    ops = r'[\^\|\<\>]+'
    #ops_search = regrep.compile(r'(\^|\||<<|>>)+')
    ops_search = regrep.compile(r'(\^|\||<|>)+')

    output_file_name = input_file_name + '.py'
    output_file = open(output_file_name,'w')
    path = os.path.abspath(output_file_name)

    output_file.write('from sympy import *\n')
    output_file.write('from sympy.galgebra.GA import *\n')
    output_file.write('from sympy.galgebra.latex_ex import *\n')
    output_file.write('from sympy.galgebra.laga import *\n')
    output_file.write('set_main(sys.modules[__name__])\n')

    found_end = False
    for line in open(input_file_name,'r'):
        if found_end:
            line = parse_all(line)
            output_file.write(line)
        else:
            found_end = "'''" in line
    output_file.close()

    theproc = subprocess.Popen([sys.executable, output_file_name])
    theproc.communicate()
    #os.system('python ' + output_file_name)  #AM  This is deprecated in 2.6
    #os.system('idle '+output_file_name)

def rank(M):                 # Return rank of matrix M.
    return len(M.rref()[1])

def printeigen(M):    # Print eigenvalues, multiplicities, eigenvectors of M.
    evects = M.eigenvects()
    for i in [0,len(evects)-1]:
        print 'Eigenvalue =', evects[i][0], '  Multiplicity =', evects[i][1],
        print 'Eigenvector =',
        result = '['
        for j in [0,len(evects)-1]:
            result += str(trigsimp(evects[i][2][0][j]).evalf(3))
            if j != len(evects)-1:
                result += ', '
        result += ']'
        print result
    print

def printGS(M):       # Print Gram-Schmidt output.
    result = '['
    for i in [0,len(M)-1]:
        result += '['
        for j in [0,len(M)-1]:
            result += str(trigsimp(M[i][j]).evalf(3))
            if j != len(M)-1:
                result += ', '
        result += ']'
        if j != len(M)-1:
            result += ' '
    result += ']'
    print result

def printrref(matrix, vars = "xyzuvwrs"):   # Print rref of matrix with variables.
    rrefmatrix = matrix.rref()[0]
    rows, cols = rrefmatrix.shape
    if len(vars) < cols-1:
        print 'Not enough variables.'
        return
    for i in range(rows):
        result = ''
        for j in range(cols-1):
            result += str(rrefmatrix[i,j]) + vars[j]
            if j != cols-2:
                result += ' + '
        result += ' = ' + str(rrefmatrix[i,cols-1])
        print result

def correlation(u,v,dec=3):    # Compute the correlation coefficient of vectors u and v.
    rows, cols = u.shape
    uave = 0
    vave = 0
    for i in range(rows):
        uave += u[i]
        vave += v[i]
    uave = uave/rows
    vave = vave/rows
    for i in range(rows):
        u[i] -= uave
        v[i] -= vave
    return u.dot(v)/(u.norm()*v.norm()).evalf(dec)

def norm(B):          # Compute |B|.
    return(sqrt(abs((B*B.rev()).scalar())))

def inv(B):           # Invert blade or versor B ONLY.
    Bnorm = B*B.rev()
    if Bnorm.is_scalar():
        invB = B.rev()/Bnorm.scalar()
        return(invB)
    else:
        sys.stderr.write('Cannot calculate inverse of '+str(B)+' since \n'+\
                         'B*Brev() = '+str(Bnorm)+' is not a scalar.\n')

def proj(B,A):        # Project blade A on blade B.
    AdotB = A<B
    invB = inv(B)
    result = (A<B)*inv(B)
    result.trigsimp()
    return( result )

def rot(itheta,A):    # Rotate blade A by angle itheta.
    theta = norm(itheta)
    i = itheta/theta
    result = (cos(theta/2) - i*sin(theta/2)) * A * (cos(theta/2) + i*sin(theta/2))
    result.trigsimp()
    return( result )

def refl(B,A):
    j = B.is_pure()
    k = A.is_pure()
    if j > -1 and k > -1:       # Reflect blade A in blade B.
        result = (-1)**(j*(k+1)) * B * A * inv(B)
        result.trigsimp()
        return( result )
    else:
        sys.stderr.write('Can only reflect blades')

def dual(M):
        return(M*MV.Iinv)

def cross(M1, M2):
    return(dual(M1^M2))

def ScalarFunction(TheFunction):
        return(MV() + TheFunction)


