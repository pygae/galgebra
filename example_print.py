import sys

prog_str = ''

def get_program():
    global prog_str
    prog_file = open(sys.argv[0],'r')
    prog_str = prog_file.read()
    prog_file.close()
    return

def print_function():
    global prog_str
    fct_name = str(sys._getframe(1).f_code.co_name)
    ifct = prog_str.find('def '+fct_name)
    iend = prog_str.find('def ',ifct+4)
    tmp_str = prog_str[ifct:iend-1]
    tmp_str = tmp_str.replace('print_function()\n','')
    tmp_str = tmp_str.replace('def '+fct_name+'():\n','')
    fct_name = fct_name.replace('_',' ')
    print '\n'+80*'*'
    print '\nCode for',fct_name+':'
    print tmp_str
    print 'Code output:\n'
    return

