from sympy import Symbol, symbols, sin, cos, Rational, expand, simplify, collect
from printer import Format, Eprint, Get_Program, Print_Function,xpdf
from ga import Ga, one, zero
from mv import Com, Nga

def F(x):
    global n,nbar
    Fx = ((x*x)*n+2*x-nbar) / 2
    return(Fx)

def make_vector(a,n = 3, ga=None):
    if isinstance(a,str):
        v = zero
        for i in range(n):
            a_i = Symbol(a+str(i+1))
            v += a_i*ga.basis[i]
        v = ga.mv(v)
        return(F(v))
    else:
        return(F(a))


def making_a_circle():
    Print_Function()
    global n, nbar
    g='1 0 0 0 0, \
       0 1 0 0 0, \
       0 0 1 0 0, \
       0 0 0 0 2, \
       0 0 0 2 0'
    c3d = Ga('e_1 e_2 e_3 n \\bar{n}',g=g)
    (e1,e2,e3,n,nbar) = c3d.mv()
    e = n + nbar

    A=make_vector(e1/2,ga=c3d)
    B=make_vector(2*e1,ga=c3d)
    C=make_vector((4.0/5)*e1 + (3.0/5)*e2,ga=c3d)

    A=8*A
    B=2*B
    C=10*C

    print 'F(a) = ',A
    print 'F(b) = ',B
    print 'F(c) = ',C
    print '#Circle through a,b,c'

    print 'A^B = ',(A^B)/3
    print '#Circle triveector:'
    print '(A^B)^C = ',A^B^C
    print 'B^C = ', (B^C)/3
    print '#The same circle trivector as before, different computation order and scaled: '
    print 'A^(B^C) = ',(A^B^C)/18
    print '#Haven\'t figured out how to format the coefficients nicely'
    print '# Are A,B,C all on the line? Let\'s check one:'
    print 'A^B^C^D =',A^B^C^C

    return


def main():
    Get_Program()

    Eprint()

    #Format()
    making_a_circle()

    #xpdf()
    print 'done'
    return

if __name__ == "__main__":
    main()
