from functools import reduce
from sympy import simplify, sqrt, Rational, Symbol, symbols
from galgebra.ga import Ga
from galgebra.mv import Mv
from galgebra.printer import Format, xtex
#Format()
(GA, e0, e1, e2, e3) = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

(x1,x2,x3) = symbols('x1 x2 x3',real=True)
x = x1*e1+x2*e2+x3*e3
print('x =',x)
a = e1+2*e2+3*e3
print('a =', a)
y = x+a
print('y =',y)
z = y.subs([x1,x2,x3],[1,2,3])
print('z =',z)
w = GA.mv('w','vector')

print('w =',w)
w = w.subs([w__0,w__1,w__2,w__3],[1,2,3,4])
w = GA.mv('w','vector')



#xtex()



"""
t = Symbol('t')
x_1 = Symbol('x_1')
        x_2 = Symbol('x_2')
        x_3 = Symbol('x_3')
        x = e_0 + x_1 * e_1 + x_2 * e_2 + x_3 * e_3

        # x(t)
        x_t = (e_0 + e_1 - 3 * e_2) + t * (e_1 + 2 * e_2 - e_3)

        d = e_1 + 2 * e_2 - e_3
        d_inv = d.inv()

        # t(x)
        #t_x = (x - (e_0 + e_1 - 3 * e_2)) / (e_1 + 2 * e_2 - e_3)
        t_x = (x - (e_0 + e_1 - 3 * e_2)) * d_inv

        for i in range(11):
            t_value = Rational(i, 10)
            x_value = x_t.subs({t: t_value})
            print('t_value =', t_value)
            print('x_value =', x_value)
            print('x_value ^ L =',x_value ^ L)
            x_c = x_value.blade_coefs()
            #_x_new = ((x_c[1]+1)*e_1+(x_c[2]-3)*e_2+x_c[3]*e_3)*d_inv
            #print(t_x_new)
            print(t_x)
            print(x_value.blade_coefs([e_1, e_2, e_3]))
            print(t_x.subs(zip([x_1, x_2, x_3], x_value.blade_coefs([e_1, e_2, e_3]))))
            self.assertEquals(x_value ^ L, 0)
            self.assertEquals(t_x.subs(zip([x_1, x_2, x_3], x_value.blade_coefs([e_1, e_2, e_3]))), t_value)
"""
