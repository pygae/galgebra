from __future__ import print_function
from sympy import symbols, sin, cos, sinh, cosh, trigsimp, S
from galgebra.printer import Format, xpdf, Print_Function
from galgebra.ga import Ga

def Product_of_Rotors():
    Print_Function()
    (na,nb,nm,alpha,th,th_a,th_b) = symbols('n_a n_b n_m alpha theta theta_a theta_b',\
                                    real = True)
    g = [[na, 0, alpha],[0, nm, 0],[alpha, 0, nb]] #metric tensor
    """
    Values of metric tensor components
    [na,nm,nb] = [+1/-1,+1/-1,+1/-1]  alpha = ea|eb
    """
    (g3d, ea, em, eb) = Ga.build('e_a e_m e_b', g=g)
    print('g =',g3d.g)
    print(r'%n_{a} = \bm{e}_{a}^{2}\;\;n_{b} = \bm{e}_{b}^{2}\;\;n_{m} = \bm{e}_{m}^{2}'+\
        r'\;\;\alpha = \bm{e}_{a}\cdot\bm{e}_{b}')
    (ca,cb,sa,sb) = symbols('c_a c_b s_a s_b',real=True)
    Ra = ca + sa*ea*em  # Rotor for ea^em plane
    Rb = cb + sb*em*eb  # Rotor for em^eb plane
    print(r'%\mbox{Rotor in }\bm{e}_{a}\bm{e}_{m}\mbox{ plane } R_{a} =',Ra)
    print(r'%\mbox{Rotor in }\bm{e}_{m}\bm{e}_{b}\mbox{ plane } R_{b} =',Rb)
    Rab = Ra*Rb  # Compound Rotor
    """
    Show that compound rotor is scalar plus bivector
    """
    print(r'%R_{a}R_{b} = S+\bm{B} =', Rab)
    Rab2 = Rab.get_grade(2)
    print(r'%\bm{B} =',Rab2)
    Rab2sq = Rab2*Rab2  # Square of compound rotor bivector part
    Ssq = (Rab.scalar())**2  # Square of compound rotor scalar part
    Bsq =  Rab2sq.scalar()
    print(r'%S^{2} =',Ssq)
    print(r'%\bm{B}^{2} =',Bsq)
    Dsq = (Ssq-Bsq).expand().simplify()
    print('%S^{2}-B^{2} =', Dsq)
    Dsq = Dsq.subs(nm**2,S(1))  # (e_m)**4 = 1
    print('%S^{2}-B^{2} =', Dsq)
    Cases = [S(-1),S(1)]  # -1/+1 squares for each basis vector
    print(r'#Consider all combinations of $\bm{e}_{a}^{2}$, $\bm{e}_{b}^{2}$'+\
          r' and $\bm{e}_{m}^2$:')
    for Na in Cases:
        for Nb in Cases:
            for Nm in Cases:
                Ba_sq = -Na*Nm
                Bb_sq = -Nb*Nm
                if Ba_sq < 0:
                    Ca_th = cos(th_a)
                    Sa_th = sin(th_a)
                else:
                    Ca_th = cosh(th_a)
                    Sa_th = sinh(th_a)
                if Bb_sq < 0:
                    Cb_th = cos(th_b)
                    Sb_th = sin(th_b)
                else:
                    Cb_th = cosh(th_b)
                    Sb_th = sinh(th_b)
                print(r'%\left [ \bm{e}_{a}^{2},\bm{e}_{b}^{2},\bm{e}_{m}^2\right ] =',\
                      [Na,Nb,Nm])
                Dsq_tmp = Dsq.subs({ca:Ca_th,sa:Sa_th,cb:Cb_th,sb:Sb_th,na:Na,nb:Nb,nm:Nm})
                print(r'%S^{2}-\bm{B}^{2} =',Dsq_tmp,' =',trigsimp(Dsq_tmp))
    print(r'#Thus we have shown that $R_{a}R_{b} = S+\bm{D} = e^{\bm{C}}$ where $\bm{C}$'+\
          r' is a bivector blade.')
    return

def main():
    Format()
    Product_of_Rotors()
    # xpdf(paper=(8.5,11))
    xpdf(pdfprog=None, paper=(8.5,11))
    return

if __name__ == "__main__":
    main()

