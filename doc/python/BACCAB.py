from galgebra.printer import Format, xpdf
from galgebra.ga import Ga

Format()
g4d = Ga('a b c d')
(a, b, c, d) = g4d.mv()

print('g_{ij} =', g4d.g)
print('\\bm{a|(b*c)} =', a | (b * c))
print('\\bm{a|(b^c)} =', a | (b ^ c))
print('\\bm{a|(b^c^d)} =', a | (b ^ c ^ d))
# FIXME:FIXED this should print 0, got blank
print('\\bm{a|(b^c)+c|(a^b)+b|(c^a)} =',
      (a | (b ^ c)) + (c | (a ^ b)) + (b | (c ^ a)))
print('\\bm{a*(b^c)-b*(a^c)+c*(a^b)} =',
      a * (b ^ c) - b * (a ^ c) + c * (a ^ b))
print('\\bm{a*(b^c^d)-b*(a^c^d)+c*(a^b^d)-d*(a^b^c)} =',
      a * (b ^ c ^ d) - b * (a ^ c ^ d) + c * (a ^ b ^ d) - d * (a ^ b ^ c))
print('\\bm{(a^b)|(c^d)} =', (a ^ b) | (c ^ d))
print('\\bm{((a^b)|c)|d} =', ((a ^ b) | c) | d)
print('\\bm{(a^b)\\times (c^d)} =', Ga.com(a ^ b, c ^ d))
xpdf(paper='letter', prog=True)
