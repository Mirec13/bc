#import numpy

#ud = numpy.random.uniform(0, 1, 1);
#ud *= 20;
#print(ud);

#from scipy.stats import levy
#num=levy.rvs(0, 100)
#print(levy.cdf(num, 0, 100));
#print(num);

from mpmath import gamma
from math import sin
print(float(gamma(0.00245)))
alpha = sin((180 * 1.3) / 2)
print(alpha)
#from scipy.stats import levy_stable
#alpha, beta = 1, 1
#r = levy_stable.rvs(alpha, beta, size=10)
#print(r);

from collections import namedtuple
MyStruct = namedtuple("MyStruct", "field1 field2 field3")
m = MyStruct(field1="foo", field2="bar", field3="baz")
print(m.field1)