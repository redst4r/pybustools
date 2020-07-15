import pyximport; pyximport.install()
# import iox
import iotest
import iotest_pyx
import gmpy2
the_int = 1944634900  # 16BP
%timeit iotest.int2base(the_int,4)
%timeit iotest_pyx.int2base(the_int,4)
%timeit gmpy2.digits(the_int, 4)
%timeit iox.int2seq(344634)
