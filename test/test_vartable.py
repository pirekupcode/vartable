from nose.tools import eq_, ok_

# might need the below
# from pyfakefs.fake_filesystem_unittest import Patcher
# from unittest.mock import patch



def func_blah_test1():
    ok_(True, "True is true")
    actual, expected = 3, 5
    eq_(3, 5, "Actual {0} is not Expected {1}".format(actual, expected))
       
    # duplicate dirs are impossible anyway, but still using Counter
def func_blah_test2():
    pass
