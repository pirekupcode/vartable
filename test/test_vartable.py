from nose.tools import eq_, ok_
from vartable.types import CDS, CompoundCDS, SimpleCDS, Result, simple_translate, compound_translate, translate_all
# might need the below
# from pyfakefs.fake_filesystem_unittest import Patcher
# from unittest.mock import patch



def func_blah_test1() -> None:
    ok_(True, "True is true")
    actual, expected = 3, 5
    eq_(3, 5, "Actual {0} is not Expected {1}".format(actual, expected))
    simple = SimpleCDS(True, 0, 3)
       
    # duplicate dirs are impossible anyway, but still using Counter
def func_blah_test2() -> None:
    pass

