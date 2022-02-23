from sets import Set
from math import *
from random import *
import sys
import numpy as np

def test():
  global_arr = [1]
  min_i = 10
  def subtest(n):
    if n == 0:
      return

    if n < min_i:
        min_i = n
    global_arr.append(1)
    print global_arr
    subtest(n - 1)

  subtest(10)


test()
