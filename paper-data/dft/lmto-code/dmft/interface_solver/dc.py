#!/bin/env python
import sys
U = float(sys.argv[1])
J = float(sys.argv[2])
n = float(sys.argv[3])

DC = U*(n-0.5)-J*(n*0.5-0.5)
print(DC)
