import numpy as np
import scipy

for ix in range(0, 11):
    for iy in range(0, 11):
        z = complex(ix, iy)
        print(scipy.special.gamma(z))
