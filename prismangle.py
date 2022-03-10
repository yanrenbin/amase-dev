import math
import numpy as np


alpha = np.arange(150)/10.+45.
refindex = 1.4563

beta = np.arcsin(np.sin(alpha*math.pi/180.)/refindex)*180./math.pi

diff = alpha-beta
angle = np.interp(20,diff,alpha)
