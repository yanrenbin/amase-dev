import numpy  as np
import matplotlib.pyplot as plt

lambda0 = 4970.
lamarr = np.arange(900)+6000.
lamarr = np.arange(600)+4600.
resolution0 = 15000.
fcoll = 200.
fiber = 0.05
alpha = np.arctan(resolution0*fiber/(2*fcoll))

nu = 2*np.sin(alpha)/lambda0

betaarr = np.arcsin(nu*lamarr-np.sin(alpha))
res = lamarr*2*fcoll*np.sin(alpha)/(fiber*lambda0*np.cos(betaarr))

plt.plot(lamarr,res)
plt.show()
