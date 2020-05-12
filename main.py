#coding=utf-8
import sys
import sync
import downlink
import matplotlib.pyplot as plt
import numpy as np

#define global buffer
IQDateBuf = [ [0] * 2 for i in range(163840)];


################################################################
#             read IQ Data
################################################################
f = open('iq.dat','r')
r = f.readlines()
f.close()

cnt = 0;
for line in r:
		iqpare = line.split('\t')
		i_p = int(iqpare[0])
		q_p = int(iqpare[1])
		
		IQDateBuf[cnt][0] = i_p
		IQDateBuf[cnt][1] = q_p
		
		cnt += 1
		if cnt > 163839:
			break


################################################################
#             Sync process
################################################################

NewIQDateBuf = [ [0] * 2 for i in range(163840)];
for i in range(76800):
	NewIQDateBuf[i] = IQDateBuf[76800+i]



CorRes = sync.corr_PSSSSS(IQDateBuf,153600,512+36,7680*5)
plt.plot(CorRes,label="corr res")
plt.legend()
plt.show()


'''
plt.plot(x,label="real")
plt.legend()
plt.show()
plt.plot(y,label="imag")
plt.legend()
plt.show()
'''


'''
x = np.linspace(0.05, 10, 1000)
y = np.cos(x)
#plt.plot(x, y, ls="-", lw=2, label="plot figure")
plt.plot(y)
plt.legend()
plt.show()

print (__name__)

import sync
print (sync.fun_a(10))

import downlink
print (downlink.fun_a(10))		
'''
