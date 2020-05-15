#coding=utf-8
import sys
import sync
import sigstream
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft,ifft

#define global variables
IQDateBuf = [0 for i in range(307200)];
SampleRate = 7680 #15360 30720


################################################################
#             read IQ Data
################################################################
#IQDateBuf = sigstream.get_sig_stream("iq.dat",0,76800)
#IQDateBuf = sigstream.get_sig_stream("arfcn1300_cellid408_7p68M_76800samples.txt",0,76800)
#IQDateBuf = sigstream.get_sig_stream("arfcn1300_cellid408_15p36M_153600samples.txt",0,153600)
IQDateBuf = sigstream.get_sig_stream("arfcn1300_cellid408_30p72M_307200samples.txt",0,307200)

'''
x=np.linspace(0,1,4096)
spos = 0
y = IQDateBuf[spos:spos+2048]

ff = fft(y)
ffabs = abs(ff)
ffabsnorm = ffabs/len(x)

fftplot = [0 for i in range(2048)];
fftplot[0:1023] = ffabsnorm[1024:2048]
fftplot[1024:2048] = ffabsnorm[0:1023]
	
plt.plot(ffabsnorm,label="fft res")
plt.show()
input()
'''

################################################################
#             Sync process
################################################################

#----Corse time Sync use pss
#CorRes = sync.Corse_Sync(IQDateBuf,153600,512+36,7680*5)
#CorRes = sync.Corse_Sync(IQDateBuf,7680*10,512+36,7680*5)
#CorRes = sync.Corse_Sync(IQDateBuf,15360*10,1024+72,15360*5)
CorRes = sync.Corse_Sync(IQDateBuf,30720*10,2048+144,30720*5)
CorseSyncPos = CorRes.index(max(CorRes))

'''
plt.plot(CorRes,label="corr res")
plt.legend()
plt.show()
'''

#----Detect NID2 use pss
CorseSyncBuf = IQDateBuf[CorseSyncPos+144:CorseSyncPos+144+2048]
NID_2 = sync.NID2_detection(CorseSyncBuf,2048)

#---Detect NID1 use sss
CorseSyncBuf = IQDateBuf[CorseSyncPos-2048:CorseSyncPos]
NID_1 = sync.NID1_detection(CorseSyncBuf,2048)


