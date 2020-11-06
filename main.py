#coding=utf-8
import sys
import dfe
import sync
import sigstream
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft,ifft

#This is master branch

#define global variables
IQDateBuf = [0 for i in range(307200)];
SampleRate = 7680 #15360 30720

print ('LTE_DL ：Starting RX Processing......')

################################################################
#             read IQ Data
################################################################
print ('LTE_DL ：Get IQ stream from file')
#IQDateBuf = sigstream.get_sig_stream("iq.dat",0,76800)
#IQDateBuf = sigstream.get_sig_stream("arfcn1300_cellid408_7p68M_76800samples.txt",0,76800)
#IQDateBuf = sigstream.get_sig_stream("arfcn1300_cellid408_15p36M_153600samples.txt",0,153600)
IQDateBuf = sigstream.get_sig_stream("arfcn1300_cellid408_30p72M_307200samples.txt",0,307200)


################################################################
#             DFE process
################################################################
print ('LTE_DL ：RX DFE processing')
print ('LTE_DL ：Decimation 30.72M->15.36M')
#30.72-->15.36
filterbuf0 = dfe.fir_filter_29taps(IQDateBuf)
filterbuf1 = filterbuf0[0::2]

print ('LTE_DL ：Decimation 15.6M->7.68M')
#15.26-->7.68	
filterbuf2 = dfe.fir_filter_29taps(filterbuf1)
filterbuf3 = filterbuf2[0::2]

print ('LTE_DL ：Decimation 7.68M->3.84M')
#7.68-->3.84
filterbuf4 = dfe.fir_filter_15taps(filterbuf3)
filterbuf5 = filterbuf4[0::2]

print ('LTE_DL ：Decimation 3.84M->1.92M')	
#3.84-->1.92
filterbuf6 = dfe.fir_filter_15tapsA(filterbuf5)
filterbuf7 = filterbuf6[0::2]

print (len(filterbuf7))

'''
x=np.linspace(0,1,4096)
spos = 0
y = filterbuf7[spos:spos+4096]

ff = fft(y)
ffabs = abs(ff)
ffabsnorm = ffabs/len(x)

	
plt.plot(ffabsnorm,label="fft res")
plt.show()
input()
'''


################################################################
#             Sync process
################################################################
print ('LTE_DL ：Coarse Time Sync processing')
#----Coarse time Sync use pss
sync_1p4m_buf = filterbuf7
CorRes = sync.Coarse_Sync(sync_1p4m_buf,1920*10,128+9,1920*5) #1.92
#CorRes = sync.Coarse_Sync(IQDateBuf,3840*10,256+18,3840*5) #3.84
#CorRes = sync.Coarse_Sync(IQDateBuf,7680*10,512+36,7680*5) #7.68
#CorRes = sync.Coarse_Sync(IQDateBuf,15360*10,1024+72,15360*5) #15.36
#CorRes = sync.Coarse_Sync(IQDateBuf,30720*10,2048+144,30720*5) #30.72
CoarseSyncPos = CorRes.index(max(CorRes))
print ('LTE_DL ：Coarse Sync Position Found:',CoarseSyncPos)

'''
plt.plot(CorRes,label="corr res")
plt.legend()
plt.show()
input()
'''

#----Detect NID2 use pss
CoarseSyncBuf = sync_1p4m_buf[CoarseSyncPos+9:CoarseSyncPos+9+128]
NID_2 = sync.NID2_detection(CoarseSyncBuf,128)


#---Detect NID1 use sss
CoarseSyncBuf = IQDateBuf[CoarseSyncPos-128:CoarseSyncPos]
NID_1 = sync.NID1_detection(CoarseSyncBuf,128)


