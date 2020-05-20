#coding=utf-8
import sys
import numpy as np
import math
from scipy.fftpack import fft,ifft
import matplotlib.pyplot as plt

#print (__name__)

def corr_pow(A,B):
	
	a0 = A.real
	b0 = A.imag
	a1 = B.real
	b1 = 0-B.imag
	
	rp = a0*a1+b0*b1
	ip = a0*b1+b0*a1
	
	res = (rp*rp+ip*ip)
	
	return res
	
def freqerr_est(tbuf0,tbuf1):
	
	fsig0 = fft(tbuf0)
	fsig1 = fft(tbuf1)
	c62_sc_0 =  np.hstack((fsig0[-31:],fsig0[:31]))
	c62_sc_1 =  np.hstack((fsig1[-31:],fsig1[:31]))
		
	cores = 0	
	for i in range(0,62):
		cores += c62_sc_0[i] * c62_sc_1[i].conjugate()
	
	at = np.arctan(cores.imag/cores.real)
	print (at/3.1415926)
	input()
	
	return 1


def poscalc(slidelen,pos):
	if pos >= slidelen:
		pos = pos - slidelen
	
	return pos


def Coarse_Sync(sig,slidelen,winlen,winitl):
	ResCorr = [0 for i in range(slidelen)]

	#calc initial value
	i = 0
	ResCorr[i]=0
	for j in range(winlen):
		ResCorr[i] += corr_pow(sig[i+j],sig[winitl+i+j]) #pss part
	
	
	#slide correlation	
	for scnt in range(1,slidelen):
		#pss part
		subv_pss = corr_pow(sig[poscalc(slidelen,scnt-1)],sig[poscalc(slidelen,scnt-1+winitl)])
		addv_pss = corr_pow(sig[poscalc(slidelen,scnt+winlen)],sig[poscalc(slidelen,scnt+winitl+winlen)])
		
		
		ResCorr[scnt] = ResCorr[scnt-1]+addv_pss-subv_pss  
		
	#scaleing
	for i in range(slidelen):
		ResCorr[i] /= 102400
	
	return ResCorr

def regen_pss(nid2):
	root_index_buf = [25,29,34]
	u = root_index_buf[nid2]
	pssseq = [0 for i in range(62)]
	
	pi = 3.1415926
	for n in range(0,62):		
		if n<=30:
			rp = np.cos(pi*u*n*(n+1)/63)
			ip = 0-np.sin(pi*u*n*(n+1)/63)
			d = complex(rp,ip)
		else :
			rp = np.cos(pi*u*(n+1)*(n+2)/63)
			ip = 0-np.sin(pi*u*(n+1)*(n+2)/63)
			d = complex(rp,ip)
		
		pssseq[n]=d

	return pssseq

def nid2_corr(rxfsig,nid2):
	#gen local pss seq
	localpss = [0 for i in range(62)]
	orgpss = regen_pss(nid2)
	
	#gen multi-bin
	shift_num = 0
	for i in range(0,62):
		pos = i+shift_num
		if pos>=62:
			pos -= 62
		if pos<0:
			pos += 62
		localpss[i] = orgpss[pos]
	
	
	
	#gen diff signal for correlation
	rxfsig_diff = [0 for i in range(62)]
	localpss_diff = [0 for i in range(62)]
	for i in range(0,61):
		rxfsig_diff[i] = rxfsig[i+1]/rxfsig[i]
		localpss_diff[i] = localpss[i+1]/localpss[i]
	
	ResCorr=0
	for i in range(62):		
		ResCorr += rxfsig_diff[i] * localpss_diff[i].conjugate()
		
	return abs(ResCorr)

def NID2_detection(tsig,fftN):
	
	#---FFT of symble include pss signal
	fsig = fft(tsig)
	ffabs = abs(fsig)
	fsig = fsig/max(ffabs)
	'''
	plt.plot(ffabs,label="fft res")
	plt.legend()
	plt.show()
	'''
	#---Get central 62 subcarriers
	c72_sc =  np.hstack((fsig[1:37] , fsig[-36:]))
	c62_sc = c72_sc[5:-5]
	
		
	#---correlation with local pss seq, and find nid2
	res_cor = [0 for i in range(3)]
	for i in range(3):
		res_cor[i] = nid2_corr(c62_sc,i)
	nid2_dect = res_cor.index(max(res_cor))
	
	
	print (res_cor)
	return nid2_dect


def regen_sss(nid1,nid2,sf):
	sss_seq = [0 for i in range(62)]
	
	q_prime = math.floor(nid1/30)
	q = math.floor(((nid1+q_prime*(q_prime+1)/2))/30)
	m_prime = nid1 + q *(q+1)/2
	m0 = (m_prime % 31)
	m0 = int(m0)
	m1 = (m0 + math.floor(m_prime/31)+1) % (31)
	
	#---generate d_even() sequence
	x_s = [0 for i in range(31)]
	x_s[4] = 1
	for i in range(0,26):
		x_s[i+5] = (x_s[i+2]+x_s[i]) % 2
	
	#Generate the sequence x_c() : x() for calculating c_tilda()
	x_c = [0 for i in range(31)]
	x_c[4] = 1
	for i in range(0,26):
		x_c[(i+5)] = (x_c[i+3]+x_c[i]) % 2
		
	#% Generate the sequence s_tilda()
	s_tilda = [0 for i in range(31)]
	for i in range(0,31):
		s_tilda[i] = 1 - 2*x_s[i]
	
	#% Generate the sequence c_tilda()
	c_tilda = [0 for i in range(31)]
	for i in range(0,31):
		c_tilda[i] = 1 - 2*x_c[i]
	
	#% Generate s0_m0_even()
	s0_m0_even = [0 for i in range(31)]
	for n in range(0,31):
		s0_m0_even[n] = s_tilda[(n+m0) % 31]
		
	#% Generate c0_even()
	c0_even = [0 for i in range(31)]
	for n in range(0,31):
		c0_even[n] = c_tilda[(n+nid2) % 31]
	
	#% Calculate d_even_sub0
	d_even_sub0 = np.multiply(s0_m0_even , c0_even)
	
	#%%%%%%%%%%%%%%%%% generate d_odd() sequence %%%%%%%%%%%%%%%%
	#% Generate the sequence x_s() : x() for calculating s_tilda()
	x_z = [0 for i in range(31)]
	x_z[4] = 1
	for i in range(0,26):
	    x_z[(i+5)] = (x_z[i+4] + x_z[i+2] + x_z[i+1]+ x_z[i]) % 2

	#% Generate the sequence z_tilda()
	z_tilda = [0 for i in range(31)]
	for i in range(0,31):
		z_tilda[i] = 1 - 2*x_z[i]
		
	#% Generate s1_m1_odd()
	s1_m1_odd = [0 for i in range(31)]
	for n in range(0,31):
		s1_m1_odd[n] = s_tilda[(n+m1) % 31]
		
	#% Generate c1_odd()
	c1_odd = [0 for i in range(31)]
	for n in range(0,31):
		c1_odd[n] = c_tilda[(n+nid2+3) % 31]

	#% Generate z1_m0_odd()
	z1_m0_odd = [0 for i in range(31)]
	for n in range(0,31):
		z1_m0_odd[n] = z_tilda[(n+(m0%8)) % 31]

	#% Calculate d_odd_sub0
	d_odd_sub0 = np.multiply(s1_m1_odd , c1_odd) #* z1_m0_odd;
	d_odd_sub0 = np.multiply(d_odd_sub0 , z1_m0_odd)
	

	
	if sf == 'SF0':
		for n in range(0,31):
			sss_seq[2*n] = d_even_sub0[n]
			sss_seq[2*n+1] = d_odd_sub0[n]
	
	if sf == 'SF5':
		for n in range(0,31):
			sss_seq[2*n] = d_odd_sub0[n]
			sss_seq[2*n+1] = d_even_sub0[n]

	
	return sss_seq



def nid1_corr(rxfsig,nid1):
	#gen local sss seq
	localsss = [0 for i in range(62)]
	localsss = regen_sss(nid1,2,'SF5')
	
	#gen diff signal for correlation
	rxfsig_diff = [0 for i in range(62)]
	localsss_diff = [0 for i in range(62)]
	
	for i in range(0,61):
		rxfsig_diff[i] = rxfsig[i+1]/rxfsig[i]
		localsss_diff[i] = localsss[i+1]/localsss[i]
	
	ResCorr=0
	for i in range(61):
		#ResCorr += localsss[i] * rxfsig[i].conjugate()
		ResCorr += localsss_diff[i] * rxfsig_diff[i].conjugate()

	return abs(ResCorr)

def NID1_detection(tsig,fftN):
	
	#---FFT of symbol include sss signal
	fsig = fft(tsig)
	ffabs = abs(fsig)
	fsig = fsig/max(ffabs)
	
	'''
	fout = open('debugout.txt', 'w')
	for i in range(len(fsig)):
		fout.write(str(fsig[i]))
		fout.write('\n')
	fout.close()
	'''
	
	
	#---Get centrul 62 subcarriers
	c72_sc =  np.hstack((fsig[-36:],fsig[1:37]))
	c62_sc = c72_sc[5:-5]
		
	
	#---correlation with local sss seq, and find nid1
	res_cor = [0 for i in range(168)]
	for i in range(0,168):
		res_cor[i] = nid1_corr(c62_sc,i)
	
	print(res_cor.index(max(res_cor)))
	
	plt.plot(res_cor,label="sss corr res")
	plt.legend()
	plt.show()
	input()
