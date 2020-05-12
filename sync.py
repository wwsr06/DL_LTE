#coding=utf-8
import sys

print (__name__)

def corr_pow(A,B):
	
	a0 = A[0]
	b0 = A[1]
	a1 = B[0]
	b1 = 0-B[1]
	
	rp = a0*a1+b0*b1
	ip = a0*b1+b0*a1
	
	res = rp*rp+ip*ip
	
	return res
	

def poscalc(pos):
	if pos >= 153600:
		pos = pos - 153600
	
	return pos


def corr_PSSSSS(sig,slidelen,winlen,winitl):
	ResCorr = [0 for i in range(slidelen)]

	#calc initial value
	sss_offset = 1648
	i = 0
	ResCorr[i]=0
	for j in range(winlen):
		ResCorr[i] += corr_pow(sig[i+j],sig[winitl+i+j]) #pss part
		#ResCorr[i] += corr_pow(sig[i+j+sss_offset],sig[winitl+i+j+sss_offset]) #sss part
	
	
	#slide correlation	
	for scnt in range(1,slidelen):
		#pss part
		subv_pss = corr_pow(sig[poscalc(scnt-1)],sig[poscalc(scnt-1+winitl)])
		addv_pss = corr_pow(sig[poscalc(scnt+winlen)],sig[poscalc(scnt+winitl+winlen)])
		
		#sss part
		subv_sss = corr_pow(sig[poscalc(scnt-1+sss_offset)],sig[poscalc(scnt-1+winitl+sss_offset)])
		addv_sss = corr_pow(sig[poscalc(scnt+winlen+sss_offset)],sig[poscalc(scnt+winitl+winlen+sss_offset)])
		
		ResCorr[scnt] = ResCorr[scnt-1]+addv_pss-subv_pss  #+addv_sss-subv_sss
		
	#scaleing
	for i in range(slidelen):
		ResCorr[i] >>= 32
	
	return ResCorr
	


