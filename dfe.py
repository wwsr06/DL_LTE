#coding=utf-8
import sys
import numpy as np

#print (__name__)


Coef_29Taps = [10,12,0,-16,-20,0,24,32,0,-52,-68,0,134,288,336,288,134,0,-68,-52,0,32,24,0,-20,-16,0,12,10];
def fir_filter_29taps(insig):
	tbuf = [0 for i in range(len(insig)+29)];
	tbuf[0:-29]=insig
	
	
	
	for n in range(1,len(insig)):
		tmpsum = 0
		for i in range(1,29):
			mm = tbuf[n+i]*Coef_29Taps[i]/1024
			tmpsum += mm
		tbuf[n]=tmpsum
	
	filtersig=tbuf[0:-29]
	#print ("Done")
	
	return filtersig


Coef_15Taps = [16,36,0,-68,-64,88,304,400,304,88,-64,-68,0,36,16];
def fir_filter_15taps(insig):
	tbuf = [0 for i in range(len(insig)+15)];
	tbuf[0:-15]=insig
	
	
	
	for n in range(1,len(insig)):
		tmpsum = 0
		for i in range(1,15):
			mm = tbuf[n+i]*Coef_15Taps[i]/1024
			tmpsum += mm
		tbuf[n]=tmpsum
	
	filtersig=tbuf[0:-15]
	#print ("Done")
	
	return filtersig


Coef_15TapsA = [28,4,-44,-64,-4,136,280,352,280,136,-4,-64,-44,4,28];
def fir_filter_15tapsA(insig):
	tbuf = [0 for i in range(len(insig)+15)];
	tbuf[0:-15]=insig
	
	
	
	for n in range(1,len(insig)):
		tmpsum = 0
		for i in range(1,15):
			mm = tbuf[n+i]*Coef_15TapsA[i]/1024
			tmpsum += mm
		tbuf[n]=tmpsum
	
	filtersig=tbuf[0:-15]
	#print ("Done")
	
	return filtersig