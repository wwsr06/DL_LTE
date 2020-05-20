#coding=utf-8
import sys

#print (__name__)

#get %length% samples of signal from iqdata file %fname%
def get_sig_stream(fname,startpos,length):
	iqbuf = [0 for i in range(307200)];
	
	f = open(fname,'r')
	r = f.readlines()
	f.close()
	
	cnt = 0;
	for line in r:
			iqpare = line.split('\t')
			i_p = int(iqpare[0])
			q_p = int(iqpare[1])
			
			iqbuf[cnt] = complex(i_p,q_p)
			#iqbuf[cnt] = i_p
			
			cnt += 1
			if cnt >= length:
				break

	return iqbuf


