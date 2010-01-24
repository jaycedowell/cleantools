#!/usr/bin/python

import os
import sys


fh = open('fourier.txt', 'r')
lines = fh.readlines()
real = lines[::2]
comp = lines[1::2]

currFeed = 0
for pair in zip(real,comp):
	rf = pair[0].split()
	cf = pair[1].split()

	if rf[0].find("FEED") != -1:
		currFeed = int(cf[0])
		print "FEED: %i" % currFeed
		continue

	sys.stdout.write("%s%i = [" % (rf[1], currFeed))
	coeffs = zip(rf[2:], cf[2:])
	
	numb = coeffs[0]
	sys.stdout.write("complex(%.4f, %.4f), $\n" % (float(numb[0]), float(numb[1])))
	
	for numb in coeffs[:-2]:
		sys.stdout.write("\tcomplex(%.4f, %.4f), $\n" % (float(numb[0]), float(numb[1])))
	
	numb = coeffs[-1]
	sys.stdout.write("\tcomplex(%.4f, %.4f)]\n" % (float(numb[0]), float(numb[1])))


