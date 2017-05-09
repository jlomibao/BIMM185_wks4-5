#File Name: read.py
#Author: John Francis Lomibao
#PID: A11591509

import os
import re
import sys

fileToRead = sys.argv[1]

with open(fileToRead, 'r') as f:
	data = [line.strip() for line in f]
for i in data:
	m = re.findall(r"(UP000[0-9]+)", i)
	for j in m:
		print j
	