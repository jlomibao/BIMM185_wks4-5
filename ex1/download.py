#File Name: download.py
#Author: John Francis Lomibao
#PID: A11591509

import urllib2
import sys

try:

	'''
	the first argument should be the file you want to download,
	and the second argument (optional) is what you want to name the file
	'''
	url = sys.argv[1]
	
	if len(sys.argv) == 3:
		fileName = sys.argv[2]
	elif len(sys.argv) == 2:
		fileName = url.split('/')[-1]

	fileToDownload = urllib2.urlopen(url)
	with open(fileName,'wb') as output:
		output.write(fileToDownload.read())

#print usage statement if error occurs
except Exception as e:
	print 'Usage: python download.py <url_to_download> <optional_name_to_give_to_file>'
	print str(e)+' at line {}'.format(sys.exc_info()[-1].tb_lineno)
	



