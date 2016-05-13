import os
import sys
import filepaths
import argparse
import json
#import logging

def main(args):
	if args.in_file == None:
		print "[error] missing input(s)"
		print "[info] Try: python formatJson.py -h"	
		exit()

	infile = args.in_file 
	outfile = args.out_file

	os.system("cat "+infile+" | python -m json.tool > "+outfile+".json" )

def run():
	parser = argparse.ArgumentParser(description='Convert RGI JSON file to Readable JSON file')
	parser.add_argument('-i','--in_file',help='input file must be in JSON format e.g Report.json')
	parser.add_argument('-o', '--out_file',  dest="out_file", default="ReportFormatted", help="Output JSON file (default=ReportFormatted)")	
	args = parser.parse_args()
	#logging.basicConfig(filename='app.log', level=logging.DEBUG, format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')
	main(args)		

if __name__ == '__main__':
	run()
