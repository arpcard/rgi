import os
import sys
import filepaths
import argparse
import json

def main(args):
	if args.in_file == None:
		print("[error] missing input(s)")
		print("[info] Try: python formatJson.py -h")
		exit()
		
	os.system("cat {} | python -m json.tool > {}.json".format(args.in_file ,args.out_file))

def run():
	parser = argparse.ArgumentParser(description='Convert RGI JSON file to Readable JSON file')
	parser.add_argument('-i','--in_file',help='input file must be in JSON format e.g Report.json', required=True)
	parser.add_argument('-o', '--out_file',  dest="out_file", default="ReportFormatted", help="Output JSON file (default=ReportFormatted)")	
	args = parser.parse_args()
	main(args)		

if __name__ == '__main__':
	run()
