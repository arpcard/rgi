""" Transform RGI json results to d3 format that we can plot """
"""
{
	name: "ARO",
	children: [
	{
		name: "Perfect",
		children: []
	},
	{
		name: "Strict",
		children: [ ]
	}
	]
}

"""
"""
import sys
import os
import json

def main(input_file, output_file):
	with open("{}".format(input_file), 'r') as jfile:
		data = json.load(jfile)

	# write d3 json
	for i in data:
		print(data[i])

if __name__ == "__main__":
	# Report.json, Out.json
    main(sys.argv[1],sys.argv[1])

"""

"""
>seq_a
GAGAGATTTTCCAATTCGACG-------CGGGGTCAGG--GAAATTT
>seq_b
GAGAGATTGGCCTTAACTACCCAACCCACGGCCTGACCGAGGTCTTC
"""
"""
A64G
A76G
G77A
C88T
T89C
T101C
C122A
A127G
A129T
T135C
G139A
A165G
C181T
G185T
C186A
A187C
A188G
A197C
G198A
A203C
C204T
C212A
T231A
G232A
G237T
C264A
C285T
C307T
A320G
T333C
C422A
A462G
T472C
C476A
T592C
T593A
G594A
T612G
G618A
C622T
C626T
C633T
T643A
C649T
A650T
A651G
T655A
C662T
G684T
C710G
G747A
C813T
C828T
G831A
C841G
C842T
G848A
G849C
C859T
G876A
G989A
G990C
C1003G
G1004C
G1008T
T1009C
T1010C
C1013T
A1024G
C1041G
G1042C
G1137C
C1140A
C1141A
G1142T
C1144G
C1145T
A1167G
T1177C
C1221G
C1222T
G1248A
C1249G
C1285A
G1296C
C1297T
G1359A
C1371T
"""
# 2 -> 1542
# 86 -> 1627
def snps():
	from Bio import AlignIO
	y=0

	alignment = AlignIO.read("fasta.fas", "fasta")
	print(alignment)
	exit()
	for r in range(0,len(alignment[1].seq)):
	    if alignment[0,r] != alignment[1,r]:
	        if alignment[0,r] != "-" and alignment[1,r] != "-":
	            y=y+1
	            # print r, alignment[0,r], alignment[1,r] #, y
	            # print alignment[0,r], r ,alignment[1,r], ": ", alignment[0].id, " : ", alignment[1].id
	            print("{}{}{}".format(alignment[0,r], r, alignment[1,r]))
	        else:
	            y=0

snps()

