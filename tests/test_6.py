import json
import os

import pytest

from app.MainBase import MainBase

inputs = "inputs/"
outputs = "outputs/"
alignment_tool = "diamond"
orf_caller = "PYRODIGAL"

working_directory = os.getcwd()


# Run all tests with
# pytest test_6.py -v -rxs --color=auto --durations=0
# or
# pytest test_6.py -v -rxs --color=auto --durations=0 -k "tet"

@pytest.fixture
def rgi():
    return MainBase(api=True)

def run_rgi(rgi, input_type, input_sequence, output_file):
    parser = rgi.main_args()
    rgi.main_run(parser.parse_args([
        '--input_type', input_type,
        '--input_sequence', input_sequence,
        '--output_file', output_file,
        '--alignment_tool', alignment_tool,
        '--orf_finder', orf_caller,
        '--clean',
        '--include_loose',
        '--include_nudge',
        '--low_quality',
        '--debug'
    ]))

def validate_results(filepath, perc_identity=0, ARO_name='', type_match=''):
    pi = ""
    name = ""
    tm = ""
    f = os.path.join("{}".format(filepath))
    if os.path.isfile(f):
        with open(f) as json_file:
            json_data = json.load(json_file)
            for i in json_data:
                if i not in ["_metadata"]:
                    for j in json_data[i]:
                        for k in json_data[i][j]:
                            pi = json_data[str(i)][str(j)]["perc_identity"]
                            name = json_data[str(i)][str(j)]["ARO_name"]
                            tm = json_data[str(i)][str(j)]["type_match"]
                        if pi == perc_identity and name == ARO_name and tm == type_match:
                            # print(pi, name, tm)
                            return True
            return False
    else:
        print("missing file: {}".format(f))
        return False

def test_rgi_nucleotide_sequence_tetM(rgi):
    filename = "Cdifficile_strain_DSM_27543_transposon.fasta"
    output_file = os.path.join(working_directory, outputs, "{}.json".format(filename))
    run_rgi(rgi, 'contig', os.path.join(working_directory, inputs, filename), output_file)

    # correct annotation is tetM
    assert validate_results(output_file, 90.3, 'tet(M)', 'Strict')




