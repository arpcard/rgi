import os
import pytest
from unittest.mock import Mock

from app.RGI import RGI
from app.HomologModel import Homolog
from app.Database import Database
from app.MainBase import MainBase

# pip3 install mock
# Run all tests with
# pytest test_2.py -v -rxs --color=auto --durations=0
# or
# pytest test_2.py -v -rxs --color=auto --durations=0 -k "fail"

@pytest.fixture
def mock_db():
	return Mock(spec=Database)

@pytest.fixture
def mock_rgi_db():
	return Mock(spec=RGI)

@pytest.fixture
def rgi():
	return MainBase(api=True)

inputs = "inputs/"
outputs = "outputs/"
alignment_tool = "diamond"
working_directory = os.getcwd()

# TODO:: finish tests
def _test_fail_not_called(mock_db):
	# print(dir(mock_db))
	mock_rgi_db.build_databases.assert_called_with() # Fail
	assert 0


# TODO:: finish tests
def _test_rgi_to_db(rgi, mock_rgi_db):
	filename = "homolog.fasta"
	parser = rgi.main_args()
	args = parser.parse_args([
		'--input_sequence', os.path.join(working_directory,inputs,filename),
		'--output_file', os.path.join(working_directory,outputs,"{}.json".format(filename)),
		'--alignment_tool', alignment_tool,
		'--clean',
		'--include_loose'
    ])
	rgi_obj = RGI(**vars(args))
	rgi_obj.run()

	# print(dir(rgi_obj))
	# print(dir(mock_rgi_db))

