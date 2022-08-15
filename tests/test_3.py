import pytest
import os, json
from app.MainBase import MainBase

inputs = "inputs/"
outputs = "outputs/"
alignment_tool = "diamond"
working_directory = os.getcwd()

# Run all tests with
# pytest test_3.py -v -rxs --color=auto --durations=0
# or
# pytest test_3.py -v -rxs --color=auto --durations=0 -k "create"

@pytest.fixture
def rgi():
	return MainBase(api=True)

def test_create_local_db(rgi):
	parser = rgi.load_args()
	f = os.path.join(working_directory,inputs,"{}".format("card.json"))
	rgi.load_run(parser.parse_args([
		'--card_json', f,
		'--local',
		'--debug'
	]))

	assert (os.path.isfile(f) and os.path.exists(f)), \
	print("add card.json to {} directory and re-run test".format(os.path.join(working_directory,inputs)))
