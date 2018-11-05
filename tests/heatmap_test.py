import pytest
import os
from app.MainBase import MainBase

"""
checks that heatmaps are created
"""

inputs = "inputs/heatmap_inputs/"
outputs = "outputs/"
working_directory = os.getcwd()
input_directory = os.path.join(working_directory, inputs)

jsons = []
files = os.listdir(input_directory)
for thing in files:
    file_path = os.path.join(input_directory, thing)
    if thing.endswith(".json") and os.path.isfile(file_path): # Check if it's a file
        jsons.append(thing)

count = len(jsons)

@pytest.fixture
def heatmap():
    return MainBase(api=True)

def validate_heatmap(output_file):
    f = output_file
    if os.path.isfile(f):
        filesize = os.stat(f).st_size
        print("filesize: {}".format(filesize))
        if "test_base-{}.png".format(count) in f:
            if filesize >= 95886:
                return True
            else:
                return False
        elif "test_category-{}.png".format(count) in f:
            if filesize >= 176612:
                return True
            else:
                return False
        elif "test_cluster-{}.png".format(count) in f:
            if filesize >= 95089:
                return True
            else:
                return False
        elif "test_frequency-{}.png".format(count) in f:
            if filesize >= 91316:
                return True
            else:
                return False
        elif "test_category_frequency-{}.png".format(count) in f:
            if filesize >= 248000:
                return True
            else:
                return False
        elif "test_category_cluster-{}.png".format(count) in f:
            if filesize >= 127000:
                return True
            else:
                return False
    else:
        print("missing file: {}".format(f))
        return False

def make_base_heatmap(heatmap, input_directory):
    parser = heatmap.heatmap_args()
    heatmap.heatmap_run(parser.parse_args([
        '--input', input_directory,
        '--output', 'test_base'
    ]))

def make_category_heatmap(heatmap, input_directory, category):
    parser = heatmap.heatmap_args()
    heatmap.heatmap_run(parser.parse_args([
        '--input', input_directory,
        '--output', 'test_category',
        '--category', category
    ]))

def make_cluster_heatmap(heatmap, input_directory):
    parser = heatmap.heatmap_args()
    heatmap.heatmap_run(parser.parse_args([
        '--input', input_directory,
        '--output', 'test_cluster',
        '--cluster', 'both'
    ]))

def make_frequency_heatmap(heatmap, input_directory):
    parser = heatmap.heatmap_args()
    heatmap.heatmap_run(parser.parse_args([
        '--input', input_directory,
        '--output', 'test_frequency',
        '--frequency'
    ]))

def make_category_frequency_heatmap(heatmap, input_directory, category):
    parser = heatmap.heatmap_args()
    heatmap.heatmap_run(parser.parse_args([
        '--input', input_directory,
        '--output', 'test_category_frequency',
        '--category', category,
        '--frequency'
    ]))

def make_category_cluster_heatmap(heatmap, input_directory, category):
    parser = heatmap.heatmap_args()
    heatmap.heatmap_run(parser.parse_args([
        '--input', input_directory,
        '--output', 'test_category_cluster',
        '--category', category,
        '--cluster', 'samples'
    ]))

def _test_base_heatmap(heatmap):
    make_base_heatmap(heatmap, input_directory)
    output_png = os.path.join(working_directory,outputs,"test_base-{}.png".format(count))
    output_eps = os.path.join(working_directory,outputs,"test_base-{}.eps".format(count))
    os.rename(os.path.join(working_directory,"test_base-{}.png".format(count)), output_png)
    os.rename(os.path.join(working_directory,"test_base-{}.eps".format(count)), output_eps)
    assert validate_heatmap(output_png) == True

def _test_category_heatmap(heatmap):
    make_category_heatmap(heatmap, input_directory, "drug_class")
    output_png = os.path.join(working_directory,outputs,"test_category-{}.png".format(count))
    output_eps = os.path.join(working_directory,outputs,"test_category-{}.eps".format(count))
    os.rename(os.path.join(working_directory,"test_category-{}.png".format(count)), output_png)
    os.rename(os.path.join(working_directory,"test_category-{}.eps".format(count)), output_eps)
    assert validate_heatmap(output_png) == True

def _test_cluster_heatmap(heatmap):
    make_cluster_heatmap(heatmap, input_directory)
    output_png = os.path.join(working_directory,outputs,"test_cluster-{}.png".format(count))
    output_eps = os.path.join(working_directory,outputs,"test_cluster-{}.eps".format(count))
    os.rename(os.path.join(working_directory,"test_cluster-{}.png".format(count)), output_png)
    os.rename(os.path.join(working_directory,"test_cluster-{}.eps".format(count)), output_eps)
    assert validate_heatmap(output_png) == True

def _test_frequency_heatmap(heatmap):
    make_frequency_heatmap(heatmap, input_directory)
    output_png = os.path.join(working_directory,outputs,"test_frequency-{}.png".format(count))
    output_eps = os.path.join(working_directory,outputs,"test_frequency-{}.eps".format(count))
    os.rename(os.path.join(working_directory,"test_frequency-{}.png".format(count)), output_png)
    os.rename(os.path.join(working_directory,"test_frequency-{}.eps".format(count)), output_eps)
    os.rename(os.path.join(working_directory,"test_frequency{}-frequency.txt".format(count)),
        os.path.join(working_directory,outputs,"test_frequency{}-frequency.txt".format(count)))
    assert validate_heatmap(output_png) == True

def _test_category_frequency_heatmap(heatmap):
    make_category_frequency_heatmap(heatmap, input_directory, "gene_family")
    output_png = os.path.join(working_directory,outputs,"test_category_frequency-{}.png".format(count))
    output_eps = os.path.join(working_directory,outputs,"test_category_frequency-{}.eps".format(count))
    os.rename(os.path.join(working_directory,"test_category_frequency-{}.png".format(count)), output_png)
    os.rename(os.path.join(working_directory,"test_category_frequency-{}.eps".format(count)), output_eps)
    os.rename(os.path.join(working_directory,"test_category_frequency{}-frequency.txt".format(count)),
        os.path.join(working_directory,outputs,"test_category_frequency{}-frequency.txt".format(count)))
    assert validate_heatmap(output_png) == True

def _test_category_cluster_heatmap(heatmap):
    make_category_cluster_heatmap(heatmap, input_directory, "resistance_mechanism")
    output_png = os.path.join(working_directory,outputs,"test_category_cluster-{}.png".format(count))
    output_eps = os.path.join(working_directory,outputs,"test_category_cluster-{}.eps".format(count))
    os.rename(os.path.join(working_directory,"test_category_cluster-{}.png".format(count)), output_png)
    os.rename(os.path.join(working_directory,"test_category_cluster-{}.eps".format(count)), output_eps)
    assert validate_heatmap(output_png) == True
