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

@pytest.fixture
def heatmap():
    return MainBase(api=True)

def validate_heatmap(output_file):
    f = output_file
    if os.path.isfile(f):
        if os.path.getsize(f) > 100000:
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

def test_base_heatmap(heatmap):
    make_base_heatmap(heatmap, input_directory)
    output_png = os.path.join(working_directory,outputs,"test_base-4.png")
    output_eps = os.path.join(working_directory,outputs,"test_base-4.eps")
    os.rename(os.path.join(working_directory,"test_base-4.png"), output_png)
    os.rename(os.path.join(working_directory,"test_base-4.eps"), output_eps)
    assert validate_heatmap(output_png) == True

def test_category_heatmap(heatmap):
    make_category_heatmap(heatmap, input_directory, "drug_class")
    output_png = os.path.join(working_directory,outputs,"test_category-4.png")
    output_eps = os.path.join(working_directory,outputs,"test_category-4.eps")
    os.rename(os.path.join(working_directory,"test_category-4.png"), output_png)
    os.rename(os.path.join(working_directory,"test_category-4.eps"), output_eps)
    assert validate_heatmap(output_png) == True

def test_cluster_heatmap(heatmap):
    make_cluster_heatmap(heatmap, input_directory)
    output_png = os.path.join(working_directory,outputs,"test_cluster-4.png")
    output_eps = os.path.join(working_directory,outputs,"test_cluster-4.eps")
    os.rename(os.path.join(working_directory,"test_cluster-4.png"), output_png)
    os.rename(os.path.join(working_directory,"test_cluster-4.eps"), output_eps)
    assert validate_heatmap(output_png) == True

def test_frequency_heatmap(heatmap):
    make_frequency_heatmap(heatmap, input_directory)
    output_png = os.path.join(working_directory,outputs,"test_frequency-4.png")
    output_eps = os.path.join(working_directory,outputs,"test_frequency-4.eps")
    os.rename(os.path.join(working_directory,"test_frequency-4.png"), output_png)
    os.rename(os.path.join(working_directory,"test_frequency-4.eps"), output_eps)
    os.rename(os.path.join(working_directory,"test_frequency4-frequency.txt"),
        os.path.join(working_directory,outputs,"test_frequency4-frequency.txt"))
    assert validate_heatmap(output_png) == True

def test_category_frequency_heatmap(heatmap):
    make_category_frequency_heatmap(heatmap, input_directory, "gene_family")
    output_png = os.path.join(working_directory,outputs,"test_category_frequency-4.png")
    output_eps = os.path.join(working_directory,outputs,"test_category_frequency-4.eps")
    os.rename(os.path.join(working_directory,"test_category_frequency-4.png"), output_png)
    os.rename(os.path.join(working_directory,"test_category_frequency-4.eps"), output_eps)
    os.rename(os.path.join(working_directory,"test_category_frequency4-frequency.txt"),
        os.path.join(working_directory,outputs,"test_category_frequency4-frequency.txt"))
    assert validate_heatmap(output_png) == True

def test_category_cluster_heatmap(heatmap):
    make_category_cluster_heatmap(heatmap, input_directory, "resistance_mechanism")
    output_png = os.path.join(working_directory,outputs,"test_category_cluster-4.png")
    output_eps = os.path.join(working_directory,outputs,"test_category_cluster-4.eps")
    os.rename(os.path.join(working_directory,"test_category_cluster-4.png"), output_png)
    os.rename(os.path.join(working_directory,"test_category_cluster-4.eps"), output_eps)
    assert validate_heatmap(output_png) == True
