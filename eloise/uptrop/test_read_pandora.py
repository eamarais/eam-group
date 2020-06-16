import read_pandora
import os

test_file_path = "test_data/Pandora101s1_Izana_L2Tot_rnvs1p1-7.txt"

def test_readpandora():
    loc, out_df = read_pandora.readpandora(test_file_path)

def test_get_column_indicies():
    groups = read_pandora.get_column_description_index(test_file_path)
    assert len(groups) == 37

def test_get_start_of_data():
    start = read_pandora.get_start_of_data(test_file_path)
    assert start == 59

