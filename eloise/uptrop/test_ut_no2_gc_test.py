from uptrop import ut_no2_gc_test

REGION = "EU"
STR_RES = "8x10"

def test_process_file():
    file_path = "test_data/ts_12_15.EU.20160603.nc"
    rolling_total = ut_no2_gc_test.ProcessedData(REGION, STR_RES)
    rolling_total.process_geochem_day(file_path)
    assert rolling_total.cloud_slice_count >0