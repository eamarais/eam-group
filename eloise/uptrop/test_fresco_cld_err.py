from fresco_cld_err import process_file, CloudVariableStore, MAX_LAT, MAX_LON, MIN_LAT, MIN_LON, DELTA_LAT, DELTA_LON
import numpy as np

# I don't like having these here, but they're stuck until I refactor the last bit of fresco_cld_err
out_lon = np.arange(MIN_LON, MAX_LON, DELTA_LON)
out_lat = np.arange(MIN_LAT, MAX_LAT, DELTA_LAT)


def test_process_file():
    X, Y = np.meshgrid(out_lon, out_lat, indexing='ij')
    test_container = CloudVariableStore(X.shape)
    fresco_path = "test_data/S5P_OFFL_L2__NO2____20191025T085808_20191025T103937_10527_01_010302_20191031T111532.nc"
    dlr_path = "test_data/S5P_OFFL_L2__CLOUD__20191025T085808_20191025T103937_10527_01_010107_20191031T081901.nc"
    process_file(dlr_path, fresco_path, test_container)
    nobs_dlr, nobs_fresco = test_container.ge
    assert test_container.nobs_dlr != 0
    assert test_container.nobs_fresco != 0
