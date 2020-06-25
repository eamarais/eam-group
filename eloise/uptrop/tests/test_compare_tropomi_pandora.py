from uptrop.compare_tropomi_pandora import get_days_since_data_start
import datetime as dt
def test_get_offset_year_day():
    assert get_days_since_data_start(dt.date(year=2019, month=6, day=1)) == 0
    assert get_days_since_data_start(dt.date(year=2019, month=7, day=1)) == 30
    assert get_days_since_data_start(dt.date(year=2020, month=5, day=31)) == 365

