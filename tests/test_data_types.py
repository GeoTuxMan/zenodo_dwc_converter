import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import data_types as dt  # noqa: E402


def test_classify_by_extension_and_name():
    assert dt.classify_name('run1_R1.fastq.gz') == 'genetic'
    assert dt.classify_name('photos/IMG_001.JPG') == 'imagery'
    assert dt.classify_name('bathymetry_50m.tif') == 'bathymetry'
    assert dt.classify_name('analysis.R') == 'documentation'
    assert dt.classify_name('CPOD_detections.csv') == 'acoustic'
    assert dt.classify_name('model.RData') == 'binary'


def test_classify_tables_by_columns():
    chem = pd.DataFrame({'Station': ['A', 'B', 'C'], 'NO3 (umol/L)': ['1', '2', '3'],
                         'PO4 (umol/L)': ['0.1', '0.2', '0.3'], 'O2 (ml/L)': ['5', '6', '7']})
    assert dt.classify_table('nutrients.csv', chem, '')[0] == 'chemistry'
    phys = pd.DataFrame({'Temperature (°C)': ['1', '2', '3'], 'Salinity': ['35', '35', '35'],
                         'Pressure (dbar)': ['1', '5', '10']})
    assert dt.classify_table('cast.csv', phys, '')[0] == 'physics'
    lookup = pd.DataFrame({'Code': ['X'], 'Name': ['Y']})
    assert dt.classify_table('book.xlsx::Country', lookup, '')[0] == 'documentation'


def test_standardise_adds_iso_dates_and_decimal_degrees():
    df = pd.DataFrame({'Date': ['12/20/2024', '12/19/2024'], 'Lat': ["45°30'N", "45°15'N"],
                       'Lon': ["29°45'E", "29°30'E"], 'T (°C)': ['9', '8'], 'empty': [None, None]})
    out, mapping, added = dt.standardise(df, dt.map_columns(list(df.columns)))
    assert 'empty' not in out.columns
    assert list(out['date_ISO8601']) == ['2024-12-20', '2024-12-19']
    assert list(out['latitude_decimal']) == [45.5, 45.25]
    nd, _, nd_added = dt.standardise(df, dt.map_columns(list(df.columns)), add_columns=False)
    assert not nd_added and 'date_ISO8601' not in nd.columns


def test_export_writes_csv_and_xlsx(tmp_path):
    rec = {'id': '1', 'meta': {'metadata': {'title': 'T', 'doi': '10.5281/zenodo.1', 'license': {'id': 'cc-by-4.0'},
                                            'creators': [{'name': 'A'}], 'publication_date': '2024-01-01'}}}
    df = pd.DataFrame({'NO3 (umol/L)': ['1', '2'], 'Station': ['A', 'B']})
    r = dt.export_table(rec, 'nut.csv', 'nut.csv', df, {}, 'chemistry', str(tmp_path))
    assert r['csv'] == 'chemistry/1/nut.csv' and r['rows'] == 2 and 'Nutrients' in r['parameters']
    sheets = pd.read_excel(tmp_path / r['xlsx'], sheet_name=None)
    assert set(sheets) == {'Data', 'Columns', 'Source'}
    assert sheets['Columns'].loc[0, 'unit (from header)'] == 'umol/L'


def test_netcdf_json_then_tables(tmp_path):
    import json
    import numpy as np
    import xarray as xr
    ds = xr.Dataset({'temp': (('time', 'depth'), np.array([[10.5, 9.0], [11.0, np.nan]], dtype='float32'),
                              {'units': 'degC', 'long_name': 'sea water temperature'})},
                    coords={'time': pd.to_datetime(['2020-01-01', '2020-01-02']), 'depth': [0.0, 10.0]},
                    attrs={'title': 'test station'})
    nc = tmp_path / 'station.nc'
    ds.to_netcdf(nc)
    rec = {'id': '9', 'meta': {'metadata': {'title': 'T', 'doi': '', 'license': {'id': 'cc-by-4.0'}}}}
    r = dt.export_netcdf(rec, 'station.nc', str(nc), 'physics', str(tmp_path / 'out'))
    doc = json.loads((tmp_path / 'out' / r['json']).read_text(encoding='utf-8'))
    assert doc['global_attributes']['title'] == 'test station'
    assert doc['variables']['temp']['data'] == [[10.5, 9.0], [11.0, None]]
    assert doc['variables']['time']['data'][0].startswith('2020-01-01')
    table = pd.read_csv(tmp_path / 'out' / r['csv'])
    assert list(table.columns) == ['time', 'depth', 'temp'] and len(table) == 4
    assert set(pd.read_excel(tmp_path / 'out' / r['xlsx'], sheet_name=None)) >= {'Data', 'Variables', 'Source'}
    grid, skipped = dt.json_tables(doc, max_rows=3)
    assert grid == [] and skipped
