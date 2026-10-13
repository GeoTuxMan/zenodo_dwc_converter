import os
import sys
import zipfile

import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import dwc_converter as dc  # noqa: E402

DATA = os.path.join(HERE, 'data')


class FakeWorms:
    """Offline stand-in for the WoRMS matcher."""
    IDS = {'Skeletonema costatum': 149074, 'Chaetoceros': 148985, 'Calanus helgolandicus': 104466,
           'Acartia clausi': 149755, 'Noctiluca scintillans': 109391, 'Mytilus galloprovincialis': 140481,
           'Nephtys hombergii': 130359, 'Abra alba': 141433, 'Pleurobrachia pileus': 106386, 'Prorocentrum cordatum': 109634}

    def match(self, names):
        return {n: ({'AphiaID': self.IDS[n], 'rank': 'Species', 'kingdom': 'Animalia', 'phylum': '',
                     'class': '', 'order': '', 'family': '', 'genus': n.split()[0], 'authority': ''}
                    if n in self.IDS else None) for n in names}


def test_normalise_and_map():
    cols = ['Station', 'Date', 'Latitude (N)', 'Longitude (E)', 'Depth (m)', 'Species', 'Abundance (cells/L)']
    m = dc.map_columns(cols)
    assert m['Latitude (N)'] == 'decimalLatitude'
    assert m['Longitude (E)'] == 'decimalLongitude'
    assert m['Species'] == 'scientificName'
    assert m['Date'] == 'eventDate'
    assert m['Station'] == 'locationID'
    assert m['Depth (m)'] == '_depth'
    assert 'Abundance (cells/L)' not in m          # measurement column -> eMoF
    assert dc.map_columns(['decimalLatitude', 'AphiaID'])['AphiaID'] == 'scientificNameID'


def test_override():
    m = dc.map_columns(['col A', 'Species'], {'col A': 'locality'})
    assert m['col A'] == 'locality'


@pytest.mark.parametrize('raw,exp', [
    ('24.03.2010', '2010-03-24'), ('2010-08-24', '2010-08-24'), ('3/24/2010', '2010-03-24'),
    ('2010', '2010'), ('2006-2015', '2006/2015'), ('40261', '2010-03-24'),
    ('2010-03-24 00:00:00', '2010-03-24'), ('12 Mar 2011', '2011-03-12'), ('garbage', None)])
def test_dates(raw, exp):
    assert dc.parse_date(raw) == exp


def test_coordinates():
    assert dc.parse_coordinate('44,25', True) == 44.25
    assert abs(dc.parse_coordinate("44°08.3'N", True) - 44.138333) < 1e-5
    assert dc.parse_coordinate("28 38 54 W", False) < 0
    assert dc.parse_coordinate('S 12.5', True) == -12.5
    assert dc.parse_coordinate('95', True) is None


def test_depth_and_taxa():
    assert dc.parse_depth('0-10') == (0, 10)
    assert dc.parse_depth('-12') == (12, 12)
    assert dc.clean_taxon('Chaetoceros sp.') == ('Chaetoceros', None, 'genus')
    assert dc.clean_taxon('cf. Acartia clausi')[0:2] == ('Acartia clausi', 'cf.')
    assert dc.clean_taxon('Total')[0] is None
    assert dc.clean_taxon('Actinoptychus senarius (Ehr.) Ehr.')[0] == 'Actinoptychus senarius'
    assert dc.clean_taxon('Aulacoseira granulata var. angustissima (O. Muller) Simonsen')[0] == \
        'Aulacoseira granulata var. angustissima'
    assert dc.clean_taxon('123')[0] is None


def test_header_detection_and_long_layout():
    tables = dc.read_tables(os.path.join(DATA, 'long_sample.csv'), DATA)
    assert len(tables) == 1
    name, df = tables[0]
    assert list(df.columns)[0] == 'Station'           # title row skipped
    m = dc.map_columns(list(df.columns))
    built = dc.build_from_table(df, m, 'T1', name, 'Phytoplankton', '')
    assert built['layout'] == 'long'
    assert len(built['occ']) == 3                     # 'Total' row dropped
    assert len(built['events']) == 2
    statuses = {o['scientificName']: o['occurrenceStatus'] for o in built['occ']}
    assert statuses['Chaetoceros'] == 'absent'
    ev = next(e for e in built['events'].values() if e['locationID'] == 'Mamaia')
    assert ev['eventDate'] == '2010-03-24' and ev['minimumDepthInMeters'] == 0 and ev['maximumDepthInMeters'] == 10
    types = {e['measurementType'] for e in built['emof']}
    assert any('per unit volume' in t for t in types)
    abund = [e for e in built['emof'] if e['measurementUnit'] == 'Number per litre']
    assert abund and abund[0]['measurementTypeID'].endswith('SDBIOL01/')


def test_wide_species_columns():
    name, df = dc.read_tables(os.path.join(DATA, 'wide_sample.csv'), DATA)[0]
    m = dc.map_columns(list(df.columns))
    built = dc.build_from_table(df, m, 'T2', name, 'Zooplankton', '')
    assert built['layout'] == 'wide'
    assert len(built['events']) == 2
    # 8 cells, 2 empty -> 6 occurrences, 1 zero -> absent
    assert len(built['occ']) == 6
    assert sum(o['occurrenceStatus'] == 'absent' for o in built['occ']) == 1
    assert {e['eventDate'] for e in built['events'].values()} == {'2015-06', '2015-07'}


def test_wide_needs_known_taxa(monkeypatch):
    df = pd.DataFrame({'Site': ['A', 'B'], 'Year': ['2010', '2011'], 'Foraging': ['1', '2'],
                       'Grouping': ['0', '3'], 'Shelter': ['2', '2']})
    m = dc.map_columns(list(df.columns))
    assert dc.detect_layout(df, m)[0] == 'wide_species_cols'          # no validator: heuristic only
    monkeypatch.setattr(dc, 'TAXON_VALIDATOR', lambda names: set())
    assert dc.detect_layout(df, m)[0] is None                         # WoRMS knows none of them


def test_measurement_from_unit_only():
    m = dc.map_columns(['Meiobenthos - taxaname', 'datecollected', 'Observed Individual Count (ind/m2)'])
    assert m['Meiobenthos - taxaname'] == 'scientificName' and m['datecollected'] == 'eventDate'
    assert 'Observed Individual Count (ind/m2)' not in m
    assert dc.measurement_kind('Observed Individual Count (ind/m2)') == 'abundance'
    assert dc.measurement_kind('ObservedWeight (mg/m2)') == 'biomass'


def test_wide_sample_columns():
    name, df = dc.read_tables(os.path.join(DATA, 'transposed_sample.csv'), DATA)[0]
    m = dc.map_columns(list(df.columns))
    built = dc.build_from_table(df, m, 'T3', name, 'Benthos', '')
    assert built['layout'] == 'wide'
    assert len(built['events']) == 4
    assert len(built['occ']) == 12


def test_missing_date_taken_from_same_sample_or_station():
    # MA12: undated row of the 0-13 m layer; MA14: the 0-14 m layer has 2 dates (19 vs 1 rows) and
    # its undated row takes the dominant one; MA20: no layer match, single station date; MA02: none
    f = ['MA12', 'MA12', 'MA12'] + ['MA14'] * 21 + ['MA20', 'MA20', 'MA02']
    y = [2024, None, 2024] + [2024] * 20 + [None, 2023, None, None]
    m = [7, None, 7] + [7] * 20 + [None, 6, None, None]
    d = [24, None, 24] + [15] * 19 + [24, None, 2, None, None]
    lat = [43.76] * 3 + [43.80] * 21 + [43.70, 43.70, 43.75]
    lon = [30.38] * 3 + [30.79] * 21 + [30.10, 30.10, 30.04]
    dmin = [0, 0, 13] + [0] * 21 + [0, 5, 0]
    dmax = [13, 13, 47] + [14] * 21 + [10, 20, 15]
    df = pd.DataFrame({'FieldNumber': f, 'YearCollected': y, 'MonthCollected': m, 'DayCollected': d,
                       'Latitude': lat, 'Longitude': lon, 'MinimumDepth': dmin, 'MaximumDepth': dmax,
                       'ScientificName': ['Acartia clausi', 'Oithona davisae'] * 13 + ['Acartia clausi'],
                       'Abundance (ind/m3)': list(range(1, 28))})
    built = dc.build_from_table(df, dc.map_columns(list(df.columns)), 'T9', 't', 'Zooplankton', '')
    ev = {(e['fieldNumber'], e['minimumDepthInMeters']): e['eventDate'] for e in built['events'].values()}
    assert ev[('MA12', 0)] == '2024-07-24' and ev[('MA14', 0)] in ('2024-07-15', '2024-07-24')
    assert ev[('MA20', 0)] == '2023-06-02' and ev[('MA02', 0)] == ''
    assert '' not in [e['eventDate'] for e in built['events'].values() if e['fieldNumber'] != 'MA02']
    assert built['stats']['filled_dates'] == 3


def test_institution_code_from_data():
    df = pd.DataFrame({'institutionCode': ['CNR-ISMAR-VE'] * 2, 'eventDate': ['2001-05-02'] * 2,
                       'decimalLatitude': [45.3, 45.3], 'decimalLongitude': [12.5, 12.5],
                       'scientificName': ['Acartia tonsa', 'Oithona nana']})
    built = dc.build_from_table(df, dc.map_columns(list(df.columns)), 'T8', 't', 'Zooplankton', '')
    assert [e['institutionCode'] for e in built['events'].values()] == ['CNR-ISMAR-VE']


def test_institution_override(tmp_path):
    f = tmp_path / 'institutions.json'
    f.write_text('{"7": {"institutionCode": "227", "institutionID": "https://edmo.seadatanet.org/report/227"}}')
    df = pd.DataFrame({'id': ['e1'], 'institutionCode': ['CNR-ISMAR-VE']})
    out = dc.add_dataset_terms(df, {'metadata': {'title': 'x'}}, '7', str(f))
    assert out.loc[0, 'institutionCode'] == '227' and out.loc[0, 'institutionID'].endswith('/227')
    assert dc.add_dataset_terms(df, {'metadata': {'title': 'x'}}, '8', str(f)).loc[0, 'institutionCode'] == 'CNR-ISMAR-VE'


def test_archive_roundtrip(tmp_path):
    name, df = dc.read_tables(os.path.join(DATA, 'long_sample.csv'), DATA)[0]
    built = dc.build_from_table(df, dc.map_columns(list(df.columns)), 'T1', name, 'Phytoplankton', '')
    res = {'n_taxa': 0, 'worms_matched_pct': None, 'issues': []}
    occ = dc.apply_worms(pd.DataFrame(built['occ']), FakeWorms(), res)
    assert res['worms_matched_pct'] == 100.0
    assert occ['scientificNameID'].str.startswith('urn:lsid:marinespecies.org:taxname:').all()
    ev = dc.tidy_frame(list(built['events'].values()), dc.EVENT_TERMS, 'eventID')
    oc = dc.tidy_frame(occ.to_dict('records'), dc.OCC_TERMS, 'occurrenceID')
    oc['id'] = oc['eventID']
    em = dc.tidy_frame(built['emof'], dc.EMOF_TERMS, 'eventID').drop(columns=['eventID'])
    zp = str(tmp_path / 'x_dwca.zip')
    meta = {'metadata': {'title': 'Test and data', 'creators': [{'name': 'Doe, Jane', 'affiliation': 'NIMRD'}],
                         'description': '<p>desc</p>', 'license': {'id': 'cc-by-4.0'},
                         'publication_date': '2026-01-01'}}
    dc.write_dwca(zp, ('event.txt', dc.ROW_EVENT, ev),
                  [('occurrence.txt', dc.ROW_OCCURRENCE, oc), ('extendedmeasurementorfact.txt', dc.ROW_EMOF, em)],
                  dc.eml_xml(meta, 'x', ('2010-03-24', '2010-08-24'), (28.6, 28.7, 44.1, 44.3)))
    assert dc.verify_dwca(zp) == []
    with zipfile.ZipFile(zp) as z:
        assert {'meta.xml', 'eml.xml', 'event.txt', 'occurrence.txt',
                'extendedmeasurementorfact.txt'} <= set(z.namelist())
    # an existing DwC-A is read back through meta.xml
    d = tmp_path / 'unz'
    dc.safe_extract(zp, str(d))
    tables, uris = dc.read_dwca_folder(str(d))
    assert tables[0][1] == dc.ROW_EVENT and len(tables[0][3]) == 2


class FakeRangeSession:
    """Serves a local file through HTTP Range semantics."""
    def __init__(self, path):
        self.data = open(path, 'rb').read()
        self.requests = 0

    def get(self, url, params=None, timeout=None, headers=None, stream=False):
        self.requests += 1
        a, b = headers['Range'].split('=')[1].split('-')
        body = self.data[int(a):int(b) + 1]

        class R:
            status_code = 206
            content = body

        return R()


def test_remote_zip_extracts_only_tables(tmp_path):
    zpath = tmp_path / 'big.zip'
    with zipfile.ZipFile(zpath, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('images/a.png', os.urandom(3_000_000))
        z.writestr('data/occ.csv', 'species,lat,lon\nAbra alba,44.1,29.2\n')
        z.writestr('README.txt', 'notes')
    sess = FakeRangeSession(zpath)
    n, k, mb, notes = dc.remote_zip_extract(sess, 'http://x', os.path.getsize(zpath), str(tmp_path / 'out'))
    assert (n, k) == (3, 1)
    assert (tmp_path / 'out' / 'data' / 'occ.csv').read_text().startswith('species')
    assert mb < 1.5  # the 3 MB image is never transferred


def test_reference_list_is_not_occurrence_data():
    df = pd.DataFrame({'Scientific_name': [f'Genus{i} species' for i in range(600)],
                       'Family': ['Fam'] * 600, 'Order': ['Ord'] * 600})
    m = dc.map_columns(list(df.columns))
    assert dc.is_reference_list('Fish_3letter_Codes_2014.xls::Sheet1', df, m)
    obs = df.assign(Latitude=44.0, Longitude=29.0)
    assert not dc.is_reference_list('Fish_3letter_Codes_2014.xls::Sheet1', obs, dc.map_columns(list(obs.columns)))


def test_taxon_by_sample_matrix_is_not_reference_list():
    df = pd.DataFrame({'Species': [f'Genus{i} species' for i in range(600)],
                       'Family': ['Fam'] * 600, 'Order': ['Ord'] * 600})
    for k in range(6):
        df[f'S{k}'] = range(600)
    assert not dc.is_reference_list('ASV_taxonomy.csv', df, dc.map_columns(list(df.columns)))


def test_dataset_terms_from_zenodo_metadata():
    meta = {'metadata': {'title': 'Black Sea  benthos', 'doi': '10.5281/zenodo.1',
                         'creators': [{'name': 'A', 'affiliation': None}, {'name': 'B', 'affiliation': 'VLIZ'}],
                         'license': {'id': 'cc-by-nc-4.0'}}}
    t = dc.dataset_terms(meta, '1')
    assert t == {'datasetName': 'Black Sea benthos', 'datasetID': 'https://doi.org/10.5281/zenodo.1',
                 'institutionCode': 'VLIZ', 'institutionID': '',
                 'license': 'http://creativecommons.org/licenses/by-nc/4.0/legalcode'}
    df = dc.add_dataset_terms(pd.DataFrame({'id': ['e1', 'e2'], 'institutionCode': ['VLIZ', '']}), meta, '1')
    assert list(df['institutionCode']) == ['VLIZ', 'VLIZ']
    assert 'institutionID' not in df.columns


def test_dataset_terms_edmo_code():
    meta = {'metadata': {'title': 'x', 'creators': [
        {'name': 'A', 'affiliation': 'University of Bucharest'},
        {'name': 'Institutul Naţional de Cercetare Dezvoltare Marină Grigore Antipa', 'affiliation': None}]}}
    t = dc.dataset_terms(meta, '5')
    assert t['institutionCode'] == '697'
    assert t['institutionID'] == 'https://edmo.seadatanet.org/report/697'
    meta['metadata']['creators'] = [{'name': 'B', 'affiliation': 'National Institute for Marine Research '
                                     'and Development "Grigore Antipa" (NIMRD): Constanta, RO'}]
    assert dc.dataset_terms(meta, '5')['institutionCode'] == '697'
    for name in ('GeoEcoMar', 'National Institute of Marine Geology and Geoecology'):
        meta['metadata']['creators'] = [{'name': name, 'affiliation': None}]
        assert dc.dataset_terms(meta, '5')['institutionCode'] == '850'
