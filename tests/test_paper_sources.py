import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import dwc_converter as dc  # noqa: E402
import paper_sources as ps  # noqa: E402


def test_coordinates_dms_and_decimal_minutes():
    pairs = ps.coord_pairs('off Elba, Italy (42°43.617 N, 10°09.598 E) and at Hornøya (70°22′N, 31°09′E, Fig. 2)')
    assert [(round(a, 3), round(b, 3)) for a, b, *_ in pairs] == [(42.727, 10.16), (70.367, 31.15)]
    area = ps.study_area([(p[0], p[1], p[4]) for p in pairs[1:]])
    assert 1000 < area['uncertainty_m'] < 2000            # whole minutes: ~1.3 km
    south_west = ps.coord_pairs('station at 33°52′S 151°12′W')[0]
    assert south_west[0] < 0 and south_west[1] < 0


def test_sampling_period_ignores_citations_units_and_page_heads():
    text = ('Sampling was conducted from July 2013 to August 2014 off the island. '
            'Samples were analysed following Smith et al. (2010) and Jones & Li (2012). '
            'Organisms between 55 and 2000 μm were collected. '
            'Frontiers in Marine Science | www.frontiersin.org 2 April 2018 | Volume 5 | samples')
    period, ev = ps.sampling_period(text, 2018)
    assert period == '2013-07/2014-08' and len(ev) == 1


def test_linked_dois_skip_data_repositories():
    meta = {'metadata': {'doi': '10.5281/zenodo.1', 'description': 'Data for doi:10.1002/ece3.5570.',
                         'related_identifiers': [
                             {'identifier': '10.5061/dryad.x1', 'relation': 'isDerivedFrom', 'scheme': 'doi'},
                             {'identifier': 'https://doi.org/10.3354/meps13101', 'relation': 'isCitedBy'},
                             {'identifier': '10.15468/dl.abc', 'relation': 'references'}]}}
    assert ps.linked_dois(meta) == ['10.3354/meps13101', '10.1002/ece3.5570']


def test_region_check():
    assert ps.in_regions(54.5, 12.0, 'Baltic Sea')
    assert not ps.in_regions(78.0, 113.0, 'Atlantic')
    assert ps.in_regions(10.0, 100.0, '')


def test_paper_values_fill_only_missing():
    events = {'a': {'eventID': 'a', 'locality': 'St. 05'}, 'b': {'eventID': 'b', 'eventDate': '2015-06-01'},
              'c': {'eventID': 'c'}}
    paper = {'paper_doi': '10.1/x', 'eventDate': '2014/2016', 'date_apply': True,
             'stations': {dc.norm_station('Station 5'): ['Station 5', 54.1, 12.2]},
             'area': {'lat': 54.0, 'lon': 12.0, 'uncertainty_m': 5000, 'n_points': 3}, 'area_apply': True}
    notes = dc.derive_event_gaps(events, {'metadata': {}}, None, paper)
    assert events['a']['eventDate'] == '2014/2016' and events['b']['eventDate'] == '2015-06-01'
    assert events['a']['decimalLatitude'] == 54.1 and 'station table' in events['a']['georeferenceRemarks']
    assert events['c']['decimalLatitude'] == 54.0 and events['c']['coordinateUncertaintyInMeters'] == 5000
    assert 'doi:10.1/x' in events['c']['georeferenceRemarks'] and len(notes) == 3
    rejected = {'x': {'eventID': 'x'}}
    dc.derive_event_gaps(rejected, {'metadata': {}}, None, dict(paper, date_apply=False, area_apply=False))
    assert 'eventDate' not in rejected['x'] and 'decimalLatitude' not in rejected['x']
