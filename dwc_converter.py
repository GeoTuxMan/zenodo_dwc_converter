"""
Zenodo -> Darwin Core Archive converter (Zenodo Inventory project, Task 3).

Takes the marine biology datasets discovered on Zenodo (Task 1, data.json), downloads their
tabular files, maps the columns to Darwin Core terms and writes one Darwin Core Archive per
dataset (Event core + Occurrence + ExtendedMeasurementOrFact, or Occurrence core when the data
have no sampling events), ready for the EMODnet Biology / EurOBIS IPT.

Pipeline per dataset:
  1. Zenodo API metadata + file list; download tabular files to a shared cache.
  2. Read every table (csv/tsv/txt, xlsx/xls/ods, zip; DwC-A zips are validated + normalised).
  3. Detect the header row, map columns through a synonym dictionary (or a user override).
  4. Detect the layout: long (one row per taxon observation) or wide (species x samples matrix).
  5. Build Event / Occurrence / eMoF tables, normalise dates, coordinates and depths.
  6. Match taxon names against WoRMS -> scientificNameID (LSID), rank and higher taxonomy.
  7. Write meta.xml + eml.xml (from the Zenodo metadata) and zip the archive.
  8. Record status, mappings and issues in conversion_results.json / .csv / .xlsx.
"""

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import os
import re
import shutil
import sys
import time
import unicodedata
import zipfile
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape as xml_escape

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

__version__ = '1.1.0'
TOOL = f'zenodo_dwc_converter v{__version__}'

ZENODO_API = 'https://zenodo.org/api/records'
WORMS_MATCH = 'https://www.marinespecies.org/rest/AphiaRecordsByMatchNames'

TABULAR_EXT = {'.csv', '.tsv', '.txt', '.tab', '.xlsx', '.xls', '.xlsm', '.ods', '.zip'}
MAX_FILE_MB = 30
MAX_FILES_PER_DATASET = 15
MAX_ROWS_PER_TABLE = 300_000
MAX_OCCURRENCES = 500_000
MAX_UNZIPPED_MB = 200
# --large-files: files above MAX_FILE_MB are read in parts, not downloaded whole
LARGE_TEXT_HEAD_MB = 50        # big text tables: only the first 50 MB (HTTP Range)
LARGE_EXCEL_MAX_MB = 150       # big spreadsheets are downloaded whole up to this size
REMOTE_MEMBER_MAX_MB = 500     # members of a remote zip larger than this are skipped

STATUS_CONVERTED = 'Converted'
STATUS_PARTIAL = 'Partial'
STATUS_FAILED = 'Failed'
STATUS_NOT_ATTEMPTED = 'Not attempted'

# Fill-rate thresholds used to decide whether a required term is "present".
REQUIRED_FILL = 0.9
SCINAMEID_FILL = 0.8

NS_DWC = 'http://rs.tdwg.org/dwc/terms/'
NS_DC = 'http://purl.org/dc/terms/'
NS_OBIS = 'http://rs.iobis.org/obis/terms/'
ROW_EVENT = NS_DWC + 'Event'
ROW_OCCURRENCE = NS_DWC + 'Occurrence'
ROW_EMOF = NS_OBIS + 'ExtendedMeasurementOrFact'

DC_TERMS = {'type', 'modified', 'language', 'license', 'rightsHolder', 'accessRights',
            'bibliographicCitation', 'references'}
OBIS_TERMS = {'measurementTypeID', 'measurementValueID', 'measurementUnitID'}


def term_uri(term):
    if term in DC_TERMS:
        return NS_DC + term
    if term in OBIS_TERMS:
        return NS_OBIS + term
    return NS_DWC + term


# ---------------------------------------------------------------------------
# Column mapping
# ---------------------------------------------------------------------------
# DwC term -> header synonyms. Headers are normalised (lower case, accents stripped, units in
# brackets removed, punctuation -> space) before comparison, so "Latitude (°N)", "LATITUDE" and
# "latitude_dd" all become "latitude" / "latitude dd".
SYNONYMS = {
    'eventID': ['eventid', 'event id', 'event', 'sample id', 'sampleid', 'sample code', 'sample',
                'sample no', 'sample number', 'sample name', 'sample label', 'haul id', 'haul no', 'haul number',
                'samplingevent', 'sampling event', 'event code', 'id sample', 'cod proba',
                'proba', 'probe', 'muestra'],
    'parentEventID': ['parenteventid', 'parent event id', 'parent event'],
    'eventDate': ['eventdate', 'event date', 'date', 'sampling date', 'sample date',
                  'collection date', 'date collected', 'datecollected', 'date sampled', 'datetime', 'date time',
                  'observation date', 'survey date', 'date of sampling', 'date of collection',
                  'sampledate', 'collectiondate', 'data', 'data recoltarii', 'datum', 'fecha',
                  'timestamp', 'date utc', 'start date', 'date start', 'dato'],
    'eventTime': ['eventtime', 'event time', 'time', 'sampling time', 'hour', 'time utc'],
    'fieldNumber': ['fieldnumber', 'field number'],
    'year': ['year', 'year collected', 'start year collected', 'year of collection', 'yr', 'yyyy', 'an', 'anul', 'anno', 'ano', 'jahr', 'annee', 'sampling year',
             'year sampled', 'survey year', 'year of sampling'],
    'month': ['month', 'month collected', 'start month collected', 'mon', 'mm', 'luna', 'mes', 'monat', 'mois', 'mese', 'sampling month'],
    'day': ['day', 'day collected', 'start day collected', 'dd', 'zi', 'ziua', 'dia', 'tag', 'jour', 'giorno', 'day of month'],
    'decimalLatitude': ['decimallatitude', 'decimal latitude', 'latitude', 'lat', 'lat dd',
                        'lat dec', 'latitude dd', 'latitude decimal', 'latitude n', 'lat n',
                        'latitudine', 'breite', 'latitud', 'y coord', 'y coordinate', 'lat wgs84',
                        'latitude wgs84', 'latdd', 'lat deg', 'latitude deg', 'y'],
    'decimalLongitude': ['decimallongitude', 'decimal longitude', 'longitude', 'lon', 'long',
                         'lng', 'lon dd', 'long dd', 'lon dec', 'longitude dd',
                         'longitude decimal', 'longitude e', 'lon e', 'long e', 'longitudine',
                         'laenge', 'lange', 'longitud', 'x coord', 'x coordinate', 'lon wgs84',
                         'longitude wgs84', 'londd', 'lon deg', 'longitude deg', 'x'],
    # start-of-transect / haul positions: used when decimalLatitude/Longitude are empty
    '_latStart': ['start latitude', 'latitude start', 'lat start', 'start lat', 'lat1',
                  'latitude 1', 'shoot latitude'],
    '_lonStart': ['start longitude', 'longitude start', 'lon start', 'start lon', 'long start',
                  'lon1', 'longitude 1', 'shoot longitude'],
    'minimumDepthInMeters': ['minimumdepthinmeters', 'minimum depth', 'min depth', 'depth min',
                             'mindepth', 'depth from', 'upper depth', 'depth top', 'top depth',
                             'depth upper', 'start depth', 'depth start', 'starting depth'],
    'maximumDepthInMeters': ['maximumdepthinmeters', 'maximum depth', 'max depth', 'depth max',
                             'maxdepth', 'depth to', 'lower depth', 'depth bottom',
                             'depth lower', 'end depth', 'depth end', 'ending depth'],
    # A single depth column fills both minimum and maximum depth (handled in build).
    '_depth': ['depth', 'sampling depth', 'sample depth', 'depth m', 'water depth sampled',
               'adancime', 'adincime', 'profundidad', 'profondeur', 'tiefe', 'profondita',
               'depth of sampling', 'depthinmeters', 'depth meters'],
    'locationID': ['locationid', 'location id', 'station', 'station id', 'stationid',
                   'station code', 'station name', 'station no', 'station number', 'st', 'stn',
                   'stat', 'statie', 'statia', 'estacion', 'site id', 'siteid', 'site code',
                   'transect', 'plot'],
    'locality': ['locality', 'location', 'site', 'site name', 'sitename', 'place', 'area',
                 'locatie', 'localitate', 'localidad', 'lieu', 'ort', 'location name',
                 'sampling site', 'sampling location', 'verbatimlocality', 'verbatim locality'],
    'waterBody': ['waterbody', 'water body', 'sea', 'basin', 'ecoregion', 'marine region'],
    'bibliographicCitation': ['bibliographiccitation', 'source', 'reference', 'references',
                              'citation', 'literature', 'bibliography'],
    'country': ['country', 'tara', 'pais', 'pays', 'land'],
    'countryCode': ['countrycode', 'country code', 'iso country'],
    'habitat': ['habitat', 'biotope', 'substrate', 'substratum', 'sediment type', 'habitat type'],
    'coordinateUncertaintyInMeters': ['coordinateuncertaintyinmeters',
                                      'coordinate uncertainty', 'coordinate precision'],
    'samplingProtocol': ['samplingprotocol', 'sampling protocol', 'sampling method', 'method',
                         'gear', 'sampling gear', 'gear type', 'instrument', 'protocol',
                         'sampler', 'metoda', 'metodo'],
    'sampleSizeValue': ['samplesizevalue', 'sample size', 'sample volume', 'volume filtered',
                        'filtered volume', 'sampled volume', 'volume sampled', 'sample area',
                        'area sampled', 'sampled area', 'surface sampled'],
    'scientificName': ['scientificname', 'scientific name', 'scientific name accepted',
                       'scientific name valid', 'valid scientific name', 'accepted scientific name',
                       'species', 'species name',
                       'speciesname', 'taxon', 'taxa', 'taxon name', 'taxonname', 'taxa name', 'taxaname',
                       'specie', 'specia', 'specii', 'organism', 'organism name', 'valid name',
                       'validname', 'accepted name', 'acceptedname', 'latin name', 'latinname',
                       'genus species', 'aphia name', 'taxonomic name', 'scientific',
                       'especie', 'espece', 'art', 'specie name', 'taxon valid', 'taxa valid',
                       'name', 'scientificname accepted', 'verbatimidentification',
                       'original name', 'originalname', 'identification'],
    'scientificNameID': ['scientificnameid', 'scientific name id', 'aphiaid', 'aphia id',
                         'aphia', 'worms id', 'wormsid', 'worms aphiaid', 'lsid',
                         'valid aphiaid', 'aphiaid accepted', 'taxonid'],
    'occurrenceID': ['occurrenceid', 'occurrence id', 'record id', 'recordid',
                     'observation id', 'observationid', 'occurrence'],
    'catalogNumber': ['catalognumber', 'catalog number', 'catalogue number', 'specimen id',
                      'voucher'],
    'individualCount': ['individualcount', 'individual count', 'count', 'counts', 'number',
                        'n', 'no', 'nr', 'numar', 'number of individuals', 'no individuals',
                        'no of individuals', 'nb individuals', 'individuals', 'ind',
                        'total count', 'total number', 'quantity', 'num', 'n individuals',
                        'abundance count', 'count n', 'number individuals', 'numbers',
                        'number of specimens', 'specimens'],
    'occurrenceStatus': ['occurrencestatus', 'occurrence status', 'presence',
                         'presence absence', 'presenceabsence', 'pa', 'p a', 'm p a', 'present', 'detected',
                         'detection'],
    'basisOfRecord': ['basisofrecord', 'basis of record'],
    'recordedBy': ['recordedby', 'recorded by', 'observer', 'observers', 'collector',
                   'collectors', 'collected by', 'sampled by', 'recorder'],
    'identifiedBy': ['identifiedby', 'identified by', 'identifier', 'determiner',
                     'determined by', 'taxonomist'],
    'lifeStage': ['lifestage', 'life stage', 'stage', 'development stage', 'developmental stage',
                  'age class'],
    'sex': ['sex', 'gender'],
    'kingdom': ['kingdom', 'regnum', 'regn'],
    'phylum': ['phylum', 'division', 'increngatura'],
    'class': ['class', 'clasa', 'classis', 'clase', 'klasse'],
    'order': ['order', 'ordin', 'ordo', 'orden', 'ordnung'],
    'family': ['family', 'familia', 'familie', 'famille'],
    'genus': ['genus', 'gen', 'genul', 'genero', 'gattung', 'genre'],
    'scientificNameAuthorship': ['scientificnameauthorship', 'scientific name author', 'authority', 'authorship',
                                 'taxon authority', 'species authority'],
    'taxonRank': ['taxonrank', 'taxon rank', 'rank'],
    'vernacularName': ['vernacularname', 'vernacular name', 'common name', 'commonname',
                       'english name'],
    'identificationQualifier': ['identificationqualifier', 'identification qualifier',
                                'qualifier'],
    'occurrenceRemarks': ['occurrenceremarks', 'occurrence remarks', 'remarks', 'remark',
                          'notes', 'note', 'comment', 'comments', 'observatii', 'observations'],
    'eventRemarks': ['eventremarks', 'event remarks', 'station remarks', 'sample remarks'],
    'organismQuantity': ['organismquantity', 'organism quantity'],
    'organismQuantityType': ['organismquantitytype', 'organism quantity type'],
}

EVENT_TERMS = ['eventID', 'parentEventID', 'fieldNumber', 'eventDate', 'eventTime', 'year', 'month', 'day',
               'locationID', 'locality', 'waterBody', 'country', 'countryCode',
               'decimalLatitude', 'decimalLongitude', 'geodeticDatum',
               'coordinateUncertaintyInMeters', 'minimumDepthInMeters', 'maximumDepthInMeters',
               'habitat', 'samplingProtocol', 'sampleSizeValue', 'sampleSizeUnit',
               'eventRemarks', 'datasetName', 'datasetID', 'institutionCode', 'license']
OCC_TERMS = ['occurrenceID', 'eventID', 'basisOfRecord', 'occurrenceStatus', 'catalogNumber',
             'scientificName', 'scientificNameID', 'scientificNameAuthorship', 'taxonRank',
             'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'vernacularName',
             'verbatimIdentification', 'identificationQualifier', 'individualCount',
             'organismQuantity', 'organismQuantityType', 'lifeStage', 'sex', 'recordedBy',
             'identifiedBy', 'occurrenceRemarks', 'bibliographicCitation']
EMOF_TERMS = ['eventID', 'occurrenceID', 'measurementType', 'measurementTypeID',
              'measurementValue', 'measurementUnit', 'measurementUnitID', 'measurementMethod',
              'measurementRemarks']

REQUIRED_EVENT = ['eventID', 'eventDate', 'decimalLatitude', 'decimalLongitude']
REQUIRED_OCC = ['occurrenceID', 'scientificName', 'scientificNameID', 'occurrenceStatus',
                'basisOfRecord']

# Measurement columns -> eMoF. NERC vocabulary IDs are only given where the concept is
# unambiguous; otherwise measurementTypeID is left empty for the data manager to assign.
P01 = 'http://vocab.nerc.ac.uk/collection/P01/current/'
P06 = 'http://vocab.nerc.ac.uk/collection/P06/current/'
UNIT_IDS = {
    'number per litre': P06 + 'UCPL/', 'number per cubic metre': P06 + 'UPMM/',
    'number per square metre': P06 + 'UPMS/', 'milligrams per cubic metre': P06 + 'UMMC/',
    'percent': P06 + 'UPCT/', 'metres': P06 + 'ULAA/', 'dimensionless': P06 + 'UUUU/',
}
MEASUREMENT_KEYWORDS = [
    # (regex on normalised header, kind)
    (r'\bbiomass\b|\bweight\b|\bbiomasa\b|\bgreutate\b|\bww\b|\bdw\b|\bafdw\b', 'biomass'),
    (r'\babundance\b|\bdensity\b|\bdensitate\b|\babundenta\b|\bdensidad\b|\bcells? per\b'
     r'|\bind(ividuals)? per\b|\bconcentration\b|\bcpue\b', 'abundance'),
    (r'\bcover\b|\bcoverage\b|\bacoperire\b', 'cover'),
    (r'\bfrequency\b|\bfrecventa\b', 'frequency'),
    (r'\blength\b|\bsize\b|\blungime\b', 'size'),
]

BASIS_OF_RECORD = ['HumanObservation', 'MachineObservation', 'MaterialSample', 'PreservedSpecimen',
                   'LivingSpecimen', 'FossilSpecimen', 'MaterialCitation', 'Occurrence']
# Columns that must hold text; a numeric column mapped to them is a false match.
TEXT_TERMS = {'scientificName', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
              'locality', 'country', 'recordedBy', 'identifiedBy', 'basisOfRecord', 'waterBody'}


def norm_basis(v):
    """Return a valid DwC basisOfRecord value or '' (then the default is used)."""
    key = re.sub(r'[^a-z]', '', str(v or '').lower())
    for b in BASIS_OF_RECORD:
        if key == b.lower():
            return b
    aliases = {'o': 'HumanObservation', 'observation': 'HumanObservation',
               'humanobs': 'HumanObservation', 's': 'PreservedSpecimen',
               'specimen': 'PreservedSpecimen', 'sample': 'MaterialSample',
               'literature': 'MaterialCitation', 'machine': 'MachineObservation'}
    return aliases.get(key, '')
ABIOTIC_KEYWORDS = [
    (r'\btemperature\b|\btemp\b|\btemperatura\b|\bsst\b', 'Temperature of the water body', 'degrees Celsius'),
    (r'\bsalinity\b|\bsalinitate\b|\bsal\b|\bpsu\b', 'Salinity of the water body', 'Dimensionless'),
    (r'\boxygen\b|\boxigen\b|\bdo\b|\bo2\b', 'Concentration of oxygen in the water body', None),
    (r'\bchlorophyll\b|\bchl a\b|\bchla\b|\bclorofila\b', 'Concentration of chlorophyll-a in the water body', None),
    (r'\bbottom depth\b|\bseafloor depth\b|\bwater depth\b|\bsounding\b', 'Seafloor depth', 'Metres'),
]

# Values in a taxon column that are clearly not taxa (totals, blanks, codes).
NON_TAXON_RE = re.compile(
    r'^(total|totals|sum|suma|grand total|nan|none|null|n/?a|na|-+|\?+|unknown|other|others|'
    r'varia|unidentified|no data|nd|empty|blank|mean|average|sd|std|min|max)$', re.I)
LATIN_COL_RE = re.compile(
    r'^[A-Z][a-z]{2,}(?:[ _.][a-z]{3,}| (?:sp|spp|cf|aff)\.?(?: [a-z]{3,})?)?(?:[ _.][a-z]{3,})?$')


def normalise_header(h):
    """Lower case, strip accents, drop text in brackets (units), punctuation -> spaces."""
    s = str(h)
    s = unicodedata.normalize('NFKD', s).encode('ascii', 'ignore').decode('ascii')
    s = re.sub(r'([a-z])([A-Z])', r'\1 \2', s)          # camelCase -> camel Case
    s = s.lower()
    s = re.sub(r'\(.*?\)|\[.*?\]|\{.*?\}', ' ', s)
    s = re.sub(r'[^a-z0-9]+', ' ', s).strip()
    return s


def header_unit(h):
    m = re.search(r'[\(\[]([^\)\]]+)[\)\]]', str(h))
    return m.group(1).strip() if m else ''


_SYN_INDEX = None


def _syn_index():
    global _SYN_INDEX
    if _SYN_INDEX is None:
        _SYN_INDEX = {}
        for term, syns in SYNONYMS.items():
            for s in [term.lower()] + syns:
                _SYN_INDEX.setdefault(normalise_header(s), term)
    return _SYN_INDEX


def map_columns(columns, override=None):
    """
    Map source columns to DwC terms. Returns {source column: term}; each term is used once.
    Exact normalised match first; then "camelCase" DwC names; then a prefix match
    (e.g. "latitude decimal degrees" -> decimalLatitude).
    """
    idx = _syn_index()
    mapping, used = {}, set()
    override = override or {}
    for col in columns:
        if col in override and override[col]:
            mapping[col] = override[col]
            used.add(override[col])
    free = [c for c in columns if c not in mapping and c not in override]

    def assign(cands):
        # lower rank (position in the synonym list) wins, e.g. "ScientificName_accepted"
        # beats "Species" (which in some schemas holds only the epithet)
        for _rank, _pos, col, t in sorted(cands):
            if col not in mapping and t not in used:
                mapping[col] = t
                used.add(t)

    # Pass 1: exact
    cands = []
    for pos, col in enumerate(free):
        n = normalise_header(col)
        t = idx.get(n) or idx.get(n.replace(' ', ''))
        if t:
            cands.append((_syn_rank(t, n), pos, col, t))
    assign(cands)
    # Pass 2: leading words (e.g. "latitude decimal degrees", "station name ro"), except
    # bookkeeping columns such as "date last modified" / "year identified"
    cands = []
    for pos, col in enumerate(free):
        if col in mapping:
            continue
        n = normalise_header(col)
        if re.search(r'\b(modified|identified|created|updated|entered|entry|end|stop|qc|flag|'
                     r'error|accuracy|precision|method|source|unit|units|duration|distance|speed|time)\b', n):
            continue
        if measurement_kind(col):
            continue
        words = n.split()
        for k in range(len(words) - 1, 0, -1):
            prefix = ' '.join(words[:k])
            t = idx.get(prefix)
            if t and len(prefix) >= 3:
                cands.append((100 + _syn_rank(t, prefix), pos, col, t))
                break
    assign(cands)
    # Pass 3: trailing words for the key terms ("Meiobenthos - taxaname", "Sampling station
    # code" ...)
    cands = []
    for pos, col in enumerate(free):
        if col in mapping:
            continue
        if measurement_kind(col):
            continue
        words = normalise_header(col).split()
        for k in range(1, len(words)):
            suffix = ' '.join(words[k:])
            t = idx.get(suffix)
            if t in ('scientificName', 'eventDate', 'eventID', 'locationID', 'decimalLatitude',
                     'decimalLongitude', 'individualCount') and len(suffix) >= 4:
                cands.append((200 + _syn_rank(t, suffix), pos, col, t))
                break
    assign(cands)
    return mapping


def _syn_rank(term, normalised):
    syns = [normalise_header(s) for s in [term] + SYNONYMS.get(term, [])]
    return syns.index(normalised) if normalised in syns else len(syns)


def measurement_kind(col):
    n = normalise_header(col)
    u = header_unit(col).lower().replace(' ', '').replace('µ', 'u')
    if re.search(r'(^|[^a-z])[umk]?g(/|per)|g\*?m-?[23]', u):          # mass per area / volume
        return 'biomass'
    for rx, kind in MEASUREMENT_KEYWORDS:
        if re.search(rx, n):
            return kind
    if re.search(r'(ind|cells?|n|no|#|individuals?)(/|per)(l|ml|m2|m3|m\^?[23])', u):
        return 'abundance'
    return None


def abiotic_kind(col):
    n = normalise_header(col)
    for rx, mtype, unit in ABIOTIC_KEYWORDS:
        if re.search(rx, n):
            return mtype, unit
    return None


def describe_measurement(kind, unit_text, context=''):
    """Return (measurementType, measurementTypeID, measurementUnit, measurementUnitID)."""
    u = (unit_text or '').lower().replace(' ', '')
    u = u.replace('^', '').replace('³', '3').replace('²', '2').replace('µ', 'u')
    per_l = bool(re.search(r'(/|per)(l|litre|liter|dm3)$|l-1$|/l\b', u))
    per_ml = bool(re.search(r'(/|per)(ml|cm3)$|ml-1$', u))
    per_m3 = bool(re.search(r'(/|per)m3$|m-3$', u))
    per_m2 = bool(re.search(r'(/|per)m2$|m-2$', u))
    ctx = context.lower()
    if kind == 'abundance' or (kind is None and (per_l or per_m3 or per_m2 or per_ml)):
        if per_m2 or 'benth' in ctx:
            if per_m2:
                return ('Abundance of biological entity specified elsewhere per unit area of the bed',
                        P01 + 'SDBIOL02/', 'Number per square metre', UNIT_IDS['number per square metre'])
        if per_l:
            return ('Abundance of biological entity specified elsewhere per unit volume of the water body',
                    P01 + 'SDBIOL01/', 'Number per litre', UNIT_IDS['number per litre'])
        if per_m3:
            return ('Abundance of biological entity specified elsewhere per unit volume of the water body',
                    P01 + 'SDBIOL01/', 'Number per cubic metre', UNIT_IDS['number per cubic metre'])
        if per_ml:
            return ('Abundance of biological entity specified elsewhere per unit volume of the water body',
                    P01 + 'SDBIOL01/', 'Number per millilitre', '')
        return ('Abundance', '', unit_text or '', '')
    if kind == 'biomass':
        words = normalise_header(ctx) + ' ' + u
        mtype = 'Wet weight biomass' if re.search(r'\b(wet|ww|fresh)\b', words) else (
            'Dry weight biomass' if re.search(r'\b(dry|dw|afdw)\b', words) else 'Biomass')
        if re.search(r'mg/m3|mgm-3|mg\*m-3', u):
            return (mtype + ' per unit volume of the water body', '', 'Milligrams per cubic metre',
                    UNIT_IDS['milligrams per cubic metre'])
        return (mtype, '', unit_text or '', '')
    if kind == 'cover':
        return ('Percentage cover', '', 'Percent', UNIT_IDS['percent'])
    if kind == 'count':
        return ('Count of biological entity specified elsewhere', P01 + 'OCOUNT01/',
                'Dimensionless', UNIT_IDS['dimensionless'])
    if kind == 'frequency':
        return ('Frequency of occurrence', '', unit_text or '', '')
    if kind == 'size':
        return ('Size', '', unit_text or '', '')
    return (kind or 'Measurement', '', unit_text or '', '')


# ---------------------------------------------------------------------------
# Value parsing
# ---------------------------------------------------------------------------
def to_number(v):
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return None if pd.isna(v) else float(v)
    s = str(v).strip().replace(' ', '').replace(' ', '')
    if not s or s.lower() in ('nan', 'na', 'n/a', '-', 'nd', 'none', 'null'):
        return None
    if re.fullmatch(r'-?\d{1,3}(\.\d{3})+,\d+', s):          # 1.234,5
        s = s.replace('.', '').replace(',', '.')
    elif re.fullmatch(r'-?\d+,\d+', s):                      # 12,5
        s = s.replace(',', '.')
    elif re.fullmatch(r'-?\d{1,3}(,\d{3})+(\.\d+)?', s):     # 1,234.5
        s = s.replace(',', '')
    s = re.sub(r'^[<>~≈]=?', '', s)
    try:
        return float(s)
    except ValueError:
        return None


DMS_RE = re.compile(
    r'^\s*(-)?\s*(\d{1,3})(?:[°º:\s]+|\s*deg\s*)(\d{1,2}(?:[.,]\d+)?)?\s*[\'′:\s]*'
    r'(\d{1,2}(?:[.,]\d+)?)?\s*(?:"|″|\'\')?\s*([NSEWnsew])?\s*$')


def parse_coordinate(v, is_lat):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    s = str(v).strip()
    if not s:
        return None
    hemi = None
    m = re.match(r'^\s*([NSEWnsew])\s*(.*)$', s)
    if m:
        hemi, s = m.group(1).upper(), m.group(2)
    num = to_number(s)
    if num is None:
        m = DMS_RE.match(s)
        if not m:
            return None
        neg, d, mi, se, h = m.groups()
        num = float(d) + (float(mi.replace(',', '.')) / 60 if mi else 0) + \
            (float(se.replace(',', '.')) / 3600 if se else 0)
        if neg:
            num = -num
        hemi = (h or hemi or '').upper() or None
    else:
        m = re.search(r'([NSEWnsew])\s*$', str(v))
        if m:
            hemi = m.group(1).upper()
    if hemi in ('S', 'W') and num > 0:
        num = -num
    limit = 90 if is_lat else 180
    if num is None or abs(num) > limit:
        return None
    return round(num, 6)


MONTHS = {m: i for i, m in enumerate(
    ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'], 1)}


def parse_date(v):
    """Return an ISO 8601 date / date range string, or None."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    if isinstance(v, (dt.date, dt.datetime, pd.Timestamp)):
        return pd.Timestamp(v).strftime('%Y-%m-%d')
    s = str(v).strip()
    if not s or s.lower() in ('nan', 'nat', 'none', 'na', 'n/a', '-'):
        return None
    m = re.fullmatch(r'(\d{4})-(\d{2})-(\d{2})(?:[ T](\d{2}):(\d{2})(?::(\d{2}))?(?:\.\d+)?)?(Z|[+-]\d{2}:?\d{2})?', s)
    if m:
        y, mo, d, hh, mi = m.group(1), m.group(2), m.group(3), m.group(4), m.group(5)
        if hh and (hh, mi) != ('00', '00'):
            return f'{y}-{mo}-{d}T{hh}:{mi}' + (m.group(7) or '')
        return f'{y}-{mo}-{d}'
    if re.fullmatch(r'\d{4}(-\d{2}(-\d{2})?)?/\d{4}(-\d{2}(-\d{2})?)?', s):
        return s
    m = re.fullmatch(r'(1[6-9]\d{2}|20\d{2})\s*[-–/]\s*(1[6-9]\d{2}|20\d{2})', s)
    if m:
        return f'{m.group(1)}/{m.group(2)}'
    if re.fullmatch(r'(1[6-9]\d{2}|20\d{2})', s):
        return s
    if re.fullmatch(r'(1[6-9]\d{2}|20\d{2})\.0', s):
        return s[:4]
    m = re.fullmatch(r'(1[6-9]\d{2}|20\d{2})[-/.](\d{1,2})', s)
    if m and 1 <= int(m.group(2)) <= 12:
        return f'{m.group(1)}-{int(m.group(2)):02d}'
    m = re.fullmatch(r'(\d{1,2})[./-](\d{1,2})[./-](\d{2,4})(?:\s+(\d{1,2}):(\d{2}).*)?', s)
    if m:
        a, b, y = int(m.group(1)), int(m.group(2)), int(m.group(3))
        if y < 100:
            y += 2000 if y < 50 else 1900
        d, mo = (a, b) if a > 12 or b <= 12 else (b, a)   # day-first unless impossible
        if b > 12 and a <= 12:
            d, mo = b, a
        try:
            return dt.date(y, mo, d).isoformat()
        except ValueError:
            return None
    m = re.fullmatch(r'(\d{4})[./](\d{1,2})[./](\d{1,2})', s)
    if m:
        try:
            return dt.date(int(m.group(1)), int(m.group(2)), int(m.group(3))).isoformat()
        except ValueError:
            return None
    m = re.fullmatch(r'(\d{1,2})[\s-]([A-Za-z]{3})[a-z]*[\s-](\d{4})', s)
    if m and m.group(2).lower() in MONTHS:
        try:
            return dt.date(int(m.group(3)), MONTHS[m.group(2).lower()], int(m.group(1))).isoformat()
        except ValueError:
            return None
    n = to_number(s)
    if n is not None and 20000 < n < 60000:                  # Excel serial date
        return (dt.date(1899, 12, 30) + dt.timedelta(days=int(n))).isoformat()
    try:
        ts = pd.to_datetime(s, dayfirst=True, errors='raise')
        if 1600 < ts.year <= dt.date.today().year:
            return ts.strftime('%Y-%m-%d')
    except Exception:
        pass
    return None


def parse_depth(v):
    """Return (min, max) depth in metres from '5', '0-10', '-12' ..."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None, None
    s = str(v).strip()
    m = re.fullmatch(r'(\d+(?:[.,]\d+)?)\s*[-–/]\s*(\d+(?:[.,]\d+)?)\s*m?', s)
    if m:
        a, b = to_number(m.group(1)), to_number(m.group(2))
        return min(a, b), max(a, b)
    n = to_number(re.sub(r'\s*m$', '', s))
    if n is None:
        return None, None
    n = abs(n)
    return n, n


def clean_taxon(name):
    """Return (scientificName, identificationQualifier, rank hint) from a raw taxon label."""
    if name is None:
        return None, None, None
    s = str(name).replace('_', ' ').replace(' ', ' ')
    s = re.sub(r'\s+', ' ', s).strip().strip('"\'')
    if not s or NON_TAXON_RE.match(s) or to_number(s) is not None or len(s) < 3:
        return None, None, None
    if not re.search(r'[A-Za-z]{3}', s):
        return None, None, None
    qualifier, rank = None, None
    m = re.search(r'\b(cf|aff)\.?\s', s + ' ')
    if m:
        qualifier = m.group(1) + '.'
        s = re.sub(r'\b(cf|aff)\.?\s*', '', s).strip()
    if re.search(r'\s(sp|spp|sp\.\s*\d*|spp\.|indet|n\.?\s*sp)\.?\s*\d*$', s, re.I):
        s = re.sub(r'\s(sp|spp|indet|n\.?\s*sp)\.?\s*\d*.*$', '', s, flags=re.I).strip()
        rank = 'genus'
    s = strip_authorship(s)
    s = re.sub(r'\s*\(.*?\)\s*$', '', s).strip() or s
    return s, qualifier, rank


def strip_authorship(s):
    """'Actinoptychus senarius (Ehr.) Ehr.' -> 'Actinoptychus senarius';
    'Aulacoseira granulata var. angustissima (O. Müller) Simonsen' keeps the variety."""
    toks = s.split(' ')
    if len(toks) < 3 or not re.match(r'^[A-Z][a-z]+$', toks[0]):
        return s
    out = [toks[0]]
    i = 1
    if i < len(toks) and re.match(r'^\([A-Z][a-z]+\)$', toks[i]):        # subgenus
        out.append(toks[i])
        i += 1
    expect_epithet = True
    while i < len(toks):
        t = toks[i]
        if expect_epithet and re.match(r'^[a-z][a-z-]+$', t):
            out.append(t)
            expect_epithet = False
        elif t in ('var.', 'subsp.', 'ssp.', 'f.', 'forma', 'var', 'subsp'):
            out.append(t)
            expect_epithet = True
        else:
            break
        i += 1
    return ' '.join(out)


def stable_id(*parts, n=12):
    return hashlib.sha1('|'.join(str(p) for p in parts).encode('utf-8')).hexdigest()[:n]


# ---------------------------------------------------------------------------
# Reading files
# ---------------------------------------------------------------------------
def detect_header_row(raw, max_scan=25):
    """Pick the header row of a header-less DataFrame (title rows above are skipped)."""
    best, best_score = 0, -1
    n = min(max_scan, len(raw))
    width = raw.shape[1]
    for i in range(n):
        row = raw.iloc[i]
        vals = [str(x).strip() for x in row if x is not None and str(x).strip() not in ('', 'nan', 'None')]
        if len(vals) < max(2, 0.4 * width):
            continue
        texty = sum(1 for x in vals if to_number(x) is None and len(x) < 80)
        score = texty / len(vals) + len(vals) / width
        if texty >= max(2, 0.6 * len(vals)) and score > best_score:
            best, best_score = i, score
            if len(vals) >= 0.8 * width:
                break
    return best


def frame_from_raw(raw):
    if raw is None or raw.empty:
        return None
    raw = raw.dropna(how='all').dropna(axis=1, how='all')
    if raw.shape[0] < 2 or raw.shape[1] < 2:
        return None
    raw = raw.reset_index(drop=True)
    h = detect_header_row(raw)
    header, seen = [], {}
    for j, x in enumerate(raw.iloc[h]):
        name = str(x).strip() if x is not None and str(x).strip() not in ('', 'nan', 'None') else f'col_{j + 1}'
        if name in seen:
            seen[name] += 1
            name = f'{name}_{seen[name]}'
        else:
            seen[name] = 0
        header.append(name)
    df = raw.iloc[h + 1:].copy()
    df.columns = header
    df = df.dropna(how='all')
    df = df.loc[:, [c for c in df.columns if not (c.startswith('col_') and df[c].isna().all())]]
    if df.shape[0] < 1 or df.shape[1] < 2:
        return None
    return df.head(MAX_ROWS_PER_TABLE).reset_index(drop=True)


def sniff_text(path):
    with open(path, 'rb') as fh:
        sample = fh.read(200_000)
    for enc in ('utf-8-sig', 'utf-8', 'cp1252', 'latin-1'):
        try:
            text = sample.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    else:
        return None, None
    if '\x00' in text[:2000]:
        return None, None
    lines = [l for l in text.splitlines()[:60] if l.strip()]
    if len(lines) < 2:
        return enc, None
    best, best_score = None, 0
    for d in ['\t', ';', ',', '|']:
        counts = [l.count(d) for l in lines[-40:]]
        if not counts or max(counts) == 0:
            continue
        mode = max(set(counts), key=counts.count)
        consistency = counts.count(mode) / len(counts)
        if mode >= 1 and consistency >= 0.6 and mode * consistency > best_score:
            best, best_score = d, mode * consistency
    if best is None:
        # whitespace-delimited tables
        counts = [len(l.split()) for l in lines[-40:]]
        mode = max(set(counts), key=counts.count)
        if mode >= 3 and counts.count(mode) / len(counts) >= 0.7:
            best = r'\s+'
    return enc, best


def read_text_table(path):
    """Read a delimited text file. Rows are padded so that title lines above the header
    (with fewer fields) do not truncate the table, which pandas.read_csv would do."""
    enc, sep = sniff_text(path)
    if not sep:
        return []
    try:
        with open(path, encoding=enc, errors='replace', newline='') as fh:
            if sep == r'\s+':
                rows = [line.split() for _, line in zip(range(MAX_ROWS_PER_TABLE + 30), fh)]
            else:
                reader = csv.reader(fh, delimiter=sep)
                rows = [r for _, r in zip(range(MAX_ROWS_PER_TABLE + 30), reader)]
    except Exception:
        return []
    rows = [r for r in rows if any(str(x).strip() for x in r)]
    if len(rows) < 2:
        return []
    width = max(len(r) for r in rows)
    raw = pd.DataFrame([[(x.strip() or None) for x in r] + [None] * (width - len(r)) for r in rows])
    df = frame_from_raw(raw)
    return [(os.path.basename(path), df)] if df is not None else []


def read_excel_tables(path):
    ext = os.path.splitext(path)[1].lower()
    engine = {'.xls': 'xlrd', '.ods': 'odf'}.get(ext, 'openpyxl')
    try:
        sheets = pd.read_excel(path, sheet_name=None, header=None, dtype=str, engine=engine)
    except Exception:
        return []
    out = []
    for name, raw in sheets.items():
        if re.search(r'readme|metadata|legend|info|description|notes|codes|abbrev', str(name), re.I) \
                and raw.shape[0] < 200:
            continue
        df = frame_from_raw(raw)
        if df is not None and df.shape[0] >= 1:
            out.append((f'{os.path.basename(path)}::{name}', df))
    return out


def safe_extract(zpath, dest):
    """Extract a zip, refusing absolute / parent paths and oversized content."""
    with zipfile.ZipFile(zpath) as z:
        total = sum(i.file_size for i in z.infolist())
        if total > MAX_UNZIPPED_MB * 1024 * 1024:
            raise ValueError(f'zip too large when extracted ({total / 1e6:.0f} MB)')
        for info in z.infolist():
            name = info.filename
            if name.startswith('/') or '..' in name.replace('\\', '/').split('/') or ':' in name:
                continue
            if '__MACOSX' in name:
                continue
            target = os.path.join(dest, name)
            if info.is_dir():
                os.makedirs(target, exist_ok=True)
                continue
            os.makedirs(os.path.dirname(target), exist_ok=True)
            with z.open(info) as src, open(target, 'wb') as dst:
                shutil.copyfileobj(src, dst)


def find_dwca_root(folder):
    for root, _dirs, files in os.walk(folder):
        if 'meta.xml' in files:
            return root
    return None


def read_tables(path, extract_root, depth=0):
    """Return [(table name, DataFrame)] for any supported file (zip handled recursively)."""
    ext = os.path.splitext(path)[1].lower()
    base = os.path.basename(path).lower()
    if re.match(r'(readme|license|licence|citation|changelog|manifest|checksum|md5)', base):
        return []
    if ext in ('.csv', '.tsv', '.txt', '.tab', '.dat'):
        return read_text_table(path)
    if ext in ('.xlsx', '.xlsm', '.xls', '.ods'):
        return read_excel_tables(path)
    if ext == '.zip' and depth < 2:
        dest = os.path.join(extract_root, re.sub(r'[^\w.-]', '_', os.path.basename(path)))
        if not os.path.isdir(dest):
            safe_extract(path, dest)
        out = []
        for root, _dirs, files in os.walk(dest):
            for f in sorted(files):
                if os.path.splitext(f)[1].lower() in TABULAR_EXT | {'.dat'}:
                    for name, df in read_tables(os.path.join(root, f), dest, depth + 1):
                        out.append((f'{os.path.basename(path)}/{name}', df))
                if len(out) > 30:
                    return out
        return out
    return []


# ---------------------------------------------------------------------------
# WoRMS
# ---------------------------------------------------------------------------
class Worms:
    def __init__(self, cache_file, session):
        self.cache_file = cache_file
        self.session = session
        self.cache = {}
        if os.path.exists(cache_file):
            try:
                with open(cache_file, encoding='utf-8') as fh:
                    self.cache = json.load(fh)
            except Exception:
                self.cache = {}

    def save(self):
        tmp = f'{self.cache_file}.{os.getpid()}.part'
        with open(tmp, 'w', encoding='utf-8') as fh:
            json.dump(self.cache, fh)
        os.replace(tmp, self.cache_file)

    def match(self, names):
        """Return {name: record dict or None}."""
        todo = [n for n in dict.fromkeys(names) if n and n not in self.cache]
        for i in range(0, len(todo), 50):
            chunk = todo[i:i + 50]
            params = [('scientificnames[]', n) for n in chunk] + [('marine_only', 'false')]
            try:
                r = self.session.get(WORMS_MATCH, params=params, timeout=90)
                if r.status_code == 204:
                    res = [[] for _ in chunk]
                elif r.status_code == 200:
                    res = r.json()
                else:
                    continue
            except Exception:
                continue
            for name, hits in zip(chunk, res):
                rec = None
                for h in hits or []:
                    if h.get('match_type') in ('exact', 'exact_genus', 'exact_subgenus', 'phonetic',
                                               'near_1', 'near_2') and h.get('AphiaID'):
                        rec = {k: h.get(k) for k in ('AphiaID', 'scientificname', 'authority', 'rank',
                                                     'kingdom', 'phylum', 'class', 'order', 'family',
                                                     'genus', 'status', 'valid_AphiaID', 'match_type',
                                                     'isMarine', 'isBrackish')}
                        break
                self.cache[name] = rec
            time.sleep(0.3)
        if todo:
            self.save()
        return {n: self.cache.get(n) for n in names}


# ---------------------------------------------------------------------------
# Building DwC tables from one source table
# ---------------------------------------------------------------------------
def looks_numeric(series, min_frac=0.7):
    vals = series.dropna().astype(str).str.strip()
    vals = vals[vals != '']
    if len(vals) == 0:
        return False
    return sum(to_number(v) is not None for v in vals.head(500)) / min(len(vals), 500) >= min_frac


# Optional callable(names) -> set of names recognised as taxa (set to WoRMS during a run).
TAXON_VALIDATOR = None


def detect_layout(df, mapping):
    """Return ('long'|'wide_species_cols'|'wide_sample_cols'|None, value columns)."""
    unmapped = [c for c in df.columns if c not in mapping]
    has_taxon = 'scientificName' in mapping.values()
    if not has_taxon:
        latin = [c for c in unmapped if LATIN_COL_RE.match(str(c).strip())
                 and not measurement_kind(c) and not abiotic_kind(c)]
        latin_num = [c for c in latin if looks_numeric(df[c], 0.5)]
        if TAXON_VALIDATOR and latin_num:
            # headers like "Foraging", "Sand", "Cod" look Latin; keep only names known to WoRMS
            names = {c: clean_taxon(c)[0] for c in latin_num}
            valid = TAXON_VALIDATOR([n for n in names.values() if n])
            ok = [c for c in latin_num if names[c] in valid]
            latin_num = ok if len(ok) >= 0.6 * len(latin_num) else []
        if len(latin_num) >= 3:
            return 'wide_species_cols', latin_num
        return None, []
    taxon_col = [c for c, t in mapping.items() if t == 'scientificName'][0]
    meas = [c for c in unmapped if measurement_kind(c) and looks_numeric(df[c], 0.5)]
    count_cols = [c for c, t in mapping.items() if t in ('individualCount', 'organismQuantity')]
    numeric_unmapped = [c for c in unmapped if c not in meas and not abiotic_kind(c)
                        and looks_numeric(df[c], 0.6)]
    if not meas and not count_cols and len(numeric_unmapped) >= 3 and \
            df[taxon_col].nunique() >= 0.8 * len(df.dropna(subset=[taxon_col])):
        return 'wide_sample_cols', numeric_unmapped
    return 'long', meas


def column_values(df, mapping, term):
    cols = [c for c, t in mapping.items() if t == term]
    return df[cols[0]] if cols else None


def build_from_table(df, mapping, record_id, table_name, context, eml_hint, default_taxon=None):
    """
    Turn one source table into rows for event / occurrence / emof.
    Returns dict with lists of dicts + diagnostics, or None if the table holds no taxon data.

    default_taxon: for single-species datasets (taxon only in the dataset title), tables with
    sampling information but no taxon column get this constant scientificName.
    """
    mapping = {c: t for c, t in mapping.items()
               if not (t in TEXT_TERMS and c in df.columns and looks_numeric(df[c], 0.8))}
    layout, value_cols = detect_layout(df, mapping)
    used_default = False
    if layout is None and default_taxon and \
            set(mapping.values()) & {'eventDate', 'year', 'decimalLatitude'}:
        df = df.copy()
        df['__taxon'] = default_taxon
        mapping = dict(mapping)
        mapping['__taxon'] = 'scientificName'
        layout, value_cols = detect_layout(df, mapping)
        used_default = layout is not None
    if layout is None:
        return None
    ctx = f'{table_name} {context}'.lower()
    dataset_kind = 'biomass' if re.search(r'biomass|weight', table_name.lower()) else 'abundance'
    count_like = bool(re.search(r'count|number|individuals|ind\b', table_name.lower()))

    # --- reshape wide -> long -------------------------------------------------
    if layout == 'wide_species_cols':
        id_cols = [c for c in df.columns if c not in value_cols]
        long = df.melt(id_vars=id_cols, value_vars=value_cols, var_name='__taxon',
                       value_name='__value')
        mapping = dict(mapping)
        mapping['__taxon'] = 'scientificName'
        meas_cols = ['__value']
        df = long
    elif layout == 'wide_sample_cols':
        id_cols = [c for c in df.columns if c not in value_cols]
        long = df.melt(id_vars=id_cols, value_vars=value_cols, var_name='__sample',
                       value_name='__value')
        mapping = dict(mapping)
        if 'eventID' in mapping.values():
            for c, t in list(mapping.items()):
                if t == 'eventID':
                    mapping[c] = 'eventRemarks'
        mapping['__sample'] = 'eventID'
        meas_cols = ['__value']
        df = long
    else:
        meas_cols = value_cols

    get = lambda term: column_values(df, mapping, term)

    taxon_raw = get('scientificName')
    n = len(df)
    rows_event, rows_occ, rows_emof = {}, [], []

    date_s, year_s, month_s, day_s = get('eventDate'), get('year'), get('month'), get('day')
    lat_s, lon_s = get('decimalLatitude'), get('decimalLongitude')
    lat0_s, lon0_s = get('_latStart'), get('_lonStart')
    dmin_s, dmax_s, depth_s = get('minimumDepthInMeters'), get('maximumDepthInMeters'), get('_depth')
    ev_s = get('eventID')
    count_s = get('individualCount')
    status_s = get('occurrenceStatus')
    occid_s = get('occurrenceID')

    simple_event = ['parentEventID', 'fieldNumber', 'eventTime', 'locationID', 'locality', 'waterBody', 'country',
                    'countryCode', 'habitat', 'coordinateUncertaintyInMeters', 'samplingProtocol',
                    'sampleSizeValue', 'eventRemarks']
    simple_occ = ['catalogNumber', 'scientificNameAuthorship', 'taxonRank', 'kingdom', 'phylum',
                  'class', 'order', 'family', 'genus', 'vernacularName', 'identificationQualifier',
                  'organismQuantity', 'organismQuantityType', 'lifeStage', 'sex', 'recordedBy',
                  'identifiedBy', 'occurrenceRemarks', 'scientificNameID', 'basisOfRecord',
                  'bibliographicCitation']
    ev_series = {t: get(t) for t in simple_event}
    occ_series = {t: get(t) for t in simple_occ}
    sample_unit = ''
    for c, t in mapping.items():
        if t == 'sampleSizeValue':
            sample_unit = header_unit(c)

    abiotic_cols = [(c, abiotic_kind(c)) for c in df.columns
                    if c not in mapping and c not in meas_cols and abiotic_kind(c)
                    and looks_numeric(df[c], 0.6)]

    basis = 'MaterialSample' if re.search(r'\bedna\b|metabarcod|environmental dna', eml_hint.lower()) \
        else 'HumanObservation'
    stats = {'bad_dates': 0, 'bad_coords': 0, 'dropped_taxa': 0}
    occ_seen = {}

    def sval(series, i):
        if series is None:
            return ''
        v = series.iat[i]
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return ''
        s = str(v).strip()
        return '' if s.lower() in ('nan', 'none', 'na', 'n/a', 'null', 'nd', '-', '--', '?') else s

    for i in range(n):
        name, qual, rank_hint = clean_taxon(taxon_raw.iat[i])
        if not name:
            stats['dropped_taxa'] += 1
            continue
        # quantity / presence
        values = {}
        for mc in meas_cols:
            num = to_number(df[mc].iat[i])
            if num is not None:
                values[mc] = num
        raw_vals_present = any(sval(df[mc], i) for mc in meas_cols)
        if layout != 'long' and not values:
            continue                      # empty cell in a matrix = not observed
        count = to_number(sval(count_s, i)) if count_s is not None else None
        status = sval(status_s, i).lower()
        if status in ('1', 'yes', 'y', 'true', 'present', 'p', 'x', '+', 'presence', 'da'):
            status = 'present'
        elif status in ('0', 'no', 'n', 'false', 'absent', 'a', '-', 'absence', 'nu'):
            status = 'absent'
        else:
            status = ''
        if not status:
            if values:
                status = 'absent' if all(v == 0 for v in values.values()) else 'present'
            elif count is not None:
                status = 'absent' if count == 0 else 'present'
            else:
                status = 'present'

        # event fields
        date = parse_date(date_s.iat[i]) if date_s is not None else None
        if date is None and year_s is not None:
            y = to_number(sval(year_s, i))
            if y and 1600 < y <= dt.date.today().year:
                date = str(int(y))
                mo = to_number(sval(month_s, i)) if month_s is not None else None
                if mo and 1 <= mo <= 12:
                    date += f'-{int(mo):02d}'
                    d = to_number(sval(day_s, i)) if day_s is not None else None
                    if d and 1 <= d <= 31:
                        try:
                            date = dt.date(int(y), int(mo), int(d)).isoformat()
                        except ValueError:
                            pass
        if date is None and ((date_s is not None and sval(date_s, i)) or (year_s is not None and sval(year_s, i))):
            stats['bad_dates'] += 1
        lat = parse_coordinate(sval(lat_s, i), True) if lat_s is not None else None
        lon = parse_coordinate(sval(lon_s, i), False) if lon_s is not None else None
        if (lat_s is not None and sval(lat_s, i) and lat is None) or (lon_s is not None and sval(lon_s, i) and lon is None):
            stats['bad_coords'] += 1
        if (lat is None or lon is None) and lat0_s is not None and lon0_s is not None:
            lat = parse_coordinate(sval(lat0_s, i), True)
            lon = parse_coordinate(sval(lon0_s, i), False)
        if lat is None or lon is None:
            lat = lon = None
        dmin = dmax = None
        if dmin_s is not None or dmax_s is not None:
            dmin = parse_depth(sval(dmin_s, i))[0] if dmin_s is not None else None
            dmax = parse_depth(sval(dmax_s, i))[1] if dmax_s is not None else None
            if dmin is None:
                dmin = dmax
            if dmax is None:
                dmax = dmin
        elif depth_s is not None:
            dmin, dmax = parse_depth(sval(depth_s, i))
        if dmin is not None and dmax is not None and dmin > dmax:   # e.g. start/end of a tow
            dmin, dmax = dmax, dmin

        ev = {t: sval(ev_series[t], i) for t in simple_event if ev_series[t] is not None}
        ev.update({'eventDate': date or '', 'decimalLatitude': '' if lat is None else lat,
                   'decimalLongitude': '' if lon is None else lon,
                   'minimumDepthInMeters': '' if dmin is None else dmin,
                   'maximumDepthInMeters': '' if dmax is None else dmax})
        if sample_unit and ev.get('sampleSizeValue'):
            ev['sampleSizeUnit'] = sample_unit
        src_ev = sval(ev_s, i) if ev_s is not None else ''
        key_parts = (src_ev,) if src_ev else (ev.get('eventDate'), ev.get('decimalLatitude'),
                                              ev.get('decimalLongitude'), ev.get('locationID', ''),
                                              ev.get('locality', ''), ev.get('minimumDepthInMeters'),
                                              ev.get('maximumDepthInMeters'), ev.get('parentEventID', ''))
        if src_ev:
            event_id = record_id + '_' + re.sub(r'\s+', '_', src_ev)
        else:
            event_id = f'{record_id}_ev_{stable_id(*key_parts)}'
        ev['eventID'] = event_id
        if event_id not in rows_event:
            rows_event[event_id] = ev
            for col, (mtype, unit) in abiotic_cols:
                v = to_number(df[col].iat[i])
                if v is not None:
                    rows_emof.append({'eventID': event_id, 'occurrenceID': '',
                                      'measurementType': mtype, 'measurementTypeID': '',
                                      'measurementValue': v,
                                      'measurementUnit': header_unit(col) or (unit or ''),
                                      'measurementUnitID': UNIT_IDS.get((unit or '').lower(), ''),
                                      'measurementMethod': '', 'measurementRemarks': f'source column: {col}'})
        else:
            # keep the first non-empty value of each field
            for k, v in ev.items():
                if v != '' and rows_event[event_id].get(k, '') == '':
                    rows_event[event_id][k] = v

        occ = {t: sval(occ_series[t], i) for t in simple_occ if occ_series[t] is not None}
        raw_name = str(taxon_raw.iat[i]).strip()
        occ.update({'eventID': event_id, 'scientificName': name,
                    'verbatimIdentification': raw_name if raw_name != name else '',
                    'occurrenceStatus': status,
                    'basisOfRecord': norm_basis(occ.get('basisOfRecord')) or basis})
        if qual and not occ.get('identificationQualifier'):
            occ['identificationQualifier'] = qual
        if rank_hint and not occ.get('taxonRank'):
            occ['taxonRank'] = rank_hint
        if count is not None:
            occ['individualCount'] = int(count) if float(count).is_integer() else ''
            if not float(count).is_integer():
                values['__count_decimal'] = count
        src_occ = sval(occid_s, i) if occid_s is not None else ''
        base_id = f'{record_id}_{src_occ}' if src_occ else \
            f'{record_id}_occ_{stable_id(event_id, name, occ.get("lifeStage", ""), occ.get("sex", ""), raw_name)}'
        k = occ_seen.get(base_id, 0)
        occ_seen[base_id] = k + 1
        occ['occurrenceID'] = base_id if k == 0 else f'{base_id}_{k + 1}'
        rows_occ.append(occ)

        for mc, v in values.items():
            if v == 0 and layout != 'long':
                continue
            if mc in ('__value',):
                kind = 'count' if (count_like and float(v).is_integer()) else dataset_kind
                unit = ''
                mtype = describe_measurement(kind, unit, ctx)
            elif mc == '__count_decimal':
                mtype = describe_measurement('abundance', '', ctx)
            else:
                mtype = describe_measurement(measurement_kind(mc), header_unit(mc), ctx + ' ' + mc)
            rows_emof.append({'eventID': event_id, 'occurrenceID': occ['occurrenceID'],
                              'measurementType': mtype[0], 'measurementTypeID': mtype[1],
                              'measurementValue': v, 'measurementUnit': mtype[2],
                              'measurementUnitID': mtype[3], 'measurementMethod': '',
                              'measurementRemarks': '' if mc.startswith('__') else f'source column: {mc}'})
        if len(rows_occ) >= MAX_OCCURRENCES:
            break

    if not rows_occ:
        return None
    return {'layout': 'long' if layout == 'long' else 'wide', 'events': rows_event, 'occ': rows_occ,
            'emof': rows_emof, 'mapping': mapping, 'stats': stats, 'default_taxon': used_default}


# ---------------------------------------------------------------------------
# DwC-A writing
# ---------------------------------------------------------------------------
def tidy_frame(rows, terms, id_term):
    df = pd.DataFrame(rows)
    cols = [t for t in terms if t in df.columns and (df[t].astype(str).str.strip() != '').any()]
    extra = [c for c in df.columns if c not in terms and not c.startswith('_')
             and (df[c].astype(str).str.strip() != '').any()]
    cols = cols + extra
    if id_term not in cols:
        cols = [id_term] + cols
    df = df[cols].copy()
    df.insert(0, 'id', df[id_term])
    return df


def clean_cell(v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ''
    if isinstance(v, float) and v.is_integer() and abs(v) < 1e15:
        return str(int(v))
    if isinstance(v, str) and re.fullmatch(r'-?\d{1,15}\.0', v):
        return v[:-2]
    return re.sub(r'[\t\r\n]+', ' ', str(v))


def write_tsv(df, path):
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(df.columns) + '\n')
        for row in df.itertuples(index=False):
            fh.write('\t'.join(clean_cell(v) for v in row) + '\n')


def meta_xml(core, extensions, uris=None):
    """core / extensions: (filename, rowType, DataFrame)."""
    uris = uris or {}

    def block(tag, fname, rowtype, df):
        idtag = 'id' if tag == 'core' else 'coreid'
        lines = [f'  <{tag} encoding="UTF-8" fieldsTerminatedBy="\\t" linesTerminatedBy="\\n" '
                 f'fieldsEnclosedBy="" ignoreHeaderLines="1" rowType="{rowtype}">',
                 '    <files>', f'      <location>{fname}</location>', '    </files>',
                 f'    <{idtag} index="0" />']
        for i, c in enumerate(df.columns[1:], start=1):
            lines.append(f'    <field index="{i}" term="{xml_escape(uris.get(c, term_uri(c)))}"/>')
        lines.append(f'  </{tag}>')
        return '\n'.join(lines)

    parts = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<archive xmlns="http://rs.tdwg.org/dwc/text/" metadata="eml.xml">',
             block('core', *core)]
    for ext in extensions:
        parts.append(block('extension', *ext))
    parts.append('</archive>')
    return '\n'.join(parts) + '\n'


LICENSE_URLS = {
    'cc-by-4.0': ('Creative Commons Attribution (CC-BY) 4.0', 'http://creativecommons.org/licenses/by/4.0/legalcode'),
    'cc-zero': ('Public Domain (CC0 1.0)', 'http://creativecommons.org/publicdomain/zero/1.0/legalcode'),
    'cc0-1.0': ('Public Domain (CC0 1.0)', 'http://creativecommons.org/publicdomain/zero/1.0/legalcode'),
    'cc-by-nc-4.0': ('Creative Commons Attribution Non Commercial (CC-BY-NC) 4.0',
                     'http://creativecommons.org/licenses/by-nc/4.0/legalcode'),
    'cc-by-sa-4.0': ('Creative Commons Attribution Share Alike (CC-BY-SA) 4.0',
                     'http://creativecommons.org/licenses/by-sa/4.0/legalcode'),
}


def eml_xml(meta, record_id, temporal=None, bbox=None):
    md = meta.get('metadata', {})
    title = xml_escape(md.get('title', f'Zenodo record {record_id}'))
    doi = md.get('doi') or meta.get('doi') or ''
    desc = re.sub(r'<[^<]+?>', ' ', md.get('description', '') or '')
    desc = xml_escape(re.sub(r'\s+', ' ', desc).strip())
    pub = md.get('publication_date', '')[:10]
    lic_id = ((md.get('license') or {}).get('id') or '').lower()
    lic_name, lic_url = LICENSE_URLS.get(lic_id, (lic_id or 'not specified', ''))

    def party(tag, c):
        name = c.get('name', '')
        if ',' in name:
            sur, given = [x.strip() for x in name.split(',', 1)]
        else:
            bits = name.split()
            given, sur = ' '.join(bits[:-1]), (bits[-1] if bits else name)
        x = [f'    <{tag}>', '      <individualName>']
        if given:
            x.append(f'        <givenName>{xml_escape(given)}</givenName>')
        x.append(f'        <surName>{xml_escape(sur)}</surName>')
        x.append('      </individualName>')
        if c.get('affiliation'):
            x.append(f'      <organizationName>{xml_escape(c["affiliation"])}</organizationName>')
        if c.get('orcid'):
            x.append(f'      <userId directory="https://orcid.org/">{xml_escape(c["orcid"])}</userId>')
        x.append(f'    </{tag}>')
        return '\n'.join(x)

    creators = md.get('creators', []) or [{'name': 'Unknown'}]
    kws = md.get('keywords', []) or []
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<eml:eml xmlns:eml="https://eml.ecoinformatics.org/eml-2.2.0" '
           'xmlns:dc="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
           'xsi:schemaLocation="https://eml.ecoinformatics.org/eml-2.2.0 '
           'https://rs.gbif.org/schema/eml-gbif-profile/1.3/eml.xsd" '
           f'packageId="zenodo-{record_id}-dwca" system="http://gbif.org" scope="system" xml:lang="eng">',
           '  <dataset>']
    if doi:
        out.append(f'    <alternateIdentifier>https://doi.org/{xml_escape(doi)}</alternateIdentifier>')
    out.append(f'    <alternateIdentifier>https://zenodo.org/records/{record_id}</alternateIdentifier>')
    out.append(f'    <title xml:lang="eng">{title}</title>')
    out += [party('creator', c) for c in creators]
    out.append(party('metadataProvider', creators[0]))
    if pub:
        out.append(f'    <pubDate>{pub}</pubDate>')
    out.append('    <language>eng</language>')
    out.append(f'    <abstract><para>{desc or title}</para></abstract>')
    if kws:
        out.append('    <keywordSet>')
        out += [f'      <keyword>{xml_escape(str(k))}</keyword>' for k in kws]
        out.append('      <keywordThesaurus>N/A</keywordThesaurus>')
        out.append('    </keywordSet>')
    if lic_url:
        out.append(f'    <intellectualRights><para>This work is licensed under a <ulink url="{lic_url}">'
                   f'<citetitle>{xml_escape(lic_name)}</citetitle></ulink>.</para></intellectualRights>')
    else:
        out.append(f'    <intellectualRights><para>License: {xml_escape(lic_name)}</para></intellectualRights>')
    out.append(f'    <distribution scope="document"><online><url function="information">'
               f'https://zenodo.org/records/{record_id}</url></online></distribution>')
    if temporal or bbox:
        out.append('    <coverage>')
        if bbox:
            w, e, s, n = bbox
            out += ['      <geographicCoverage>',
                    '        <geographicDescription>Bounding box of the occurrence records</geographicDescription>',
                    '        <boundingCoordinates>',
                    f'          <westBoundingCoordinate>{w}</westBoundingCoordinate>',
                    f'          <eastBoundingCoordinate>{e}</eastBoundingCoordinate>',
                    f'          <northBoundingCoordinate>{n}</northBoundingCoordinate>',
                    f'          <southBoundingCoordinate>{s}</southBoundingCoordinate>',
                    '        </boundingCoordinates>', '      </geographicCoverage>']
        if temporal:
            out += ['      <temporalCoverage><rangeOfDates>',
                    f'        <beginDate><calendarDate>{temporal[0]}</calendarDate></beginDate>',
                    f'        <endDate><calendarDate>{temporal[1]}</calendarDate></endDate>',
                    '      </rangeOfDates></temporalCoverage>']
        out.append('    </coverage>')
    out.append(party('contact', creators[0]))
    out.append('  </dataset>')
    out.append('  <additionalMetadata><metadata><gbif>'
               f'<dateStamp>{dt.datetime.now().isoformat(timespec="seconds")}</dateStamp>'
               f'<citation>{xml_escape("; ".join(c.get("name", "") for c in creators))} ({pub[:4]}). {title}. '
               f'Zenodo. {("https://doi.org/" + xml_escape(doi)) if doi else ""}</citation>'
               f'<resourceLogoUrl></resourceLogoUrl><hierarchyLevel>dataset</hierarchyLevel>'
               '</gbif></metadata></additionalMetadata>')
    out.append('</eml:eml>')
    return '\n'.join(out) + '\n'


def write_dwca(zip_path, core, extensions, eml, uris=None):
    """core/extensions: (filename, rowType, DataFrame)."""
    tmpdir = zip_path + '_tmp'
    os.makedirs(tmpdir, exist_ok=True)
    try:
        for fname, _rt, df in [core] + list(extensions):
            write_tsv(df, os.path.join(tmpdir, fname))
        with open(os.path.join(tmpdir, 'meta.xml'), 'w', encoding='utf-8') as fh:
            fh.write(meta_xml(core, extensions, uris))
        with open(os.path.join(tmpdir, 'eml.xml'), 'w', encoding='utf-8') as fh:
            fh.write(eml)
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
            for f in os.listdir(tmpdir):
                z.write(os.path.join(tmpdir, f), f)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def verify_dwca(zip_path):
    """Check that meta.xml parses, listed files exist, header widths match and core ids are unique."""
    problems = []
    with zipfile.ZipFile(zip_path) as z:
        names = set(z.namelist())
        if 'meta.xml' not in names:
            return ['meta.xml missing']
        root = ET.fromstring(z.read('meta.xml'))
        ns = {'d': 'http://rs.tdwg.org/dwc/text/'}
        for tag in ('core', 'extension'):
            for el in root.findall(f'd:{tag}', ns):
                loc = el.find('d:files/d:location', ns).text
                if loc not in names:
                    problems.append(f'{loc} listed in meta.xml but missing')
                    continue
                with z.open(loc) as fh:
                    header = fh.readline().decode('utf-8').rstrip('\n').split('\t')
                    n_fields = len(el.findall('d:field', ns)) + 1
                    if len(header) != n_fields:
                        problems.append(f'{loc}: {len(header)} columns vs {n_fields} in meta.xml')
                    if tag == 'core':
                        ids = [line.decode('utf-8').split('\t', 1)[0] for line in fh]
                        if len(ids) != len(set(ids)):
                            problems.append(f'{loc}: duplicate core ids')
        if 'eml.xml' in names:
            try:
                ET.fromstring(z.read('eml.xml'))
            except ET.ParseError as e:
                problems.append(f'eml.xml does not parse: {e}')
    return problems


# ---------------------------------------------------------------------------
# Existing DwC-A (validate + normalise)
# ---------------------------------------------------------------------------
def read_dwca_folder(folder):
    root = ET.parse(os.path.join(folder, 'meta.xml')).getroot()
    ns = {'d': 'http://rs.tdwg.org/dwc/text/'}
    tables, uris = [], {}
    for tag in ('core', 'extension'):
        for el in root.findall(f'd:{tag}', ns):
            loc = el.find('d:files/d:location', ns).text
            sep = el.get('fieldsTerminatedBy', ',').encode().decode('unicode_escape')
            quote = el.get('fieldsEnclosedBy', '"').encode().decode('unicode_escape') or None
            skip = int(el.get('ignoreHeaderLines', '0'))
            idx = el.find('d:id', ns) if tag == 'core' else el.find('d:coreid', ns)
            cols = {}
            for f in el.findall('d:field', ns):
                if f.get('index') is not None:
                    term = f.get('term').rstrip('/').split('/')[-1]
                    cols[int(f.get('index'))] = term
                    uris[term] = f.get('term')
            path = os.path.join(folder, loc)
            df = pd.read_csv(path, sep=sep, header=None, skiprows=skip, dtype=str,
                             quotechar=quote or '"', quoting=csv.QUOTE_NONE if not quote else csv.QUOTE_MINIMAL,
                             encoding=el.get('encoding', 'UTF-8'), keep_default_na=False,
                             on_bad_lines='skip', engine='python')
            id_index = int(idx.get('index')) if idx is not None else 0
            named = pd.DataFrame({'id': df[id_index]})
            for i, term in sorted(cols.items()):
                if i < df.shape[1]:
                    named[term] = df[i]
            tables.append((tag, el.get('rowType'), os.path.basename(loc), named))
    return tables, uris


# ---------------------------------------------------------------------------
# Zenodo access + cache
# ---------------------------------------------------------------------------
def robust_session():
    s = requests.Session()
    s.headers['User-Agent'] = f'{TOOL} (EMODnet Biology, NIMRD)'
    retries = Retry(total=3, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504],
                    respect_retry_after_header=True, raise_on_status=False)
    s.mount('https://', HTTPAdapter(max_retries=retries))
    return s


# Circuit breaker: after 3 consecutive server failures Zenodo is considered down for 5 minutes
# and datasets are skipped quickly (marked "retry later") instead of waiting on time-outs.
_BREAKER = {'fails': 0, 'until': 0.0}


def _breaker_check():
    if time.time() < _BREAKER['until']:
        raise RuntimeError('Zenodo unavailable (circuit open)')


def _breaker_fail():
    _BREAKER['fails'] += 1
    if _BREAKER['fails'] >= 3:
        _BREAKER['until'] = time.time() + 300
        _BREAKER['fails'] = 0


def zenodo_record(session, record_id, cache_dir, token=''):
    meta_path = os.path.join(cache_dir, record_id, '_record.json')
    if os.path.exists(meta_path):
        with open(meta_path, encoding='utf-8') as fh:
            return json.load(fh)
    params = {'access_token': token} if token else {}
    _breaker_check()
    try:
        r = session.get(f'{ZENODO_API}/{record_id}', params=params, timeout=(15, 60))
    except requests.RequestException:
        _breaker_fail()
        raise
    time.sleep(1)
    if r.status_code != 200:
        if r.status_code >= 500 or r.status_code == 429:
            _breaker_fail()
        raise RuntimeError(f'Zenodo API {r.status_code}')
    _BREAKER['fails'] = 0
    data = r.json()
    os.makedirs(os.path.dirname(meta_path), exist_ok=True)
    tmp = f'{meta_path}.{os.getpid()}.part'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(data, fh)
    os.replace(tmp, meta_path)
    return data


def file_entries(record):
    files = record.get('files') or []
    if isinstance(files, dict):
        files = files.get('entries', [])
        if isinstance(files, dict):
            files = list(files.values())
    out = []
    for f in files:
        key = f.get('key') or f.get('filename')
        links = f.get('links', {})
        url = links.get('self') or links.get('content') or links.get('download')
        if url and not url.endswith('/content') and '/api/records/' in url and '/files/' in url:
            url = url + '/content'
        out.append({'key': key, 'size': f.get('size') or f.get('filesize') or 0, 'url': url})
    return out


def download(session, url, dest, token=''):
    if os.path.exists(dest):
        return dest
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = f'{dest}.{os.getpid()}.part'
    params = {'access_token': token} if token else {}
    _breaker_check()
    try:
        with session.get(url, params=params, stream=True, timeout=(15, 90)) as r:
            if r.status_code != 200:
                if r.status_code >= 500 or r.status_code == 429:
                    _breaker_fail()
                raise RuntimeError(f'download HTTP {r.status_code}')
            with open(tmp, 'wb') as fh:
                for chunk in r.iter_content(1 << 16):
                    fh.write(chunk)
    except requests.RequestException:
        _breaker_fail()
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    _BREAKER['fails'] = 0
    try:
        os.replace(tmp, dest)
    except OSError:
        if os.path.exists(dest):
            os.remove(tmp)
        else:
            raise
    time.sleep(1)
    return dest


class HttpRangeFile(io.RawIOBase):
    """Read-only, seekable view of a remote file through HTTP Range requests, so that
    zipfile can list and extract single members of a multi-GB archive without
    downloading it. Blocks of 1 MB are cached (16 at most)."""
    BLOCK = 1 << 20

    def __init__(self, session, url, size, token=''):
        super().__init__()
        self.session, self.url, self.size = session, url, size
        self.params = {'access_token': token} if token else {}
        self.pos = 0
        self.blocks = {}
        self.fetched = 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        base = {0: 0, 1: self.pos, 2: self.size}[whence]
        self.pos = max(0, min(self.size, base + offset))
        return self.pos

    def _block(self, n):
        if n not in self.blocks:
            start = n * self.BLOCK
            end = min(self.size, start + self.BLOCK) - 1
            last = None
            for attempt in range(5):
                try:
                    r = self.session.get(self.url, params=self.params, timeout=(15, 90),
                                         headers={'Range': f'bytes={start}-{end}'})
                except requests.RequestException as e:
                    last = e
                else:
                    if r.status_code == 206:
                        break
                    last = RuntimeError(f'range request HTTP {r.status_code}')
                    if r.status_code == 200:
                        raise last  # Range ignored: do not pull the whole file
                time.sleep(10 * (attempt + 1))
            else:
                raise last
            if len(self.blocks) >= 16:
                self.blocks.pop(next(iter(self.blocks)))
            self.blocks[n] = r.content
            self.fetched += len(r.content)
        return self.blocks[n]

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        out = bytearray()
        while n > 0 and self.pos < self.size:
            b = self._block(self.pos // self.BLOCK)
            off = self.pos % self.BLOCK
            chunk = b[off:off + n]
            if not chunk:
                break
            out += chunk
            self.pos += len(chunk)
            n -= len(chunk)
        return bytes(out)

    def readinto(self, buf):
        data = self.read(len(buf))
        buf[:len(data)] = data
        return len(data)


def _trim_to_last_line(path):
    """Drop the incomplete last line of a truncated text file."""
    with open(path, 'rb+') as fh:
        fh.seek(0, 2)
        size = fh.tell()
        fh.seek(max(0, size - (1 << 20)))
        tail = fh.read()
        cut = tail.rfind(b'\n')
        if cut >= 0:
            fh.truncate(size - len(tail) + cut + 1)


def download_head(session, url, dest, nbytes, token=''):
    """Download only the first nbytes of a large text table."""
    if os.path.exists(dest):
        return dest
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = f'{dest}.{os.getpid()}.part'
    params = {'access_token': token} if token else {}
    with session.get(url, params=params, stream=True, timeout=(15, 90),
                     headers={'Range': f'bytes=0-{nbytes - 1}'}) as r:
        if r.status_code not in (200, 206):
            raise RuntimeError(f'download HTTP {r.status_code}')
        got = 0
        with open(tmp, 'wb') as fh:
            for chunk in r.iter_content(1 << 16):
                fh.write(chunk[:nbytes - got])
                got += len(chunk)
                if got >= nbytes:
                    break
    _trim_to_last_line(tmp)
    os.replace(tmp, dest)
    return dest


def remote_zip_extract(session, url, size, dest, token=''):
    """Extract only the tabular members of a large remote zip (HTTP Range).
    Returns (n members, n members extracted, MB transferred, notes)."""
    rf = HttpRangeFile(session, url, size, token)
    notes = []
    with zipfile.ZipFile(io.BufferedReader(rf, buffer_size=1 << 16)) as z:
        infos = [i for i in z.infolist() if not i.is_dir() and '__MACOSX' not in i.filename]
        is_dwca = any(os.path.basename(i.filename).lower() == 'meta.xml' for i in infos)
        want = (TABULAR_EXT - {'.zip'}) | {'.dat'} | ({'.xml'} if is_dwca else set())
        cand = []
        for i in infos:
            name = i.filename
            if name.startswith('/') or '..' in name.replace('\\', '/').split('/') or ':' in name:
                continue
            ext = os.path.splitext(name)[1].lower()
            base = os.path.basename(name).lower()
            if ext not in want or re.match(r'(readme|license|licence|checksum|md5)', base):
                continue
            if i.file_size > REMOTE_MEMBER_MAX_MB * 1024 * 1024:
                notes.append(f'{base}: {i.file_size / 1e6:.0f} MB member skipped')
                continue
            cand.append(i)
        prio = {'.xml': -1, '.csv': 0, '.tsv': 0, '.xlsx': 1, '.xls': 1, '.ods': 1,
                '.txt': 2, '.tab': 2, '.dat': 2}
        cand.sort(key=lambda i: (prio.get(os.path.splitext(i.filename)[1].lower(), 9), i.file_size))
        cand = cand[:MAX_FILES_PER_DATASET * 2]
        head = LARGE_TEXT_HEAD_MB * 1024 * 1024
        for i in cand:
            target = os.path.join(dest, i.filename)
            if len(target) > 230:  # Windows MAX_PATH: flatten deep members
                h = hashlib.md5(i.filename.encode()).hexdigest()[:8]
                target = os.path.join(dest, f'{h}_{os.path.basename(i.filename)[-60:]}')
            os.makedirs(os.path.dirname(target), exist_ok=True)
            text = os.path.splitext(i.filename)[1].lower() in ('.csv', '.tsv', '.txt', '.tab', '.dat')
            cut = text and i.file_size > head
            with z.open(i) as src, open(target, 'wb') as dst:
                if cut:
                    dst.write(src.read(head))
                else:
                    shutil.copyfileobj(src, dst)
            if cut:
                _trim_to_last_line(target)
                notes.append(f'{os.path.basename(i.filename)}: first {LARGE_TEXT_HEAD_MB} MB read')
    return len(infos), len(cand), rf.fetched / 1e6, notes


def fetch_large(session, f, rec_dir, extract_root, token=''):
    """--large-files: get the usable part of a file above MAX_FILE_MB.
    Returns (local path, note). For a zip the path is that of the zip (never downloaded),
    whose extraction folder has been filled; read_tables then walks that folder."""
    key, size = f['key'], f['size'] or 0
    ext = os.path.splitext(key)[1].lower()
    safe = re.sub(r'[\\/:*?"<>|]', '_', key)
    if ext in ('.csv', '.tsv', '.txt', '.tab'):
        dest = os.path.join(rec_dir, '_head', safe)
        download_head(session, f['url'], dest, LARGE_TEXT_HEAD_MB * 1024 * 1024, token)
        return dest, f'{key}: {size / 1e6:.0f} MB, only the first {LARGE_TEXT_HEAD_MB} MB read'
    if ext == '.zip':
        d = os.path.join(extract_root, re.sub(r'[^\w.-]', '_', key))
        marker = os.path.join(d, '_remote_zip.json')
        if os.path.exists(marker):
            with open(marker, encoding='utf-8') as fh:
                info = json.load(fh)
        else:
            n, k, mb, notes = remote_zip_extract(session, f['url'], size, d, token)
            info = {'members': n, 'extracted': k, 'mb': round(mb, 1), 'notes': notes}
            os.makedirs(d, exist_ok=True)
            with open(marker, 'w', encoding='utf-8') as fh:
                json.dump(info, fh)
        note = (f'{key}: {size / 1e6:.0f} MB zip read remotely - {info["extracted"]} of '
                f'{info["members"]} files are tables ({info["mb"]} MB transferred)')
        if info['notes']:
            note += '; ' + '; '.join(info['notes'][:3])
        return os.path.join(rec_dir, safe), note
    if ext in ('.xlsx', '.xlsm', '.xls', '.ods') and size <= LARGE_EXCEL_MAX_MB * 1024 * 1024:
        dest = os.path.join(rec_dir, safe)
        download(session, f['url'], dest, token)
        return dest, f'{key}: {size / 1e6:.0f} MB spreadsheet downloaded whole'
    raise ValueError(f'{size / 1e6:.0f} MB, cannot be read in parts')


# ---------------------------------------------------------------------------
# Dataset conversion
# ---------------------------------------------------------------------------
def record_id_of(item):
    return str(item.get('Metadata URL', '')).rstrip('/').split('/')[-1]


def load_override(mapping_dir, record_id):
    if not mapping_dir:
        return {}
    for ext in ('.json', '.yaml', '.yml'):
        p = os.path.join(mapping_dir, f'{record_id}{ext}')
        if os.path.exists(p):
            with open(p, encoding='utf-8') as fh:
                if ext == '.json':
                    return json.load(fh)
                import yaml
                return yaml.safe_load(fh) or {}
    return {}


def new_result(item, record_id):
    doi = item.get('Data URL/DOI', '')
    return {'record_id': record_id, 'title': item.get('Data title', ''), 'doi': doi,
            'url': f'https://zenodo.org/records/{record_id}', 'taxa': item.get('Taxa covered', ''),
            'regions': item.get('Geographic coverage', ''), 'source_files': [], 'input_format': '',
            'layout': 'unknown', 'core': None, 'status': STATUS_NOT_ATTEMPTED, 'n_events': 0,
            'n_occurrences': 0, 'n_emof': 0, 'n_taxa': 0, 'worms_matched_pct': None,
            'mapped_fields': {}, 'unmapped_columns': [], 'missing_required': [], 'issues': [],
            'dwca_file': None}


def fill_rate(df, col):
    if df is None or col not in df.columns or len(df) == 0:
        return 0.0
    return float((df[col].astype(str).str.strip() != '').mean())


def apply_worms(occ, worms, res):
    names = sorted(set(occ['scientificName']))
    matches = worms.match(names)
    matched = {n: m for n, m in matches.items() if m}
    res['n_taxa'] = len(names)
    res['worms_matched_pct'] = round(100 * len(matched) / len(names), 1) if names else None
    unmatched = [n for n in names if n not in matched]
    if unmatched:
        res['issues'].append(f'{len(unmatched)} taxon names not matched in WoRMS (e.g. '
                             f'{", ".join(unmatched[:5])})')
    existing = occ['scientificNameID'] if 'scientificNameID' in occ.columns else pd.Series([''] * len(occ))

    def lsid(row_name, current):
        cur = str(current or '').strip()
        if re.fullmatch(r'\d+', cur):
            return f'urn:lsid:marinespecies.org:taxname:{cur}'
        if cur.startswith('urn:lsid:marinespecies.org'):
            return cur
        m = matched.get(row_name)
        return f'urn:lsid:marinespecies.org:taxname:{m["AphiaID"]}' if m else ''

    occ['scientificNameID'] = [lsid(n, c) for n, c in zip(occ['scientificName'], existing)]
    for term, key in [('taxonRank', 'rank'), ('kingdom', 'kingdom'), ('phylum', 'phylum'),
                      ('class', 'class'), ('order', 'order'), ('family', 'family'), ('genus', 'genus'),
                      ('scientificNameAuthorship', 'authority')]:
        cur = occ[term] if term in occ.columns else pd.Series([''] * len(occ), index=occ.index)
        occ[term] = [c if str(c).strip() else ((matched.get(n) or {}).get(key) or '')
                     for n, c in zip(occ['scientificName'], cur)]
    if 'taxonRank' in occ.columns:
        occ['taxonRank'] = occ['taxonRank'].astype(str).str.lower().replace('nan', '')
    return occ


def convert_dwca_source(folder, record_id, meta, res, worms, out_dir):
    tables, uris = read_dwca_folder(folder)
    core = [t for t in tables if t[0] == 'core']
    if not core:
        res['status'] = STATUS_FAILED
        res['issues'].append('DwC-A without core in meta.xml')
        return res
    _tag, rowtype, fname, core_df = core[0]
    res['core'] = 'Event' if rowtype.endswith('Event') else ('Occurrence' if rowtype.endswith('Occurrence') else rowtype.split('/')[-1])
    occ_df = core_df if res['core'] == 'Occurrence' else next(
        (t[3] for t in tables if t[1].endswith('Occurrence')), None)
    emof_df = next((t[3] for t in tables if 'MeasurementOrFact' in t[1]), None)
    ev_df = core_df if res['core'] == 'Event' else None
    res['n_events'] = len(ev_df) if ev_df is not None else 0
    res['n_occurrences'] = len(occ_df) if occ_df is not None else 0
    res['n_emof'] = len(emof_df) if emof_df is not None else 0
    res['mapped_fields'] = {f'{t[2]}:{c}': c for t in tables for c in t[3].columns if c != 'id'}
    if core_df['id'].duplicated().any():
        res['issues'].append('duplicate ids in core file')
    if occ_df is not None and 'scientificName' in occ_df.columns:
        occ_df = occ_df.copy()
        occ_df['scientificName'] = occ_df['scientificName'].astype(str).str.strip()
        occ_df = apply_worms(occ_df, worms, res)
    missing = []
    loc = ev_df if ev_df is not None else occ_df
    for t in REQUIRED_EVENT[1:]:
        if fill_rate(loc, t) < REQUIRED_FILL:
            missing.append(t)
    for t in REQUIRED_OCC:
        if occ_df is None or fill_rate(occ_df, t) < (SCINAMEID_FILL if t == 'scientificNameID' else REQUIRED_FILL):
            missing.append(t)
    res['missing_required'] = missing
    # normalised archive: same tables, enriched occurrence taxonomy
    exts = []
    for tag, rt, fn, df in tables:
        if tag == 'core':
            continue
        if rt.endswith('Occurrence') and occ_df is not None and res['core'] != 'Occurrence':
            df = occ_df
        exts.append((fn, rt, df))
    core_out = occ_df if res['core'] == 'Occurrence' and occ_df is not None else core_df
    bbox, temporal = bbox_temporal(loc)
    zip_path = os.path.join(out_dir, 'dwca', f'{record_id}_dwca.zip')
    for c in list(core_out.columns) + [c for _f, _r, d in exts for c in d.columns]:
        uris.setdefault(c, term_uri(c))
    write_dwca(zip_path, (fname, rowtype, core_out), exts, eml_xml(meta, record_id, temporal, bbox), uris)
    res['dwca_file'] = f'dwca/{record_id}_dwca.zip'
    res['status'] = STATUS_PARTIAL if missing else STATUS_CONVERTED
    return res


def title_taxon(meta, worms):
    """Single Latin name in the title / keywords (single-species datasets), WoRMS-validated."""
    md = meta.get('metadata', {})
    text = ' '.join([md.get('title', '')] + [str(k) for k in md.get('keywords', []) or []])
    cands = list(dict.fromkeys(m.group(0) for m in re.finditer(r'\b[A-Z][a-z]{3,} [a-z]{4,}\b', text)))
    cands = [c for c in cands if not re.match(r'(The|This|Data|Dataset|Black|North|Baltic|Sea|Marine|'
                                              r'Northern|Southern|Western|Eastern|Long|Annual)\b', c)]
    if not cands or len(cands) > 6:
        return None
    hits = [(c, m) for c, m in worms.match(cands).items()
            if m and m.get('match_type') == 'exact'
            and str(m.get('rank', '')).lower() in ('species', 'subspecies')]
    return hits[0][0] if len(hits) == 1 else None


def bbox_temporal(df):
    bbox = temporal = None
    if df is None:
        return None, None
    if 'decimalLatitude' in df.columns and 'decimalLongitude' in df.columns:
        lat = pd.to_numeric(df['decimalLatitude'], errors='coerce').dropna()
        lon = pd.to_numeric(df['decimalLongitude'], errors='coerce').dropna()
        if len(lat) and len(lon):
            bbox = (round(lon.min(), 4), round(lon.max(), 4), round(lat.min(), 4), round(lat.max(), 4))
    if 'eventDate' in df.columns:
        d = df['eventDate'].astype(str).str.extract(r'^(\d{4}(?:-\d{2}(?:-\d{2})?)?)')[0].dropna()
        d2 = df['eventDate'].astype(str).str.extract(r'/(\d{4}(?:-\d{2}(?:-\d{2})?)?)$')[0].dropna()
        allv = pd.concat([d, d2])
        if len(allv):
            temporal = (allv.min(), allv.max())
    return bbox, temporal


REFERENCE_NAME_RE = re.compile(r'code|lookup|look_up|checklist|species.?list|taxa.?list|taxonomy|'
                               r'classification|glossary|dictionary|synonym', re.I)


def is_reference_list(name, df, mapping):
    """True for taxon lists (code lists, checklists, taxonomies) that carry no observations:
    one row per taxon, higher classification columns, no coordinates or stations."""
    cols = [c for c, t in mapping.items() if t == 'scientificName']
    if not cols or len(df) < 50:
        return False
    terms = set(mapping.values())
    if terms & {'decimalLatitude', 'decimalLongitude', 'locationID', 'locality', 'eventID',
                'individualCount', 'organismQuantity'}:
        return False
    taxa = df[cols[0]].dropna().astype(str).str.strip()
    if len(taxa) < 50 or taxa.nunique() / len(taxa) < 0.9:
        return False
    numeric = sum(pd.to_numeric(df[c], errors='coerce').notna().mean() > 0.8 for c in df.columns)
    if numeric >= 5:  # taxon x sample matrix (counts, reads, biomass): real data
        return False
    higher = terms & {'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'vernacularName'}
    return bool(REFERENCE_NAME_RE.search(str(name))) or (len(higher) >= 2 and len(taxa) >= 500)


def convert_dataset(item, session, worms, cache_dir, out_dir, token='', mapping_dir=None,
                    large_files=False):
    record_id = record_id_of(item)
    res = new_result(item, record_id)
    if str(item.get('Access conditions', '')).lower() != 'open':
        res['issues'].append(f'Access conditions: {item.get("Access conditions")}')
        return res
    try:
        meta = zenodo_record(session, record_id, cache_dir, token)
    except Exception as e:
        res['issues'].append(f'Zenodo API unavailable, retry later: {e}')
        return res
    files = [f for f in file_entries(meta) if f['key']
             and os.path.splitext(f['key'])[1].lower() in TABULAR_EXT]
    if not files:
        res['issues'].append('no tabular files in the record')
        return res
    small = [f for f in files if (f['size'] or 0) <= MAX_FILE_MB * 1024 * 1024]
    large = [f for f in files if (f['size'] or 0) > MAX_FILE_MB * 1024 * 1024]
    if large and large_files:
        small += sorted(large, key=lambda f: f['size'])
    elif not small:
        res['issues'].append(f'all tabular files larger than {MAX_FILE_MB} MB (use --large-files)')
        return res
    elif large:
        res['issues'].append(f'{len(large)} file(s) > {MAX_FILE_MB} MB skipped (use --large-files)')
    prio = {'.csv': 0, '.tsv': 0, '.xlsx': 1, '.xls': 1, '.ods': 1, '.txt': 2, '.tab': 2, '.zip': 3}
    small.sort(key=lambda f: (prio.get(os.path.splitext(f['key'])[1].lower(), 9), f['size'] or 0))
    small = small[:MAX_FILES_PER_DATASET]

    override = load_override(mapping_dir, record_id)
    rec_dir = os.path.join(cache_dir, record_id)
    extract_root = os.path.join(rec_dir, '_dwc_extract')
    tables, formats = [], set()
    for f in small:
        dest = os.path.join(rec_dir, re.sub(r'[\\/:*?"<>|]', '_', f['key']))
        try:
            if (f['size'] or 0) > MAX_FILE_MB * 1024 * 1024:
                dest, note = fetch_large(session, f, rec_dir, extract_root, token)
                res['issues'].append(note)
            else:
                download(session, f['url'], dest, token)
        except Exception as e:
            res['issues'].append(f'{f["key"]}: download failed ({str(e)[:160]})')
            continue
        res['source_files'].append(f['key'])
        ext = os.path.splitext(f['key'])[1].lower().lstrip('.')
        formats.add(ext)
        if ext == 'zip':
            try:
                d = os.path.join(extract_root, re.sub(r'[^\w.-]', '_', f['key']))
                if not os.path.isdir(d):
                    safe_extract(dest, d)
                droot = find_dwca_root(d)
                if droot:
                    res['input_format'] = 'dwca'
                    res['layout'] = 'dwca'
                    try:
                        return convert_dwca_source(droot, record_id, meta, res, worms, out_dir)
                    except Exception as e:
                        res['status'] = STATUS_FAILED
                        res['issues'].append(f'DwC-A could not be read: {e}')
                        return res
            except Exception as e:
                res['issues'].append(f'{f["key"]}: {e}')
                continue
        try:
            tables += read_tables(dest, extract_root)
        except Exception as e:
            res['issues'].append(f'{f["key"]}: unreadable ({str(e)[:120]})')
    res['input_format'] = '|'.join(sorted(formats))
    if not res['source_files']:
        res['issues'].append('no file could be downloaded (Zenodo unavailable?), retry later')
        return res
    if not tables:
        res['status'] = STATUS_FAILED
        res['issues'].append('no readable data tables found')
        return res

    description = re.sub(r'<[^<]+?>', ' ', meta.get('metadata', {}).get('description', '') or '')
    hint = f'{item.get("Data title", "")} {item.get("Type of data", "")} {description}'
    context = f'{item.get("Taxa covered", "")} {item.get("Data title", "")}'

    # station lookup tables: location key + coordinates, no taxa
    station_coords = {}
    for name, df in tables:
        m = map_columns(list(df.columns), override.get(name) or override.get('*'))
        terms = set(m.values())
        if 'scientificName' not in terms and {'decimalLatitude', 'decimalLongitude'} <= terms \
                and terms & {'locationID', 'locality'}:
            key_col = next(c for c, t in m.items() if t in ('locationID', 'locality'))
            lat_c = next(c for c, t in m.items() if t == 'decimalLatitude')
            lon_c = next(c for c, t in m.items() if t == 'decimalLongitude')
            for _, r in df.iterrows():
                la, lo = parse_coordinate(r[lat_c], True), parse_coordinate(r[lon_c], False)
                if la is not None and lo is not None:
                    station_coords[str(r[key_col]).strip()] = (la, lo)

    global TAXON_VALIDATOR
    tables_kept = []
    for name, df in tables:
        if is_reference_list(name, df, map_columns(list(df.columns), override.get(name) or override.get('*'))):
            res['issues'].append(f'{name}: taxon reference list (codes / checklist), not occurrences - skipped')
        else:
            tables_kept.append((name, df))
    tables = tables_kept
    TAXON_VALIDATOR = lambda names: {n for n, m in worms.match(names).items() if m}
    default_taxon = title_taxon(meta, worms)
    events, occs, emofs, layouts = {}, [], [], set()
    for name, df in tables:
        m = map_columns(list(df.columns), override.get(name) or override.get('*'))
        built = build_from_table(df, m, record_id, name, context, hint, default_taxon)
        if not built:
            continue
        if built.get('default_taxon'):
            res['issues'].append(f'{name}: no taxon column - scientificName "{default_taxon}" '
                                 f'taken from the dataset title')
        layouts.add(built['layout'])
        for c, t in built['mapping'].items():
            if not c.startswith('__'):
                t = {'_depth': 'minimum/maximumDepthInMeters', '_latStart': 'decimalLatitude (start)',
                     '_lonStart': 'decimalLongitude (start)'}.get(t, t)
                key = c if res['mapped_fields'].get(c, t) == t else f'{name}:{c}'
                res['mapped_fields'][key] = t
        used = set(built['mapping'])
        for c in df.columns:
            if c not in used and not measurement_kind(c) and not abiotic_kind(c) \
                    and not (built['layout'] == 'wide' and LATIN_COL_RE.match(str(c))):
                if c not in res['unmapped_columns']:
                    res['unmapped_columns'].append(c)
        st = built['stats']
        if st['bad_dates']:
            res['issues'].append(f'{name}: {st["bad_dates"]} unparseable dates')
        if st['bad_coords']:
            res['issues'].append(f'{name}: {st["bad_coords"]} invalid coordinates')
        for eid, ev in built['events'].items():
            if eid in events:
                for k, v in ev.items():
                    if v != '' and events[eid].get(k, '') == '':
                        events[eid][k] = v
            else:
                events[eid] = ev
        occs += built['occ']
        emofs += built['emof']
        if len(occs) >= MAX_OCCURRENCES:
            res['issues'].append(f'truncated at {MAX_OCCURRENCES} occurrences')
            break

    res['unmapped_columns'] = res['unmapped_columns'][:60]
    if not occs:
        res['status'] = STATUS_FAILED
        res['issues'].append('no taxon column / species matrix detected in any table')
        return res

    if station_coords:
        filled = 0
        for ev in events.values():
            if ev.get('decimalLatitude', '') == '':
                key = str(ev.get('locationID') or ev.get('locality') or '').strip()
                if key in station_coords:
                    ev['decimalLatitude'], ev['decimalLongitude'] = station_coords[key]
                    filled += 1
        if filled:
            res['issues'].append(f'coordinates of {filled} events taken from a station table')

    res['layout'] = 'wide' if 'wide' in layouts else 'long'
    ev_df = tidy_frame(list(events.values()), EVENT_TERMS, 'eventID')
    has_event_info = any(fill_rate(ev_df, t) > 0 for t in
                         ('eventDate', 'decimalLatitude', 'locationID', 'locality'))
    ev_df['geodeticDatum'] = ['EPSG:4326' if str(la).strip() != '' else ''
                              for la in ev_df.get('decimalLatitude', pd.Series([''] * len(ev_df)))]
    occ_df = pd.DataFrame(occs).fillna('')
    emof_df = pd.DataFrame(emofs).fillna('') if emofs else None
    # Tables repeated in several files / sheets produce identical rows: keep one copy.
    n0 = len(occ_df)
    occ_df = occ_df.astype(str).drop_duplicates().reset_index(drop=True)
    if emof_df is not None:
        emof_df = emof_df.astype(str).drop_duplicates().reset_index(drop=True)
    if len(occ_df) < n0:
        res['issues'].append(f'{n0 - len(occ_df)} duplicate occurrence rows removed (repeated tables)')
    dup = occ_df['occurrenceID'].duplicated(keep=False)
    if dup.any():
        k = occ_df.groupby('occurrenceID').cumcount()
        occ_df.loc[dup & (k > 0), 'occurrenceID'] = occ_df['occurrenceID'] + '_dup' + k.astype(str)
        res['issues'].append(f'{int((dup & (k > 0)).sum())} conflicting occurrenceIDs made unique')
    occ_df = apply_worms(occ_df, worms, res)

    if has_event_info and len(ev_df) < len(occ_df):
        res['core'] = 'Event'
        core_df = ev_df
        occ_out = tidy_frame(occ_df.to_dict('records'), OCC_TERMS, 'occurrenceID')
        occ_out['id'] = occ_out['eventID']
        exts = [('occurrence.txt', ROW_OCCURRENCE, occ_out)]
        if emof_df is not None:
            em = tidy_frame(emof_df.to_dict('records'), EMOF_TERMS, 'eventID')
            em['id'] = em['eventID']
            em = em.drop(columns=['eventID'])
            exts.append(('extendedmeasurementorfact.txt', ROW_EMOF, em))
        core = ('event.txt', ROW_EVENT, core_df)
        loc_df = ev_df
    else:
        res['core'] = 'Occurrence'
        merged = occ_df.merge(ev_df.drop(columns=['id']), on='eventID', how='left') if len(ev_df) else occ_df
        if not has_event_info:
            merged = merged.drop(columns=[c for c in ['eventID'] if c in merged.columns])
        core_df = tidy_frame(merged.to_dict('records'), OCC_TERMS + EVENT_TERMS, 'occurrenceID')
        core = ('occurrence.txt', ROW_OCCURRENCE, core_df)
        exts = []
        if emof_df is not None:
            em = emof_df[emof_df['occurrenceID'] != ''].copy()
            if len(em):
                em = tidy_frame(em.to_dict('records'), EMOF_TERMS, 'occurrenceID')
                em = em.drop(columns=[c for c in ['eventID'] if c in em.columns])
                exts.append(('extendedmeasurementorfact.txt', ROW_EMOF, em))
            dropped = len(emof_df) - len(em)
            if dropped:
                res['issues'].append(f'{dropped} event-level measurements dropped (no events)')
        loc_df = core_df
        occ_df = core_df

    res['n_events'] = len(ev_df) if res['core'] == 'Event' else 0
    res['n_occurrences'] = len(occ_df)
    res['n_emof'] = sum(len(e[2]) for e in exts if e[1] == ROW_EMOF)

    missing = []
    for t in REQUIRED_EVENT[1:]:
        if fill_rate(loc_df, t) < REQUIRED_FILL:
            missing.append(t)
    for t in REQUIRED_OCC:
        thr = SCINAMEID_FILL if t == 'scientificNameID' else REQUIRED_FILL
        if fill_rate(occ_df, t) < thr:
            missing.append(t)
    res['missing_required'] = missing
    if res['n_emof'] == 0:
        res['issues'].append('no quantitative measurements (presence-only data)')

    bbox, temporal = bbox_temporal(loc_df)
    zip_path = os.path.join(out_dir, 'dwca', f'{record_id}_dwca.zip')
    os.makedirs(os.path.dirname(zip_path), exist_ok=True)
    write_dwca(zip_path, core, exts, eml_xml(meta, record_id, temporal, bbox))
    problems = verify_dwca(zip_path)
    if problems:
        res['issues'] += [f'archive check: {p}' for p in problems]
    res['dwca_file'] = f'dwca/{record_id}_dwca.zip'
    res['status'] = STATUS_PARTIAL if (missing or problems) else STATUS_CONVERTED
    return res


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------
def summarise(results):
    c = lambda s: sum(1 for r in results if r['status'] == s)
    return {'total': len(results), 'converted': c(STATUS_CONVERTED), 'partial': c(STATUS_PARTIAL),
            'failed': c(STATUS_FAILED), 'not_attempted': c(STATUS_NOT_ATTEMPTED),
            'total_events': sum(r['n_events'] for r in results),
            'total_occurrences': sum(r['n_occurrences'] for r in results),
            'total_emof': sum(r['n_emof'] for r in results)}


def write_reports(results, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    payload = {'generated': dt.datetime.now().isoformat(timespec='seconds'), 'tool': TOOL,
               'summary': summarise(results), 'datasets': results}
    with open(os.path.join(out_dir, 'conversion_results.json'), 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    flat = []
    for r in results:
        row = {k: v for k, v in r.items() if k not in ('mapped_fields',)}
        for k in ('source_files', 'unmapped_columns', 'missing_required', 'issues'):
            row[k] = '; '.join(map(str, r[k]))
        row['mapped_fields'] = '; '.join(f'{a} -> {b}' for a, b in r['mapped_fields'].items())
        flat.append(row)
    df = pd.DataFrame(flat)
    df.to_csv(os.path.join(out_dir, 'conversion_results.csv'), index=False, encoding='utf-8-sig')
    maps = [{'record_id': r['record_id'], 'title': r['title'], 'source column': a, 'DwC term': b}
            for r in results for a, b in r['mapped_fields'].items()]
    issues = [{'record_id': r['record_id'], 'title': r['title'], 'status': r['status'], 'issue': i}
              for r in results for i in r['issues'] + [f'missing required: {m}' for m in r['missing_required']]]
    summ = pd.DataFrame([{'metric': k, 'value': v} for k, v in payload['summary'].items()] +
                        [{'metric': 'generated', 'value': payload['generated']},
                         {'metric': 'tool', 'value': TOOL}])
    with pd.ExcelWriter(os.path.join(out_dir, 'Conversion_Report_DwC.xlsx'), engine='openpyxl') as w:
        summ.to_excel(w, sheet_name='Summary', index=False)
        df.to_excel(w, sheet_name='Datasets', index=False)
        pd.DataFrame(maps, columns=['record_id', 'title', 'source column', 'DwC term']).to_excel(
            w, sheet_name='Field mappings', index=False)
        pd.DataFrame(issues, columns=['record_id', 'title', 'status', 'issue']).to_excel(
            w, sheet_name='Issues', index=False)
    return payload


# ---------------------------------------------------------------------------
# Selection + main
# ---------------------------------------------------------------------------
def select_datasets(inventory, ids=None, qc_results=None, include_review=False):
    by_id = {record_id_of(r): r for r in inventory}
    if ids:
        return [by_id[i] for i in ids if i in by_id]
    if qc_results:
        with open(qc_results, encoding='utf-8') as fh:
            qc = json.load(fh)
        ok = [d for d in qc.get('datasets', []) if d.get('qc_status') in ('Pass', 'Warnings')]
        ok.sort(key=lambda d: -(d.get('data_score') or 0))
        return [by_id[str(d['record_id'])] for d in ok if str(d.get('record_id')) in by_id]
    sel = []
    for r in inventory:
        if r.get('Status') != 'Included' and not include_review:
            continue
        exts = set(str(r.get('File type', '')).split(', '))
        if r.get('Access conditions') == 'Open' and exts & TABULAR_EXT:
            sel.append(r)

    def prio(r):
        exts = set(str(r.get('File type', '')).split(', '))
        black = 'Black Sea' in r.get('Geographic coverage', '')
        direct = bool(exts & {'.csv', '.tsv', '.xlsx', '.xls'})
        text = (r.get('Data title', '') + ' ' + r.get('Review notes', '')).lower()
        occ = bool(re.search(r'occurrence|abundance|survey|monitoring|distribution|biomass', text))
        return (not black, not direct, not occ, len(exts))
    sel.sort(key=prio)
    return sel


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description='Convert Zenodo marine biology datasets to Darwin Core Archives')
    ap.add_argument('--inventory', default=os.path.join(here, '..', '..', 'Web', 'zenodo_inventory', 'data.json'),
                    help='Task 1 inventory JSON (data.json)')
    ap.add_argument('--cache', default=os.path.join(here, '..', 'zenodo_cache'), help='download cache folder')
    ap.add_argument('--out', default=os.path.join(here, 'output'), help='output folder')
    ap.add_argument('--ids', nargs='*', help='Zenodo record ids to convert')
    ap.add_argument('--limit', type=int, help='max number of datasets')
    ap.add_argument('--mapping-dir', help='folder with <record_id>.json/.yaml column mapping overrides')
    ap.add_argument('--qc-results', help='Task 2 qc_results.json: convert only Pass/Warnings datasets')
    ap.add_argument('--include-review', action='store_true', help='also take "Manual review" records')
    ap.add_argument('--retry-passes', type=int, default=2,
                    help='extra passes over datasets skipped because Zenodo was unavailable')
    ap.add_argument('--retry-wait', type=int, default=600, help='seconds to wait before a retry pass')
    ap.add_argument('--merge', action='store_true',
                    help='update the existing conversion_results.json in --out instead of replacing it '
                         '(for re-runs of a subset with --ids)')
    ap.add_argument('--large-files', action='store_true',
                    help=f'also process files > {MAX_FILE_MB} MB: zips are read remotely (HTTP Range, '
                         f'tables only), text tables from their first {LARGE_TEXT_HEAD_MB} MB')
    args = ap.parse_args(argv)

    token = os.environ.get('ZENODO_TOKEN', '')
    with open(args.inventory, encoding='utf-8') as fh:
        inventory = json.load(fh)
    todo = select_datasets(inventory, args.ids, args.qc_results, args.include_review)
    if args.limit:
        todo = todo[:args.limit]
    os.makedirs(args.cache, exist_ok=True)
    os.makedirs(os.path.join(args.out, 'dwca'), exist_ok=True)
    session = robust_session()
    worms = Worms(os.path.join(args.cache, '_worms_cache.json'), session)
    results = {}
    previous = {}
    prev_path = os.path.join(args.out, 'conversion_results.json')
    if args.merge and os.path.exists(prev_path):
        with open(prev_path, encoding='utf-8') as fh:
            previous = {r['record_id']: r for r in json.load(fh)['datasets']}

    def merged():
        return list({**previous, **results}.values())

    print(f'{TOOL}: {len(todo)} datasets')

    def unavailable(res):
        return res['status'] == STATUS_NOT_ATTEMPTED and any('retry later' in i for i in res['issues'])

    queue = todo
    for attempt in range(args.retry_passes + 1):
        if attempt:
            print(f'Retry pass {attempt}: {len(queue)} datasets, Zenodo was unavailable; '
                  f'waiting {args.retry_wait}s', flush=True)
            time.sleep(args.retry_wait)
        for k, item in enumerate(queue, 1):
            rid = record_id_of(item)
            t0 = time.time()
            try:
                res = convert_dataset(item, session, worms, args.cache, args.out, token, args.mapping_dir,
                                      args.large_files)
            except Exception as e:
                res = new_result(item, rid)
                res['status'] = STATUS_FAILED
                res['issues'].append(f'unexpected error: {type(e).__name__}: {str(e)[:200]}')
            results[rid] = res
            print(f'[{k}/{len(queue)}] {rid} {res["status"]:<13} ev={res["n_events"]} occ={res["n_occurrences"]} '
                  f'emof={res["n_emof"]} ({time.time() - t0:.0f}s) {item.get("Data title", "")[:60]}', flush=True)
            if k % 10 == 0:
                write_reports(merged(), args.out)
        queue = [it for it in queue if unavailable(results[record_id_of(it)])]
        if not queue:
            break
    payload = write_reports(merged(), args.out)
    print(json.dumps(payload['summary'], indent=2))
    return payload


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
