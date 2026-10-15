"""
Data types of the Zenodo records and export of the non-biological tables (Task 3, other data).

The converter (dwc_converter.py) turns the biological tables into Darwin Core Archives. The same
records also hold chemistry, physical oceanography, bathymetry, sediment, tracking, acoustic,
genetic and other data. This script

  1. lists every file of every record selected by the converter (open, Included), including the
     members of zip archives (listed remotely through HTTP Range when the zip is not cached);
  2. classifies each file / table into a data type (file extension, file and sheet names,
     column headers; the record title when the names say nothing);
  3. exports every non-biological table as a standardised CSV (UTF-8) and XLSX:
       - original columns kept, empty rows / columns dropped;
       - ISO 8601 dates and decimal-degree coordinates added when date / coordinate columns
         are recognised (date_ISO8601, latitude_decimal, longitude_decimal);
       - a "Columns" sheet (original header, unit, parameter group, Darwin Core term, fill rate)
         and a "Source" sheet (record, DOI, license, creators, original file);
     small NetCDF files are flattened to a table as well;
  4. writes output/data_types.json and Data_Types_Report.xlsx with the counts per data type:
       found      files (zip members included) of that type,
       processed  tables read and standardised (biology: tables in a Darwin Core Archive),
       saved      CSV + XLSX files written (biology: the DwC-A).

    python data_types.py                 # all records of conversion_results.json
    python data_types.py --ids 3998049   # some records
"""
import argparse
import gc
import json
import os
import re
import shutil
import sys
import time
import zipfile
from datetime import date

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dwc_converter import (TOOL, MAX_FILE_MB, MAX_ROWS_PER_TABLE, HttpRangeFile, detect_layout,  # noqa: E402
                           file_entries, header_unit, map_columns, normalise_header, parse_coordinate,
                           parse_date, read_tables, robust_session, to_number)

TYPES = {  # key: (label, colour on the web page)
    'biology': ('Biological occurrences / abundance', '#2e7d32'),
    'biology_other': ('Biological measurements without taxon column', '#66bb6a'),
    'traits': ('Traits / morphometrics / diet', '#7cb342'),
    'chemistry': ('Chemistry / biogeochemistry', '#c62828'),
    'physics': ('Physical oceanography / meteorology', '#1565c0'),
    'bathymetry': ('Bathymetry / topography', '#283593'),
    'sediment': ('Sediment / geology', '#8d6e63'),
    'tracking': ('Tracking / telemetry', '#ef6c00'),
    'acoustic': ('Acoustic', '#6a1b9a'),
    'genetic': ('Genetic / eDNA / sequences', '#00897b'),
    'imagery': ('Imagery / video', '#f9a825'),
    'spatial': ('GIS / rasters / model output', '#546e7a'),
    'fisheries': ('Fisheries / effort / human activities', '#ad1457'),
    'binary': ('R / MATLAB / other binary data', '#9e9e9e'),
    'documentation': ('Documentation / code', '#bdbdbd'),
    'other': ('Other tables', '#757575'),
}
EXPORT_TYPES = {'biology_other', 'traits', 'chemistry', 'physics', 'bathymetry', 'sediment', 'tracking', 'acoustic',
                'spatial', 'fisheries', 'other', 'genetic'}

EXT_TYPES = {
    'genetic': {'.fastq', '.fq', '.fasta', '.fas', '.fa', '.fna', '.faa', '.bam', '.sam', '.vcf', '.trees',
                '.tre', '.newick', '.nwk', '.nex', '.nexus', '.phy', '.aln', '.biom', '.gb', '.gbk'},
    'imagery': {'.jpg', '.jpeg', '.png', '.bmp', '.gif', '.heic', '.mp4', '.avi', '.mov', '.mkv', '.mts'},
    'acoustic': {'.wav', '.flac', '.mp3', '.aif', '.aiff', '.cp1', '.cp3', '.fp1', '.fp3'},
    'spatial': {'.nc', '.nc4', '.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj', '.gpkg', '.geojson', '.kmz',
                '.tif', '.tiff', '.grd', '.gri', '.asc', '.img', '.hdf', '.he5', '.grib', '.grb'},
    'binary': {'.rdata', '.rda', '.rds', '.mat', '.parquet', '.h5', '.hdf5', '.npy', '.npz', '.sav', '.pkl',
               '.feather', '.sqlite', '.db', '.mdb', '.accdb'},
    'documentation': {'.pdf', '.doc', '.docx', '.md', '.rtf', '.html', '.htm', '.r', '.rmd', '.py', '.m',
                      '.ipynb', '.rproj', '.json', '.out', '.log', '.xml', '.yml', '.yaml', '.sh', '.jl',
                      '.c', '.cpp', '.js', '.css', '.svg', '.ppt', '.pptx', '.odt', '.tex', '.bib',
                      '.qmd', '.sas', '.do', '.stan', '.bug', '.jags', '.cfg', '.ini', '.toml', '.lock'},
}
TABLE_EXT = {'.csv', '.tsv', '.txt', '.tab', '.dat', '.xlsx', '.xls', '.xlsm', '.ods'}
DOC_NAME_RE = re.compile(r'readme|licen[cs]e|citation|changelog|manifest|checksum|md5|dictionary|codebook|'
                         r'legend|metadata|description|abbreviation|column.?def|variables?.?list|glossary|'
                         r'lookup|look.?up|vocabular|changelog', re.I)

# keyword patterns on the normalised file / sheet name and column headers
TYPE_RE = {
    'biology_other': r'biomass|abundance|\bdensity\b|zooplankton|phytoplankton|benthos|benthic|macrofauna|'
                     r'meiofauna|plankton|\bfish|species richness|diversity|shannon|\bcpue\b|\bcount\b|'
                     r'\bcounts\b|individuals|\bind\b|cover|occurrence|presence|sighting|seagrass|macroalga',
    'chemistry': r'nutrient|nitrat|nitrit|phosphat|silicat|silic|ammoni|\bnh4|\bno3|\bno2|\bpo4|\bsio[24]?\b|'
                 r'oxygen|\bo2\b|\bdo\b|\bph\b|alkalin|\bdic\b|\bta\b|pco2|\bco2|carbonate|\btoc\b|\bdoc\b|'
                 r'\bpoc\b|\bpon\b|\btn\b|\btp\b|chlorophyl|\bchl|pigment|metal|\bcu\b|\bzn\b|\bpb\b|\bcd\b|'
                 r'\bhg\b|mercury|\bfe\b|\bmn\b|\bni\b|\bcr\b|\bas\b|\bpah|\bpcb|pollut|contamin|isotop|'
                 r'd13c|d15n|d18o|13c|15n|lipid|fatty acid|biomarker|chemi|methane|\bch4\b|sulfide|sulphide|'
                 r'\bh2s\b|nitrogen|phosphor|carbon|\bc n\b|chloride|calcium|magnesium|toxin|microplastic',
    'physics': r'\bctd\b|temperat|\btemp\b|salinit|\bsal\b|\bpsu\b|conductiv|densit|sigma|turbid|current|'
               r'velocity|\bwave|tide|tidal|sea level|\bpar\b|irradian|\blight\b|\bwind|\bsst\b|mixed layer|'
               r'\bmld\b|\bice\b|pressure|meteo|weather|precip|air temp|stratific|upwelling|secchi|'
               r'transparen|fluoresc|beam|attenuat',
    'bathymetry': r'bathym|multibeam|backscatter|\bdem\b|elevation|topograph|seafloor|sea bed|seabed|slope|'
                  r'rugosity|depth model|sounding|bottom depth|water depth',
    'sediment': r'sediment|grain size|granulom|\bmud\b|\bsand\b|gravel|silt|\bclay\b|porosit|organic matter|'
                r'\bloi\b|loss on|\bcore\b|cores|geolog|mineral|\bd50\b|median grain|sorting',
    'acoustic': r'acoust|hydrophon|\bc ?pod\b|\bf ?pod\b|click|\bdpm\b|sound|noise|\bspl\b|echo|sonar|'
                r'vocali|whistle|detection positive',
    'tracking': r'\bgps\b|argos|track|telemetr|\btag\b|\btags\b|\bptt\b|\bdive|movement|migrat|geolocat|'
                r'acceleromet|\btdr\b|location class|\blc\b|haul out|foraging trip',
    'traits': r'trait|morpho|\blength\b|\bweight\b|\bsize\b|biometr|\bage\b|growth|condition|fecund|'
              r'physiolog|respirat|metabol|feeding|\bdiet\b|stomach|\bgut\b|prey|body mass|\btl\b|\bsl\b|'
              r'wingspan|otolith|egg|clutch|breeding|reproduct',
    'fisheries': r'\bcatch|landing|effort|vessel|fleet|fisher|logbook|bycatch|\bgear\b|aquacult|shipping|'
                 r'tourism|\bais\b|trawling intensity|fishing',
    'spatial': r'raster|\bmodel|simulat|predict|scenario|projection|\bgrid\b|habitat suitab|\bsdm\b|'
               r'covariate|environmental layer',
}
TYPE_RE = {k: re.compile(v) for k, v in TYPE_RE.items()}

PARAM_GROUPS = [  # (regex on normalised header, parameter group) for the Columns sheet
    (r'\b(lat|latitude|lon|long|longitude|lng|x|y|easting|northing|coord)', 'Position'),
    (r'\b(date|time|year|month|day|datetime|timestamp|hour|doy|julian)\b', 'Date / time'),
    (r'\b(depth|pressure|press|pres|altitude|elevation|bathym)', 'Depth / pressure / elevation'),
    (r'nitrat|nitrit|phosphat|silic|ammoni|\bnh4|\bno3|\bno2|\bpo4|\bsio|\bdin\b|\bdip\b|\bsrp\b', 'Nutrients'),
    (r'oxygen|\bo2\b|\bdo\b', 'Dissolved oxygen'),
    (r'\bph\b|alkalin|\bdic\b|pco2|\bco2|carbonate|\bta\b', 'Carbonate system'),
    (r'chlorophyl|\bchl|pigment|phaeo|fluoresc', 'Pigments / fluorescence'),
    (r'\btoc\b|\bdoc\b|\bpoc\b|\bpon\b|\btn\b|\btp\b|organic|carbon|nitrogen|\bc n\b|\bloi\b', 'Organic matter / elemental'),
    (r'd13c|d15n|d18o|13c|15n|isotop', 'Stable isotopes'),
    (r'\bcu\b|\bzn\b|\bpb\b|\bcd\b|\bhg\b|\bfe\b|\bmn\b|\bni\b|\bcr\b|\bas\b|metal|mercury', 'Metals'),
    (r'\bpah|\bpcb|pollut|contamin|pesticid|microplastic|toxin', 'Contaminants'),
    (r'temperat|\btemp\b|\bsst\b', 'Temperature'),
    (r'salinit|\bsal\b|\bpsu\b|conductiv', 'Salinity / conductivity'),
    (r'turbid|secchi|transparen|attenuat|\bpar\b|irradian|light', 'Light / turbidity'),
    (r'current|velocity|\bwave|tide|tidal|\bwind', 'Currents / waves / wind'),
    (r'grain|\bmud\b|\bsand\b|gravel|silt|\bclay\b|\bd50\b|porosit', 'Grain size / sediment'),
    (r'\blength\b|\bweight\b|\bsize\b|\bmass\b|\bage\b|\btl\b|\bsl\b', 'Size / mass / age'),
    (r'\b(species|taxon|taxa|genus|family|scientific)', 'Taxonomy'),
    (r'\b(id|code|station|site|sample|replicate|transect|individual|tag|ptt)\b', 'Identifier / sampling unit'),
]
PARAM_GROUPS = [(re.compile(rx), g) for rx, g in PARAM_GROUPS]

MAX_ZIP_LIST = 5000
MAX_NC_MB = MAX_FILE_MB
XLSX_MAX_ROWS = 1_048_575


def ext_of(name):
    n = name.lower()
    if n.endswith('.gz') or n.endswith('.bz2'):
        n = n.rsplit('.', 1)[0]
    return os.path.splitext(n)[1]


def name_text(*parts):
    return ' ' + ' '.join(normalise_header(re.sub(r'[_\-./]+', ' ', str(p))) for p in parts if p) + ' '


def keyword_scores(text):
    return {k: len(rx.findall(text)) for k, rx in TYPE_RE.items()}


def classify_name(name, record_text=''):
    """Data type of a file from its name / extension (tables: preliminary)."""
    ext = ext_of(name)
    base = os.path.basename(name)
    nt = name_text(base)
    if ext in ('.tif', '.tiff') and re.search(r'photo|image|img|picture|scan|plate|fig', nt):
        return 'imagery'
    if ext in ('.tif', '.tiff', '.asc', '.grd', '.gri', '.xyz') and TYPE_RE['bathymetry'].search(nt):
        return 'bathymetry'
    if ext in ('.kml', '.gpx'):
        return 'tracking' if TYPE_RE['tracking'].search(nt) else 'spatial'
    if ext == '.cnv':
        return 'physics'
    for t, exts in EXT_TYPES.items():
        if ext in exts:
            return t
    if ext in TABLE_EXT or ext == '.zip':
        if DOC_NAME_RE.search(base) and ext not in ('.xlsx', '.xls'):
            return 'documentation'
        s = keyword_scores(nt)
        best = max(s, key=s.get)
        return best if s[best] else None
    return 'other' if ext else 'documentation'


def classify_table(name, df, record_text):
    """Data type of a table from its name and columns; biology when it holds taxa."""
    try:
        mapping = map_columns(list(df.columns))
    except Exception:
        mapping = {}
    terms = set(mapping.values())
    layout, _ = (detect_layout(df, mapping) if len(df.columns) else (None, []))
    if DOC_NAME_RE.search(name.split('::')[-1]) and len(df) < 500:
        return 'documentation', mapping
    if len(df) < 3:
        return 'documentation', mapping                 # one-row lookup / metadata sheets
    s = keyword_scores(name_text(name) + ' ' + ' '.join(name_text(c) for c in df.columns))
    rs = keyword_scores(record_text)
    s = {k: 3 * v + min(rs[k], 2) if v else 0 for k, v in s.items()}
    best = max(s, key=s.get)
    if layout or 'scientificName' in terms:
        # a species table that also has a few environmental columns is still biology,
        # a table of traits / diet per species is not an occurrence table
        if s['traits'] >= 3 and s['traits'] > s[best] - 1 and not any(t in terms for t in
                                                                       ('eventDate', 'decimalLatitude')):
            return 'traits', mapping
        return 'biology', mapping
    return (best if s[best] else 'other'), mapping


# ---------------------------------------------------------------------------
# standardised export
# ---------------------------------------------------------------------------
def param_group(col):
    n = normalise_header(col)
    return next((g for rx, g in PARAM_GROUPS if rx.search(n)), '')


def standardise(df, mapping, add_columns=True):
    df = df.copy()
    df.columns = [re.sub(r'\s+', ' ', str(c)).strip() or f'column_{i + 1}' for i, c in enumerate(df.columns)]
    seen = {}
    cols = []
    for c in df.columns:                                # unique headers
        seen[c] = seen.get(c, 0) + 1
        cols.append(c if seen[c] == 1 else f'{c}_{seen[c]}')
    mapping = {cols[list(df.columns).index(k)] if k in df.columns else k: v for k, v in mapping.items()}
    df.columns = cols
    df = df.replace(r'^\s*$', None, regex=True).dropna(how='all').dropna(axis=1, how='all')
    added = []
    if not add_columns:
        return df, mapping, added
    inv = {v: k for k, v in mapping.items() if k in df.columns}
    if 'eventDate' in inv:
        iso = df[inv['eventDate']].map(lambda v: parse_date(v) if v is not None else None)
        same = (iso.astype(str) == df[inv['eventDate']].astype(str)).mean()
        if iso.notna().mean() >= 0.5 and same < 0.95:
            df['date_ISO8601'] = iso
            added.append(('date_ISO8601', f'ISO 8601 date parsed from "{inv["eventDate"]}"'))
    for term, new, is_lat in (('decimalLatitude', 'latitude_decimal', True),
                              ('decimalLongitude', 'longitude_decimal', False)):
        if term in inv:
            vals = df[inv[term]].map(lambda v: parse_coordinate(v, is_lat) if v is not None else None)
            same = (vals == df[inv[term]].map(to_number)).mean()
            if vals.notna().mean() >= 0.5 and same < 0.95:
                df[new] = vals
                added.append((new, f'decimal degrees (WGS84 assumed) parsed from "{inv[term]}"'))
    return df, mapping, added


def columns_sheet(df, mapping, added):
    notes = dict(added)
    rows = []
    for c in df.columns:
        s = df[c]
        nn = s.notna().sum()
        num = s.dropna().head(500).map(lambda v: to_number(v) is not None).mean() if nn else 0
        rows.append({'column': c, 'unit (from header)': header_unit(c), 'parameter group': param_group(c),
                     'Darwin Core term': mapping.get(c, ''), 'type': 'numeric' if num >= 0.8 else 'text',
                     'non-empty values': int(nn), 'fill rate %': round(100 * nn / max(len(df), 1), 1),
                     'example': str(s.dropna().iloc[0])[:80] if nn else '',
                     'note': notes.get(c, '')})
    return pd.DataFrame(rows)


def source_sheet(rec, file_name, table_name, dtype, n_rows, truncated):
    md = rec['meta'].get('metadata', {})
    creators = '; '.join(c.get('name', '') for c in md.get('creators', []) or [])
    lic = (md.get('license') or {}).get('id', '')
    rows = [('Zenodo record', f'https://zenodo.org/records/{rec["id"]}'),
            ('Title', md.get('title', '')), ('DOI', f'https://doi.org/{md["doi"]}' if md.get('doi') else ''),
            ('Creators', creators), ('License', lic), ('Publication date', md.get('publication_date', '')),
            ('Original file', file_name), ('Table / sheet', table_name),
            ('Data type', TYPES[dtype][0]), ('Rows', n_rows),
            ('Note', '; '.join(x for x in (
                'XLSX truncated at 1,048,575 rows; the CSV holds all rows' if truncated else '',
                'NoDerivatives license: format change only, no columns added' if no_derivatives(rec) else '')
                if x)),
            ('Processed', f'{date.today().isoformat()} with {TOOL} (data_types.py)'),
            ('Citation', f'Please cite the original dataset: {creators} ({str(md.get("publication_date", ""))[:4]}). '
                         f'{md.get("title", "")}. Zenodo. https://doi.org/{md.get("doi", "")}')]
    return pd.DataFrame(rows, columns=['field', 'value'])


def safe_name(s, n=90):
    s = re.sub(r'[^\w.\-]+', '_', s).strip('_')
    return s[:n] or 'table'


def no_derivatives(rec):
    lic = str(((rec['meta'].get('metadata', {}) or {}).get('license') or {}).get('id', '')).lower()
    return '-nd' in lic


def export_table(rec, file_name, table_name, df, mapping, dtype, out_root):
    # CC BY-ND: change of format only (allowed), no added columns (an adaptation)
    df, mapping, added = standardise(df, mapping, add_columns=not no_derivatives(rec))
    if df.empty or df.shape[1] == 0:
        return None
    d = os.path.join(out_root, dtype, rec['id'])
    os.makedirs(d, exist_ok=True)
    base = safe_name(re.sub(r'(?i)\.(csv|tsv|txt|tab|dat|xlsx|xls|xlsm|ods|nc)(?=::|$)', '',
                            table_name.split('/')[-1]).replace('::', '__'))
    csv_path, xlsx_path = os.path.join(d, base + '.csv'), os.path.join(d, base + '.xlsx')
    df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    truncated = len(df) > XLSX_MAX_ROWS
    with pd.ExcelWriter(xlsx_path, engine='xlsxwriter',
                        engine_kwargs={'options': {'strings_to_urls': False, 'strings_to_formulas': False,
                                                   'strings_to_numbers': True}}) as xw:
        df.head(XLSX_MAX_ROWS).to_excel(xw, sheet_name='Data', index=False)
        columns_sheet(df, mapping, added).to_excel(xw, sheet_name='Columns', index=False)
        source_sheet(rec, file_name, table_name, dtype, len(df), truncated).to_excel(xw, sheet_name='Source',
                                                                                     index=False)
    return {'csv': os.path.relpath(csv_path, out_root).replace('\\', '/'),
            'xlsx': os.path.relpath(xlsx_path, out_root).replace('\\', '/'),
            'rows': int(len(df)), 'columns': int(df.shape[1]),
            'csv_mb': round(os.path.getsize(csv_path) / 1e6, 2), 'xlsx_mb': round(os.path.getsize(xlsx_path) / 1e6, 2),
            'parameters': sorted({param_group(c) for c in df.columns} - {''}),
            'added': [a for a, _ in added]}


MAX_NC_JSON_VALUES = 50_000_000   # values of one variable written to the JSON
MAX_NC_TABLE_ROWS = 100_000       # larger variable groups (gridded fields) stay in the JSON only


def _jsonable(v):
    import numpy as np
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return None if np.isnan(v) else float(v)
    if isinstance(v, np.ndarray):
        return [_jsonable(x) for x in v.tolist()]
    if isinstance(v, (bytes, bytearray)):
        return v.decode('utf-8', 'ignore')
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, float) and v != v:
        return None
    return v if isinstance(v, (int, float, str, bool)) or v is None else str(v)


def _array_json(values):
    """numpy array -> nested lists for JSON: NaN -> null, float32 to 7 significant digits,
    datetimes as ISO 8601 strings."""
    import numpy as np
    a = np.asarray(values)
    if a.dtype.kind == 'M':
        flat = [None if pd.isna(v) else pd.Timestamp(v).isoformat() for v in a.ravel()]
    elif a.dtype.kind == 'f':
        fmt = '%.7g' if a.dtype.itemsize <= 4 else '%.15g'
        flat = [None if s == 'nan' else float(s) for s in np.char.mod(fmt, a.ravel())]
    elif a.dtype.kind in 'iub':
        flat = a.ravel().tolist()
    else:
        flat = [_jsonable(v) for v in a.ravel()]
    if a.ndim <= 1:
        return flat if a.ndim == 1 else flat[0]
    return np.array(flat, dtype=object).reshape(a.shape).tolist()


def netcdf_to_json(path, rec, file_name):
    """NetCDF -> JSON (CF-JSON-like): source, global attributes, dimensions and every variable
    with its dimensions, attributes, shape and data (nested arrays in dimension order)."""
    import xarray as xr
    try:
        ds = xr.open_dataset(path)
    except Exception:
        ds = xr.open_dataset(path, decode_times=False)
    with ds:
        md = rec['meta'].get('metadata', {})
        doc = {'source': {'zenodo_record': f'https://zenodo.org/records/{rec["id"]}', 'title': md.get('title', ''),
                          'doi': md.get('doi', ''), 'license': (md.get('license') or {}).get('id', ''),
                          'file': file_name, 'converted_by': f'{TOOL} (data_types.py)',
                          'converted': date.today().isoformat()},
               'global_attributes': {k: _jsonable(v) for k, v in ds.attrs.items()},
               'dimensions': {str(k): int(v) for k, v in ds.sizes.items()},
               'variables': {}}
        for name, var in list(ds.coords.items()) + list(ds.data_vars.items()):
            v = {'dimensions': [str(d) for d in var.dims], 'shape': [int(s) for s in var.shape],
                 'dtype': str(var.dtype), 'role': 'coordinate' if name in ds.coords else 'data',
                 'attributes': {k: _jsonable(x) for k, x in var.attrs.items()}}
            if var.size <= MAX_NC_JSON_VALUES:
                v['data'] = _array_json(var.values)
            else:
                v['note'] = f'{var.size:,} values: more than {MAX_NC_JSON_VALUES:,}, data not written'
            doc['variables'][str(name)] = v
    return doc


def json_tables(doc, max_rows=MAX_NC_TABLE_ROWS):
    """Long tables built from a NetCDF JSON document, one per group of data variables sharing the
    same dimensions (columns: the dimension coordinates, then the variables). Groups with more
    than max_rows rows (gridded fields) are left in the JSON only."""
    import itertools
    import numpy as np
    variables = doc['variables']
    groups = {}
    for name, v in variables.items():
        if v['role'] == 'data' and 'data' in v:
            groups.setdefault(tuple(v['dimensions']), []).append(name)
    tables, skipped = [], []
    if len(groups) > 1:
        groups.pop((), None)                             # scalar variables: in the JSON only
    for dims, names in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        n = int(np.prod([doc['dimensions'][d] for d in dims])) if dims else 1
        label = '_'.join(dims) or 'scalars'
        if n > max_rows:
            skipped.append(f'{label} ({n:,} cells)')
            continue
        axes = []
        for d in dims:
            c = variables.get(d)
            axes.append(c['data'] if c and 'data' in c and c['dimensions'] == [d]
                        else list(range(doc['dimensions'][d])))
        cols = {d: [] for d in dims}
        for combo in itertools.product(*axes):
            for d, x in zip(dims, combo):
                cols[d].append(x)
        for name in names:
            cols[name] = np.array(variables[name]['data'], dtype=object).reshape(-1).tolist() if dims \
                else [variables[name]['data']]
        tables.append((label, pd.DataFrame(cols)))
    return tables, skipped


def export_netcdf(rec, file_name, path, dtype, out_root):
    doc = netcdf_to_json(path, rec, file_name)
    d = os.path.join(out_root, dtype, rec['id'])
    os.makedirs(d, exist_ok=True)
    base = safe_name(re.sub(r'(?i)\.nc4?$', '', file_name.split('/')[-1]))
    json_path = os.path.join(d, base + '.json')
    with open(json_path, 'w', encoding='utf-8') as fh:
        json.dump(doc, fh, ensure_ascii=False, allow_nan=False, separators=(',', ':'))
    tables, skipped = json_tables(doc)
    rel = lambda p: os.path.relpath(p, out_root).replace('\\', '/')  # noqa: E731
    data_vars = [k for k, v in doc['variables'].items() if v['role'] == 'data']
    out = {'json': rel(json_path), 'json_mb': round(os.path.getsize(json_path) / 1e6, 2),
           'rows': sum(len(t) for _, t in tables), 'columns': sum(t.shape[1] for _, t in tables),
           'variables': len(data_vars), 'dimensions': doc['dimensions'],
           'parameters': sorted({param_group(c) for c in doc['variables']} - {''})}
    if not tables:
        out['note'] = 'NetCDF converted to JSON; gridded field (' + '; '.join(skipped) + '): JSON only'
        return out
    csvs = []
    for i, (name, t) in enumerate(tables):
        p = os.path.join(d, f'{base}.csv' if i == 0 else f'{base}__{safe_name(name, 40)}.csv')
        t.to_csv(p, index=False, encoding='utf-8-sig')
        csvs.append(rel(p))
    xlsx_path = os.path.join(d, base + '.xlsx')
    var_rows = [{'variable': k, 'role': v['role'], 'dimensions': ', '.join(v['dimensions']), 'dtype': v['dtype'],
                 'units': v['attributes'].get('units', ''), 'long_name': v['attributes'].get('long_name', ''),
                 'standard_name': v['attributes'].get('standard_name', ''), 'parameter group': param_group(k)}
                for k, v in doc['variables'].items()]
    with pd.ExcelWriter(xlsx_path, engine='xlsxwriter',
                        engine_kwargs={'options': {'strings_to_urls': False, 'strings_to_formulas': False}}) as xw:
        for i, (name, t) in enumerate(tables[:20]):
            t.to_excel(xw, sheet_name=('Data' if i == 0 else f'Data_{safe_name(name, 20)}')[:31], index=False)
        pd.DataFrame(var_rows).to_excel(xw, sheet_name='Variables', index=False)
        pd.DataFrame([{'attribute': k, 'value': str(v)[:32000]} for k, v in doc['global_attributes'].items()]
                     ).to_excel(xw, sheet_name='Global attributes', index=False)
        source_sheet(rec, file_name, ' + '.join(n for n, _ in tables), dtype, out['rows'], False).to_excel(
            xw, sheet_name='Source', index=False)
    out.update(csv=csvs[0], csv_extra=csvs[1:], xlsx=rel(xlsx_path),
               csv_mb=round(sum(os.path.getsize(os.path.join(out_root, c)) for c in csvs) / 1e6, 2),
               xlsx_mb=round(os.path.getsize(xlsx_path) / 1e6, 2),
               note='NetCDF converted to JSON, then to CSV / XLSX' +
                    (f'; gridded fields JSON only: {"; ".join(skipped)}' if skipped else ''))
    return out


# ---------------------------------------------------------------------------
# records
# ---------------------------------------------------------------------------
def zip_members(session, rid, f, rec_dir, token):
    """Member names of a zip: local file when cached, otherwise its central directory over HTTP."""
    cache = os.path.join(rec_dir, '_zip_listing.json')
    listing = {}
    if os.path.exists(cache):
        with open(cache, encoding='utf-8') as fh:
            listing = json.load(fh)
    if f['key'] in listing:
        return listing[f['key']]
    local = os.path.join(rec_dir, f['key'])
    names = None
    try:
        if os.path.exists(local):
            with zipfile.ZipFile(local) as z:
                names = [(i.filename, i.file_size) for i in z.infolist() if not i.is_dir()]
        elif f.get('url'):
            with zipfile.ZipFile(HttpRangeFile(session, f['url'], f['size'], token)) as z:
                names = [(i.filename, i.file_size) for i in z.infolist() if not i.is_dir()]
            time.sleep(0.5)
    except Exception as e:
        names = {'error': f'{type(e).__name__}: {str(e)[:120]}'}
    if isinstance(names, list):
        names = [n for n in names if not re.search(r'(^|/)(__MACOSX|\.DS_Store|\._)', n[0])][:MAX_ZIP_LIST]
    listing[f['key']] = names
    os.makedirs(rec_dir, exist_ok=True)
    with open(cache, 'w', encoding='utf-8') as fh:
        json.dump(listing, fh)
    return names


def process_record(d, args, session, token):
    rid = str(d['record_id'])
    rec_dir = os.path.join(args.cache, rid)
    try:
        with open(os.path.join(rec_dir, '_record.json'), encoding='utf-8') as fh:
            meta = json.load(fh)
    except OSError:
        return []
    md = meta.get('metadata', {})
    rec = {'id': rid, 'meta': meta, 'nc': 0}
    for t in TYPES:                                       # outputs of an earlier run of this record
        old = os.path.join(args.other_out, t, rid)
        if os.path.isdir(old):
            shutil.rmtree(old)
    record_text = name_text(md.get('title', ''), ' '.join(md.get('keywords', []) or []))
    has_dwca = bool(d.get('dwca_file'))
    out = []
    extract_root = os.path.join(rec_dir, '_dwc_extract')

    def entry(file_name, dtype, size, **kw):
        e = {'record_id': rid, 'title': md.get('title', '')[:200], 'file': file_name, 'type': dtype,
             'ext': ext_of(file_name) or '(none)', 'size_mb': round((size or 0) / 1e6, 3),
             'processed': False, 'saved': False, 'note': ''}
        e.update(kw)
        out.append(e)
        return e

    pending = []

    def handle_tables(file_name, tables, size, prelim):
        if not tables:
            doc = DOC_NAME_RE.search(os.path.basename(file_name))
            entry(file_name, 'documentation' if doc else (prelim or 'other'), size,
                  note='' if doc else 'no readable table')
            return
        for tname, df in tables:
            dtype, mapping = classify_table(tname, df, record_text)
            e = entry(file_name, dtype, size if tname == tables[0][0] else 0, table=tname,
                      rows=int(len(df)), columns=int(df.shape[1]))
            pending.append((e, df, mapping))

    for f in file_entries(meta):
        key, size = f['key'], f['size'] or 0
        ext = ext_of(key)
        local = os.path.join(rec_dir, key)
        prelim = classify_name(key, record_text)
        if ext == '.zip':
            members = zip_members(session, rid, f, rec_dir, token)
            if not isinstance(members, list):
                entry(key, prelim or 'other', size, note=f'zip not listed: {members.get("error", "")}')
                continue
            readable = os.path.exists(local) and size <= args.max_mb * 1e6
            tables_by_member = {}
            if readable:
                try:
                    for tname, df in read_tables(local, extract_root):
                        member = tname.split('/', 1)[1].split('::')[0] if '/' in tname else tname
                        tables_by_member.setdefault(os.path.basename(member), []).append(
                            (f'{key}/{tname.split("/", 1)[-1]}', df))
                except Exception:
                    pass
            for mname, msize in members:
                mext = ext_of(mname)
                full = f'{key}/{mname}'
                mprelim = classify_name(mname, record_text)
                tabs = tables_by_member.pop(os.path.basename(mname), None)
                if tabs:
                    handle_tables(full, tabs, msize, mprelim)
                elif mext in ('.nc', '.nc4') and readable:
                    _nc(rec, full, os.path.join(extract_root, re.sub(r'[^\w.-]', '_', key), mname), msize,
                        mprelim, entry, args, zip_src=(local, mname))
                else:
                    note = '' if mext not in TABLE_EXT else (
                        'zip larger than the size limit: members listed, not read' if not readable else
                        'no readable table')
                    entry(full, mprelim or ('other' if mext in TABLE_EXT else 'documentation'), msize, note=note)
            continue
        if ext in TABLE_EXT and not (prelim == 'documentation' and ext not in ('.xlsx', '.xls')):
            if not os.path.exists(local) or size > args.max_mb * 1e6:
                entry(key, prelim or 'other', size,
                      note=f'larger than {args.max_mb} MB: not downloaded' if size > args.max_mb * 1e6
                      else 'not downloaded')
                continue
            try:
                tables = read_tables(local, extract_root)
            except Exception:
                tables = []
            handle_tables(key, tables, size, prelim)
            continue
        if ext in ('.nc', '.nc4'):
            _nc(rec, key, local, size, prelim, entry, args, f=f, session=session, token=token)
            continue
        entry(key, prelim or 'other', size)

    # tables of a relational biological dataset (cruise, station, sample) without other content
    has_bio = any(e['type'] == 'biology' for e, _, _ in pending)
    rs = keyword_scores(record_text)
    rec_type = max(rs, key=rs.get) if max(rs.values()) else None
    for e, df, mapping in pending:
        col_text = ' '.join(name_text(c) for c in df.columns)
        if has_bio and e['type'] == 'other' or (has_bio and e['type'] == 'biology_other'
                                                and not keyword_scores(col_text)['biology_other']):
            e['type'], e['note'] = 'biology', 'station / sample table of the biological dataset'
        elif e['type'] == 'other' and rec_type:
            e['type'] = rec_type                         # names and columns say nothing: record title
        if e['type'] == 'biology':
            e['processed'] = e['saved'] = has_dwca
            e['note'] = e['note'] or (f'in the Darwin Core Archive ({d.get("status")})' if has_dwca else
                                      'biological table not converted: ' + '; '.join(d.get('issues', [])[:1]))
            e['dwca'] = d.get('dwca_file', '')
        elif e['type'] in EXPORT_TYPES:
            try:
                r = export_table(rec, e['file'], e['table'], df, mapping, e['type'], args.other_out)
            except Exception as ex:
                r, e['note'] = None, f'export failed: {type(ex).__name__}: {str(ex)[:120]}'
            if r:
                e.update(r, processed=True, saved=True)
    pending.clear()
    gc.collect()
    return out


def _nc(rec, name, path, size, prelim, entry, args, f=None, session=None, token='', zip_src=None):
    """NetCDF: converted to JSON, then to CSV / XLSX."""
    e = entry(name, prelim or 'spatial', size)
    if size > args.max_nc_mb * 1e6:
        e['note'] = f'larger than {args.max_nc_mb:g} MB: listed only'
        return
    if args.max_nc_per_record and rec['nc'] >= args.max_nc_per_record:
        e['note'] = f'listed only (more than {args.max_nc_per_record} NetCDF files in the record)'
        return
    rec['nc'] += 1
    try:
        if not os.path.exists(path) and zip_src:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with zipfile.ZipFile(zip_src[0]) as z, z.open(zip_src[1]) as src, open(path, 'wb') as dst:
                shutil.copyfileobj(src, dst)
        if not os.path.exists(path) and f and f.get('url'):
            from dwc_converter import download
            download(session, f['url'], path, token)
        import xarray as xr
        with xr.open_dataset(path, decode_times=False) as ds:
            names = ' '.join(name_text(v) + name_text(ds[v].attrs.get('long_name', '')) for v in ds.variables)
        s = keyword_scores(names)
        best = max(s, key=s.get)
        dtype = best if s[best] and best in EXPORT_TYPES else e['type']
        r = export_netcdf(rec, name, path, dtype, args.other_out)
        e.update(r, type=dtype, processed=True, saved=True)
    except Exception as ex:
        e['note'] = f'NetCDF not read: {type(ex).__name__}: {str(ex)[:100]}'


def summarise(files, n_records):
    summary = {}
    for t, (label, colour) in TYPES.items():
        fs = [f for f in files if f['type'] == t]
        if not fs:
            continue
        summary[t] = {'label': label, 'colour': colour,
                      'records': len({f['record_id'] for f in fs}),
                      'found': len({(f['record_id'], f['file']) for f in fs}),
                      'tables': sum(1 for f in fs if f.get('table')),
                      'processed': sum(1 for f in fs if f['processed']),
                      'saved': (len({(f['record_id'], f.get('dwca')) for f in fs if f['saved']}) if t == 'biology'
                                else sum(bool(f.get('csv')) + bool(f.get('xlsx')) + bool(f.get('json')) +
                                         len(f.get('csv_extra') or []) for f in fs if f['saved'])),
                      'size_mb': round(sum(f['size_mb'] for f in fs), 1)}
    return {'records': n_records, 'files': len({(f['record_id'], f['file']) for f in files}),
            'tables_exported': sum(1 for f in files if f['saved'] and f['type'] != 'biology'),
            'types': summary}


def main(argv=None):
    ap = argparse.ArgumentParser(description='Data types of the Zenodo records; export of non-biological tables')
    ap.add_argument('--cache', default=os.path.join(HERE, '..', 'zenodo_cache'))
    ap.add_argument('--out', default=os.path.join(HERE, 'output'))
    ap.add_argument('--ids', nargs='*')
    ap.add_argument('--max-mb', type=float, default=MAX_FILE_MB, help='largest table file read (MB)')
    ap.add_argument('--max-nc-mb', type=float, default=MAX_NC_MB, help='largest NetCDF file converted (MB)')
    ap.add_argument('--max-nc-per-record', type=int, default=0,
                    help='NetCDF files converted per record (0 = all)')
    args = ap.parse_args(argv)
    args.other_out = os.path.join(args.out, 'other_data')
    with open(os.path.join(args.out, 'conversion_results.json'), encoding='utf-8') as fh:
        conv = json.load(fh)['datasets']
    path = os.path.join(args.out, 'data_types.json')
    previous = []
    if args.ids:
        conv = [d for d in conv if str(d['record_id']) in set(args.ids)]
        if os.path.exists(path):
            with open(path, encoding='utf-8') as fh:
                previous = [f for f in json.load(fh)['files'] if f['record_id'] not in set(args.ids)]
    session = robust_session()
    token = os.environ.get('ZENODO_TOKEN', '')
    files = list(previous)
    for k, d in enumerate(conv, 1):
        t0 = time.time()
        try:
            fs = process_record(d, args, session, token)
        except Exception as e:
            fs = [{'record_id': str(d['record_id']), 'title': d.get('title', ''), 'file': '', 'type': 'other',
                   'ext': '', 'size_mb': 0, 'processed': False, 'saved': False,
                   'note': f'unexpected error: {type(e).__name__}: {str(e)[:150]}'}]
        files += fs
        types = {}
        for f in fs:
            types[f['type']] = types.get(f['type'], 0) + 1
        print(f'[{k}/{len(conv)}] {d["record_id"]} files={len(fs)} saved={sum(f["saved"] for f in fs)} '
              f'({time.time() - t0:.0f}s) {types}', flush=True)
    n_records = len({f['record_id'] for f in files})
    payload = {'generated': date.today().isoformat(), 'tool': TOOL, 'summary': summarise(files, n_records),
               'files': files}
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=1, ensure_ascii=False)
    s = payload['summary']
    rows = [{'data type': v['label'], 'records': v['records'], 'files found': v['found'], 'tables': v['tables'],
             'processed': v['processed'], 'saved files': v['saved'], 'size MB': v['size_mb']}
            for v in s['types'].values()]
    with pd.ExcelWriter(os.path.join(args.out, 'Data_Types_Report.xlsx'), engine='xlsxwriter') as xw:
        pd.DataFrame(rows).to_excel(xw, sheet_name='Summary', index=False)
        pd.DataFrame(files).drop(columns=['parameters', 'added'], errors='ignore').to_excel(
            xw, sheet_name='Files', index=False)
    print(json.dumps({k: (v['found'], v['processed'], v['saved']) for k, v in s['types'].items()}, indent=1))


if __name__ == '__main__':
    main()
