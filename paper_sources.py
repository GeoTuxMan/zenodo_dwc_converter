"""
Sampling dates and coordinates stated in the papers linked to the Zenodo records (Task 3, step 3).

Many Zenodo uploads are the data behind a paper ("Data from: ..."): the station list, the
coordinates and the sampling period are in the paper, not in the files. For every record that has
a Darwin Core Archive this script

  1. collects the DOIs of the linked papers (Zenodo related identifiers and the description),
     keeping only scholarly works (OpenAlex type article / review / preprint / book chapter);
  2. gets an open-access full text: Europe PMC XML first, then an open-access PDF found by
     OpenAlex (text read with PyMuPDF); closed papers are reported and skipped;
  3. extracts from the Methods / Study area part of the text
       - station tables (station name + latitude + longitude),
       - every coordinate pair (DMS or decimal with hemisphere) -> study-area centroid and radius,
       - the sampling period (dates in sentences about sampling / collection / surveys);
  4. writes output/paper_values.json (read by dwc_converter.py) and Paper_Values_Review.xlsx.

Values are applied by the converter only to events that have no date / coordinates in the data,
and are labelled as derived from the paper. Set "date_apply" / "area_apply" to false in
paper_values.json to reject a value; low-confidence values are written with false already.

    python paper_sources.py                      # all records with an archive
    python paper_sources.py --ids 4970001 15057989
"""
import argparse
import json
import math
import os
import re
import sys
import time

import pandas as pd
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dwc_converter import (parse_coordinate, robust_session, norm_station,  # noqa: E402
                           MAX_DERIVED_UNCERTAINTY_M)

OPENALEX = 'https://api.openalex.org/works/doi:{}'
EPMC_SEARCH = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
EPMC_XML = 'https://www.ebi.ac.uk/europepmc/webservices/rest/{}/fullTextXML'
XYLOOKUP = 'https://api.obis.org/xylookup'
PAPER_TYPES = {'article', 'review', 'preprint', 'book-chapter', 'letter'}
SKIP_RELATIONS = {'isversionof', 'hasversion', 'ispartof', 'haspart', 'isidenticalto', 'isnewversionof',
                  'ispreviousversionof', 'isderivedfrom', 'issourceof', 'isvariantformof'}
# DOI prefixes of data repositories / services (never a paper)
DATA_DOI_RE = re.compile(r'^10\.(5281|5061|15468|48670|24381|5067|5285|1594|25921|34892|12770|24400|'
                         r'18150|17882|6084|5880|14284|7289|4060)/', re.I)
DOI_RE = re.compile(r'10\.\d{4,9}/[^\s"<>,;)\]]+')
MAX_PAPERS = 3
# compilations: the period / area of the authors' own fieldwork is not that of the records
COMPILATION_RE = re.compile(r'databases?|data ?base|compil|checklist|literature|review|meta-?analys|synthes|'
                            r'harmoni[sz]ed|global|museum|collections?\b|historical records|set of', re.I)
# generous boxes (lat min, lat max, lon min, lon max) of the regions of the Task 1 inventory: a study
# area outside every declared region comes from another place named in the paper (not applied)
REGION_BOXES = {'atlantic': (-60, 75, -100, 30), 'baltic sea': (53, 66.5, 9, 31),
                'mediterranean sea': (30, 46, -6, 37), 'north sea': (50, 62, -5, 13),
                'arctic': (60, 90, -180, 180), 'black sea': (40, 47.5, 27, 42)}
MAX_DATE_SPAN_YEARS = 10

DEG = r'[°º˚◦∘o]'
MIN = r"['′’ʹ`´]"
SEC = r'["″”ʺ]|\'\'|′′'
NUM = r'\d{1,3}(?:[.,]\d+)?'
COORD_RE = re.compile(
    rf'(?:(?P<h1>[NSEW])\s*(?P<d1>\d{{1,3}}(?:[.,]\d+)?)\s*{DEG}\s*(?:(?P<m1>\d{{1,2}}(?:[.,]\d+)?)\s*{MIN}?\s*'
    rf'(?:(?P<s1>\d{{1,2}}(?:[.,]\d+)?)\s*(?:{SEC})?)?)?)'
    rf'|(?:(?P<d2>\d{{1,3}}(?:[.,]\d+)?)\s*{DEG}\s*(?:(?P<m2>\d{{1,2}}(?:[.,]\d+)?)\s*{MIN}?\s*'
    rf'(?:(?P<s2>\d{{1,2}}(?:[.,]\d+)?)\s*(?:{SEC})?\s*)?)?\s*(?P<h2>[NSEW])\b)'
    rf'|(?:(?P<d3>\d{{1,3}}[.,]\d{{2,}})\s*(?P<h3>[NSEW])\b)')

MONTHS = ['january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september',
          'october', 'november', 'december']
MONTH_RE = r'(jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|june?|july?|aug(?:ust)?|sep(?:t(?:ember)?)?|' \
           r'oct(?:ober)?|nov(?:ember)?|dec(?:ember)?)\.?'
YEAR_RE = r'(1[89]\d{2}|20[0-4]\d)'
SAMPLING_RE = re.compile(r'sampl|collect|survey|cruise|expedition|field ?work|field ?trip|campaign|'
                         r'monitor|trawl|caught|captured|deploy|tagged|\bdives?\b|transects?\b|\bgrabs?\b|'
                         r'\bhauls?\b|census|conducted|carried out|visited', re.I)
PDF_META_RE = re.compile(rb'<meta[^>]+name=["\']citation_pdf_url["\'][^>]+content=["\']([^"\']+)', re.I)
PDF_META_RE2 = re.compile(rb'<meta[^>]+content=["\']([^"\']+)["\'][^>]+name=["\']citation_pdf_url', re.I)
# author-year citations: "(Smith et al. 2010; Jones & Li, 2012)", "Olaso (1998)", "Smith et al. 2010"
CITATION_RE = re.compile(r'\(\s*(?:19|20)\d{2}[a-z]?\s*\)|\((?=[^()]*\d{4})(?![^()]*\b(?i:jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]*\.?\s'
                         r'*\d)[^()]*[A-Za-z][^()]*\)|'
                         r'(?!(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec|Spring|Summer|Autumn|Fall|'
                         r'Winter|In|Since|From|During|Until|And|To)[a-z]*\.?\s)'
                         r'[A-Z][^\s,;()]+(?:\s+et\s+al\.?|\s+(?:and|&)\s+[A-Z][^\s,;()]+)?,?\s*\(?(?:19|20)\d{2}[a-z]?\)?'
                         r'(?=[\s,;.)]|$)')
# running heads / footers of PDF pages ("Frontiers in ... | www.frontiersin.org | 2 April 2018")
PAGE_FURNITURE_RE = re.compile(r'www\.|https?:|\bvolume\b|\bvol\.\s*\d|article\s+\d|\bdoi\b|received|accepted|'
                               r'published|copyright|©|licen[cs]e', re.I)
METHODS_HEAD_RE = re.compile(r'(?im)^\s*(?:\d+(?:\.\d+)*\.?\s*)?(materials?\s+and\s+methods|methods?|'
                             r'study\s+(?:area|site)s?|sampling|field\s+work|data\s+collection)\b.{0,60}$')
END_HEAD_RE = re.compile(r'(?im)^\s*(?:\d+(?:\.\d+)*\.?\s*)?(results?|discussion|references|literature\s+cited|'
                         r'acknowledge?ments?)\s*$')


# ---------------------------------------------------------------------------
# paper discovery and full text
# ---------------------------------------------------------------------------
def linked_dois(meta):
    md = meta.get('metadata', {})
    out = []
    for r in md.get('related_identifiers', []) or []:
        if str(r.get('relation', '')).lower() in SKIP_RELATIONS:
            continue
        out += DOI_RE.findall(str(r.get('identifier', '')))
    out += DOI_RE.findall(re.sub(r'<[^<]+?>', ' ', md.get('description', '') or ''))
    own = str(md.get('doi', '')).lower()
    seen = []
    for d in out:
        d = d.rstrip('.').lower()
        if d != own and not DATA_DOI_RE.match(d) and d not in seen:
            seen.append(d)
    return seen


class Papers:
    def __init__(self, cache_dir, session):
        self.dir = os.path.join(cache_dir, '_papers')
        os.makedirs(self.dir, exist_ok=True)
        self.session = session

    def _path(self, doi, ext):
        return os.path.join(self.dir, re.sub(r'[^\w.-]+', '_', doi) + ext)

    def _get_json(self, url, path, **kw):
        if os.path.exists(path):
            with open(path, encoding='utf-8') as fh:
                return json.load(fh)
        try:
            r = self.session.get(url, timeout=60, **kw)
            data = r.json() if r.status_code == 200 else {}
        except (requests.RequestException, ValueError):
            return {}
        with open(path, 'w', encoding='utf-8') as fh:
            json.dump(data, fh)
        time.sleep(0.2)
        return data

    def work(self, doi):
        return self._get_json(OPENALEX.format(doi), self._path(doi, '.openalex.json'))

    def on_land(self, lat, lon):
        """True when OBIS xylookup puts the point inside the coastline and > 5 m above sea level
        (the centre of a coastal study area can fall on land); None when unknown."""
        path = os.path.join(self.dir, f'_xy_{lat}_{lon}.json')
        if os.path.exists(path):
            with open(path, encoding='utf-8') as fh:
                r = json.load(fh)
        else:
            try:
                res = self.session.post(XYLOOKUP, json={'points': [[lon, lat]], 'areas': False, 'grids': True,
                                                        'shoredistance': True}, timeout=60)
                r = res.json()[0] if res.status_code == 200 else None
            except (requests.RequestException, ValueError, IndexError):
                r = None
            if r is None:
                return None
            with open(path, 'w', encoding='utf-8') as fh:
                json.dump(r, fh)
        bathy = (r.get('grids') or {}).get('bathymetry')
        return (r.get('shoredistance') or 0) < 0 and (bathy is None or bathy < -5)

    def fulltext(self, doi, work):
        """(kind, text, tables, url); tables = list of row lists (Europe PMC XML only)."""
        epmc = self._get_json(EPMC_SEARCH, self._path(doi, '.epmc.json'),
                              params={'query': f'DOI:"{doi}"', 'format': 'json', 'resultType': 'lite'})
        pmcid = next((x.get('pmcid') for x in (epmc.get('resultList') or {}).get('result', [])
                      if x.get('pmcid') and x.get('isOpenAccess') == 'Y'), None)
        if pmcid:
            xml = self._fetch(EPMC_XML.format(pmcid), self._path(doi, '.xml'))
            if xml and b'<article' in xml[:5000]:
                text, tables = parse_jats(xml)
                return 'Europe PMC full text', text, tables, f'https://europepmc.org/article/PMC/{pmcid}'
        urls, pages = [], []
        for loc in [work.get('best_oa_location') or {}] + (work.get('locations') or []):
            if loc and loc.get('pdf_url') and loc['pdf_url'] not in urls:
                urls.append(loc['pdf_url'])
            if loc and loc.get('landing_page_url') and loc['landing_page_url'] not in pages:
                pages.append(loc['landing_page_url'])
        tried = 0
        for url in urls + pages:
            if tried >= 6:
                break
            tried += 1
            data = self._fetch(url, self._path(doi, f'.{tried}.bin'))
            if data[:5] != b'%PDF-':                     # landing page: follow its citation_pdf_url
                links = [m.group(1) for rx in (PDF_META_RE, PDF_META_RE2) for m in rx.finditer(data)]
                links = sorted(dict.fromkeys(links), key=lambda u: not re.search(rb'pdf', u, re.I))
                for k, link in enumerate(links[:3]):
                    link = requests.compat.urljoin(url, link.decode('utf-8', 'ignore').replace('&amp;', '&'))
                    data = self._fetch(link, self._path(doi, f'.{tried}.{k}.pdf'))
                    if data[:5] == b'%PDF-':
                        url = link
                        break
            if data[:5] == b'%PDF-':
                text = pdf_text(data)
                if len(text) > 2000:
                    return 'open-access PDF', text, [], url
        return None, '', [], ''

    def _fetch(self, url, path):
        if os.path.exists(path):
            with open(path, 'rb') as fh:
                return fh.read()
        try:
            r = self.session.get(url, timeout=90, headers={'Accept': 'application/pdf,application/xml,*/*'})
            data = r.content if r.status_code == 200 else b''
        except requests.RequestException:
            data = b''
        with open(path, 'wb') as fh:
            fh.write(data)
        time.sleep(0.5)
        return data


def parse_jats(xml):
    from lxml import etree
    root = etree.fromstring(xml, parser=etree.XMLParser(recover=True, huge_tree=True))
    parts = []
    for sec in root.iter('sec'):
        if sec.getparent() is not None and sec.getparent().tag == 'sec':
            continue                                   # top-level sections only (children included)
        title = ''.join(sec.findtext('title') or '')
        parts.append(f'\n{title}\n' + ' '.join(t.strip() for t in sec.itertext()))
    tables = []
    for tw in root.iter('table'):
        rows = []
        for tr in tw.iter('tr'):
            rows.append([' '.join(''.join(c.itertext()).split()) for c in tr if c.tag in ('td', 'th')])
        if rows:
            tables.append(rows)
    return '\n'.join(parts), tables


def pdf_text(data):
    import fitz
    try:
        doc = fitz.open(stream=data, filetype='pdf')
        return '\n'.join(page.get_text() for page in doc)
    except Exception:
        return ''


def methods_part(text):
    """Methods / study area text (up to Results); the whole text without references otherwise."""
    m = METHODS_HEAD_RE.search(text)
    if m:
        end = END_HEAD_RE.search(text, m.end())
        part = text[m.start(): end.start() if end and end.start() - m.start() > 500 else m.start() + 40000]
        if len(part) > 500:
            return part
    end = re.search(r'(?im)^\s*(references|literature cited)\s*$', text)
    return text[: end.start()] if end else text


# ---------------------------------------------------------------------------
# extraction
# ---------------------------------------------------------------------------
def _num(s):
    return float(s.replace(',', '.')) if s else 0.0


def coord_matches(text):
    """[(pos, end, value, axis)] with axis 'lat' / 'lon'."""
    out = []
    for m in COORD_RE.finditer(text):
        g = m.groupdict()
        for i in ('1', '2', '3'):
            if g.get('h' + i):
                d, mi, s, h = _num(g['d' + i]), _num(g.get('m' + i)), _num(g.get('s' + i)), g['h' + i]
                last = g.get('s' + i) or g.get('m' + i) or g['d' + i]
                unit = 3600 if g.get('s' + i) else 60 if g.get('m' + i) else 1
                break
        if mi >= 60 or s >= 60:
            continue
        v = d + mi / 60 + s / 3600
        axis = 'lat' if h in 'NS' else 'lon'
        if (axis == 'lat' and v > 90) or v > 180:
            continue
        decimals = len(re.split(r'[.,]', last)[1]) if re.search(r'[.,]', last) else 0
        precision = 1 / unit / 10 ** decimals / 2 * 111320 * 1.414        # half a step, both axes
        out.append((m.start(), m.end(), -v if h in 'SW' else v, axis, precision))
    return out


def coord_pairs(text, max_gap=60):
    """Latitude/longitude pairs written next to each other: [(lat, lon, start, end)]."""
    ms, pairs, i = coord_matches(text), [], 0
    while i < len(ms) - 1:
        a, b = ms[i], ms[i + 1]
        if a[3] != b[3] and b[0] - a[1] <= max_gap:
            lat, lon = (a[2], b[2]) if a[3] == 'lat' else (b[2], a[2])
            pairs.append((round(lat, 5), round(lon, 5), a[0], b[1], max(a[4], b[4])))
            i += 2
        else:
            i += 1
    return pairs


def station_tables(tables):
    """{normalised station: (label, lat, lon)} from tables with station + latitude + longitude columns."""
    out = {}
    for rows in tables:
        head_i = next((i for i, r in enumerate(rows[:3]) if any(re.search(r'lat', c, re.I) for c in r)), None)
        if head_i is None:
            continue
        head = [c.lower() for c in rows[head_i]]
        lat_c = next((i for i, c in enumerate(head) if 'lat' in c), None)
        lon_c = next((i for i, c in enumerate(head) if re.search(r'lon|lng', c)), None)
        st_c = next((i for i, c in enumerate(head) if re.search(r'station|site|stn|location|locality|area|'
                                                              r'name|code|transect|id\b', c)), 0)
        if lat_c is None or lon_c is None or st_c in (lat_c, lon_c):
            continue
        for r in rows[head_i + 1:]:
            if len(r) <= max(lat_c, lon_c, st_c):
                continue
            lat, lon = parse_coordinate(r[lat_c], True), parse_coordinate(r[lon_c], False)
            if 'w' in r[lon_c].lower() and lon is not None and lon > 0:
                lon = -lon
            if 's' in r[lat_c].lower() and lat is not None and lat > 0:
                lat = -lat
            if lat is not None and lon is not None and r[st_c].strip():
                out[norm_station(r[st_c])] = (r[st_c].strip()[:60], round(lat, 5), round(lon, 5))
    return out


def station_lines(text):
    """Station rows of tables flattened into PDF text: '<label> <lat> <lon>' on one line."""
    out = {}
    for line in text.splitlines():
        for lat, lon, s, e, _p in coord_pairs(line, max_gap=25):
            label = line[:s].strip(' \t:|-–')
            if 1 <= len(label) <= 25 and re.search(r'[A-Za-z0-9]', label) and len(label.split()) <= 3:
                out[norm_station(label)] = (label, lat, lon)
            break
    return out if len(out) >= 2 else {}


def study_area(pairs):
    """Centroid and radius (m) of all coordinate pairs [(lat, lon, precision_m)]; the radius
    includes the precision of the coordinates as written (whole minutes: about 1.3 km)."""
    if not pairs:
        return None
    lats, lons = [p[0] for p in pairs], [p[1] for p in pairs]
    lat0, lon0 = (min(lats) + max(lats)) / 2, (min(lons) + max(lons)) / 2
    dlat = (max(lats) - min(lats)) * 111320 / 2
    dlon = (max(lons) - min(lons)) * 111320 / 2 * max(0.05, math.cos(math.radians(lat0)))
    prec = max((p[2] for p in pairs if len(p) > 2), default=100)
    return {'lat': round(lat0, 5), 'lon': round(lon0, 5),
            'uncertainty_m': int(max(math.hypot(dlat, dlon) + prec, 100)),
            'n_points': len({(p[0], p[1]) for p in pairs})}


def _month(s):
    return next(i + 1 for i, m in enumerate(MONTHS) if m.startswith(s.lower().rstrip('.')[:3]))


def sentence_dates(sentence):
    """[(year, month or None)] mentioned in one sentence (citations removed)."""
    s = CITATION_RE.sub(' ', sentence)
    s = re.sub(r'\b(?:doi|https?)\S+', ' ', s)
    found, used = [], []
    for m in re.finditer(rf'(?i)\b{MONTH_RE}\s*(?:[-–/]|to|and|until)\s*{MONTH_RE}\s*,?\s*{YEAR_RE}\b', s):
        y = int(m.group(3))
        a, b = _month(m.group(1)), _month(m.group(2))
        found += [(y if a <= b else y - 1, a), (y, b)]
        used.append(m.span())
    for m in re.finditer(rf'(?i)\b(?:\d{{1,2}}(?:st|nd|rd|th)?\s+)?{MONTH_RE}\s*(?:\d{{1,2}}(?:st|nd|rd|th)?,?\s*)?'
                         rf'{YEAR_RE}\b', s):
        if not any(a <= m.start() < b for a, b in used):
            found.append((int(m.group(2)), _month(m.group(1))))
            used.append(m.span())
    for m in re.finditer(rf'\b{YEAR_RE}\b(?!\s*(?:[µμu]m|mm|cm|km|m|ml|l|g|kg|mg|ind|individuals|samples|'
                         rf'specimens|stations|sites|reads|sequences|m2|m3|%|°|ha|psu)\b)', s):
        if not any(a <= m.start() < b for a, b in used):
            found.append((int(m.group(1)), None))
    return found


def sampling_period(text, pub_year):
    """(ISO interval, evidence sentences) from the sentences about sampling."""
    text = re.sub(r'-\n(?=[a-z])', '', text)
    sentences = re.split(r'(?<=[.;])\s+(?=[A-Z(])', re.sub(r'\s+', ' ', text))
    dates, evidence = [], []
    for s in sentences:
        if not SAMPLING_RE.search(s) or len(s) > 600 or PAGE_FURNITURE_RE.search(s):
            continue
        # data are collected before the paper: bare years of the publication year onward are dropped
        ds = [d for d in sentence_dates(s) if 1900 <= d[0] <= (pub_year or 2100)
              and (d[1] or d[0] < (pub_year or 2100))]
        if ds:
            dates += ds
            evidence.append(s.strip()[:220])
    if not dates:
        return None, []
    lo, hi = min(dates, key=lambda d: (d[0], d[1] or 1)), max(dates, key=lambda d: (d[0], d[1] or 12))
    fmt = lambda d: f'{d[0]}-{d[1]:02d}' if d[1] else str(d[0])  # noqa: E731
    period = fmt(lo) if lo == hi else f'{fmt(lo)}/{fmt(hi)}'
    return period, evidence[:3]


def period_years(period):
    ys = [int(y) for y in re.findall(r'\d{4}', period or '')]
    return (max(ys) - min(ys)) if ys else None


# ---------------------------------------------------------------------------
def in_regions(lat, lon, regions):
    boxes = [REGION_BOXES[r.strip().lower()] for r in str(regions or '').split(',')
             if r.strip().lower() in REGION_BOXES]
    return not boxes or any(a <= lat <= b and c <= lon <= d for a, b, c, d in boxes)


def record_values(rid, meta, papers, regions=''):
    out = {'record_id': rid, 'papers': [], 'status': 'no linked paper'}
    for doi in linked_dois(meta)[:8]:
        work = papers.work(doi)
        if not work or work.get('type') not in PAPER_TYPES:
            continue
        info = {'doi': doi, 'title': (work.get('title') or '')[:200], 'year': work.get('publication_year'),
                'open_access': bool((work.get('open_access') or {}).get('is_oa'))}
        out['papers'].append(info)
        if len(out['papers']) > MAX_PAPERS:
            break
        kind, text, tables, url = papers.fulltext(doi, work)
        if not kind:
            info['fulltext'] = 'not available (closed or not downloadable)'
            continue
        info['fulltext'], info['url'] = kind, url
        part = methods_part(text)
        stations = station_tables(tables) or station_lines(part)
        pairs = [(p[0], p[1], p[4]) for p in coord_pairs(part)] or [(v[1], v[2], 100) for v in stations.values()]
        area = study_area(pairs)
        period, ev = sampling_period(part, work.get('publication_year'))
        if not (stations or area or period):
            info['fulltext'] += ': no dates or coordinates found'
            continue
        out.update(status='values found', paper_doi=doi, source=kind, url=url)
        compilation = bool(COMPILATION_RE.search(meta.get('metadata', {}).get('title', '') or ''))
        if period:
            span = period_years(period)
            out.update(eventDate=period, date_evidence=' | '.join(ev),
                       date_apply=bool(span is not None and span <= MAX_DATE_SPAN_YEARS) and not compilation)
        if area:
            inside = in_regions(area['lat'], area['lon'], regions)
            land = papers.on_land(area['lat'], area['lon']) if inside and not compilation else None
            out.update(area=area, area_apply=area['uncertainty_m'] <= MAX_DERIVED_UNCERTAINTY_M and inside
                       and not compilation and not land)
            if land:
                out['area_note'] = 'centre of the study area on land: not applied'
            elif compilation:
                out['area_note'] = 'compilation of records: not applied'
            elif not inside:
                out['area_note'] = f'outside the declared region ({regions})'
            elif area['uncertainty_m'] > MAX_DERIVED_UNCERTAINTY_M:
                out['area_note'] = 'study area too large'
        if compilation:
            out['date_note'] = 'compilation of records: not applied'
        if stations:
            out['stations'] = {k: list(v) for k, v in stations.items()}
        return out
    if out['papers']:
        out['status'] = ('no open-access full text' if all('not available' in p.get('fulltext', '')
                                                             for p in out['papers'])
                         else 'full text read, no dates or coordinates found')
    return out


def write_review(values, conv, path):
    rows = []
    titles = {str(d['record_id']): d for d in conv}
    for rid, v in values.items():
        d = titles.get(rid, {})
        a = v.get('area') or {}
        rows.append({'record_id': rid, 'title': d.get('title', ''), 'conversion_status': d.get('status', ''),
                     'missing_required': ', '.join(d.get('missing_required') or []),
                     'status': v['status'], 'paper_doi': v.get('paper_doi', ''), 'source': v.get('source', ''),
                     'url': v.get('url', ''), 'eventDate': v.get('eventDate', ''),
                     'date_apply': v.get('date_apply', ''), 'date_evidence': v.get('date_evidence', ''),
                     'area_lat': a.get('lat', ''), 'area_lon': a.get('lon', ''),
                     'area_uncertainty_m': a.get('uncertainty_m', ''), 'area_points': a.get('n_points', ''),
                     'area_apply': v.get('area_apply', ''), 'area_note': v.get('area_note', ''),
                     'stations': len(v.get('stations') or {}),
                     'papers': '; '.join(f"{p['doi']} ({p.get('fulltext', '')})" for p in v['papers'])})
    pd.DataFrame(rows).to_excel(path, index=False, sheet_name='Paper values')


def main(argv=None):
    ap = argparse.ArgumentParser(description='Sampling dates and coordinates from the papers linked to '
                                             'the Zenodo records')
    ap.add_argument('--cache', default=os.path.join(HERE, '..', 'zenodo_cache'))
    ap.add_argument('--out', default=os.path.join(HERE, 'output'))
    ap.add_argument('--ids', nargs='*')
    args = ap.parse_args(argv)
    with open(os.path.join(args.out, 'conversion_results.json'), encoding='utf-8') as fh:
        conv = json.load(fh)['datasets']
    todo = [str(d['record_id']) for d in conv if d.get('dwca_file')]
    regions = {str(d['record_id']): d.get('regions', '') for d in conv}
    if args.ids:
        todo = [r for r in todo if r in set(args.ids)] or list(args.ids)
    path = os.path.join(args.out, 'paper_values.json')
    values = {}
    if os.path.exists(path):
        with open(path, encoding='utf-8') as fh:
            values = json.load(fh)
    papers = Papers(args.cache, robust_session())
    for k, rid in enumerate(todo, 1):
        try:
            with open(os.path.join(args.cache, rid, '_record.json'), encoding='utf-8') as fh:
                meta = json.load(fh)
        except OSError:
            continue
        old = values.get(rid, {})
        v = record_values(rid, meta, papers, regions.get(rid, ''))
        for key in ('date_apply', 'area_apply'):         # keep the data manager's decisions
            if key in old and old.get('eventDate' if key == 'date_apply' else 'area') == \
                    v.get('eventDate' if key == 'date_apply' else 'area') and key in v:
                v[key] = old[key]
        values[rid] = v
        print(f'[{k}/{len(todo)}] {rid} {v["status"]:<45} {v.get("eventDate", "")}'
              f' {v.get("area", {}).get("uncertainty_m", "")} st={len(v.get("stations") or {})}', flush=True)
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(values, fh, indent=1, ensure_ascii=False)
    write_review(values, conv, os.path.join(args.out, 'Paper_Values_Review.xlsx'))
    st = {}
    for v in values.values():
        st[v['status']] = st.get(v['status'], 0) + 1
    print(json.dumps(st, indent=1))


if __name__ == '__main__':
    main()
