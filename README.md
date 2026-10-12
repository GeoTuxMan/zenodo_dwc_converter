# Zenodo → Darwin Core Converter — Task 3

Converts the marine biology datasets discovered on Zenodo into **Darwin Core Archives (DwC-A)**
that can be published through the EMODnet Biology / EurOBIS IPT.

This is **Task 3** of the *Zenodo Inventory* project:

| Task | Repository | What it does |
|---|---|---|
| 1. Inventory | [`zenodo_inventory`](https://github.com/GeoTuxMan/zenodo_inventory) | Searches Zenodo, filters for marine relevance, writes Excel + JSON |
| 2. Quality Check | [`zenodo_qc`](https://github.com/GeoTuxMan/zenodo_qc) | Checks the datasets against the EMODnet Biology data quality criteria |
| 3. Darwin Core conversion | **this repo** | Converts the datasets to Darwin Core Archives (Event / Occurrence / eMoF) |

Results dashboard: <https://sdn2.rmri.ro/zenodo_inventory_new/>

---

## 1. Why Darwin Core?

EMODnet Biology (via EurOBIS / OBIS) only ingests data in the **Darwin Core** standard, packaged
as a Darwin Core Archive: tab-delimited data files + `meta.xml` (describes the files and the
terms) + `eml.xml` (dataset metadata). Datasets on Zenodo come in whatever format their authors
chose — CSV, Excel workbooks with several sheets, species × station matrices, zipped folders —
so every dataset must be restructured before it can be mobilised. This tool automates the
repetitive part of that work and documents what still needs a data manager.

The target structure follows the OBIS-ENV-DATA model, as in the EurOBIS dataset
*Romanian Black Sea Phytoplankton data 2006–2015* (NIMRD):

```
event.txt  (core)           one row per sampling event: date, place, depth, protocol
 ├── occurrence.txt          one row per taxon per event: name, WoRMS LSID, presence/absence
 └── extendedmeasurementorfact.txt
                             abundance, biomass, counts, cover... (+ abiotic event measurements)
meta.xml  eml.xml
```

When a dataset has no sampling-event information (no date, place or station), an
**Occurrence core** archive is written instead.

## 2. Supported inputs

| Input | Handling |
|---|---|
| `.csv` `.tsv` `.txt` `.tab` | encoding sniffing (UTF-8, CP1252, Latin-1), delimiter detection (`,` `;` tab `\|` whitespace), title rows above the header skipped |
| `.xlsx` `.xlsm` `.xls` `.ods` | every sheet read; README / metadata / legend sheets skipped |
| `.zip` | extracted safely (no path traversal, ≤ 200 MB uncompressed), tables read recursively |
| existing **DwC-A** (zip with `meta.xml`) | validated and normalised: read through `meta.xml`, taxonomy enriched with WoRMS, re-packaged |

Files larger than 30 MB are skipped (and reported) unless `--large-files` is given;
README/licence files are ignored.

### Large files (`--large-files`)
Large Zenodo deposits are usually zips of images, sound, rasters or model output that hold
a few small tables. Instead of downloading them whole, the converter reads them **remotely
through HTTP Range requests** (Zenodo answers `206 Partial Content`):

| File > 30 MB | Handling |
|---|---|
| `.zip` | only the zip directory is fetched to list the members; only tabular members (`.csv`, `.txt`, `.xlsx`…, plus `meta.xml`/`eml.xml` of a DwC-A) are extracted, each ≤ 500 MB, text members cut to their first 50 MB; deep paths are flattened (Windows MAX_PATH) |
| `.csv` `.tsv` `.txt` `.tab` | only the first 50 MB are downloaded, cut at the last complete line |
| `.xlsx` `.xls` `.ods` | downloaded whole up to 150 MB |

Example: a 21 GB zip of satellite rasters is inspected with ~1 MB of traffic. What was read
is reported in `issues` (`… zip read remotely - 6 of 171 files are tables (3.1 MB transferred)`).

### Layouts

| Layout | Detection | Conversion |
|---|---|---|
| **long** | a taxon column is mapped | one occurrence per row; measurement columns → eMoF |
| **wide – species as columns** | no taxon column, ≥ 3 numeric columns whose headers look like Latin names | melted; each non-empty cell → one occurrence |
| **wide – samples as columns** | taxon column + ≥ 3 numeric unmapped columns, no measurement columns, unique taxa | melted; column header → `eventID` |

In matrices an empty cell means "not recorded" and is dropped; `0` becomes
`occurrenceStatus = absent`.

## 3. Column mapping

Headers are normalised (lower case, accents removed, units in brackets removed, punctuation →
space) and looked up in a synonym dictionary (`SYNONYMS` in `dwc_converter.py`, English plus
common Romanian, Spanish, French, German and Italian words). Excerpt:

| DwC term | Example source headers |
|---|---|
| `eventID` | sample, sample id, event id, haul, proba |
| `eventDate` | date, sampling date, datetime, data, datum, fecha |
| `year` / `month` / `day` | year, an, anno / month, luna / day, zi |
| `decimalLatitude` / `decimalLongitude` | lat, latitude (°N), lat_dd, y / lon, long, longitude E, x |
| `minimumDepthInMeters` / `maximumDepthInMeters` | depth from, min depth / depth to, max depth; a single *depth* column fills both, `0-10` → 0 / 10 |
| `locationID` / `locality` | station, stn, statie / site, location, locality |
| `scientificName` | species, taxon, scientific name, valid name, specie, organism |
| `scientificNameID` | AphiaID, WoRMS id, LSID |
| `individualCount` | count, number of individuals, N, nr |
| `occurrenceStatus` | presence, presence/absence, PA (1/0, yes/no, +/−) |
| `samplingProtocol` | method, gear, sampling method, instrument |
| `recordedBy` / `identifiedBy` | observer, collector / identified by, taxonomist |
| `lifeStage`, `sex`, `habitat`, `kingdom` … `genus`, `occurrenceRemarks` | direct equivalents |

**Value normalisation**

* dates → ISO 8601 (`24.03.2010`, `3/24/2010`, `12 Mar 2011`, Excel serials, `2006-2015` →
  `2006/2015`; year/month/day columns are combined);
* coordinates → decimal degrees WGS84 (`44,25`, `44°08.3'N`, `28 38 54 W`), out-of-range values
  rejected; `geodeticDatum = EPSG:4326`;
* depths → positive metres;
* taxon labels cleaned (`Chaetoceros sp.` → `Chaetoceros` + `taxonRank = genus`;
  `cf.`/`aff.` → `identificationQualifier`; the raw label kept in `verbatimIdentification`);
  totals, codes and blanks dropped.

### Measurements (eMoF)

Columns containing *abundance, density, biomass, wet/dry weight, cover, CPUE, frequency,
length…* become eMoF records. The unit is read from the header brackets and mapped to the
NERC vocabularies where the concept is unambiguous, e.g.

| Header | measurementType | measurementTypeID | measurementUnit (ID) |
|---|---|---|---|
| `Abundance (cells/L)` | Abundance of biological entity specified elsewhere per unit volume of the water body | P01 `SDBIOL01` | Number per litre (P06 `UCPL`) |
| `Density (ind/m3)` | idem | P01 `SDBIOL01` | Number per cubic metre (P06 `UPMM`) |
| `Density (ind/m2)` | Abundance … per unit area of the bed | P01 `SDBIOL02` | Number per square metre (P06 `UPMS`) |
| `Biomass (mg/m3)` | Wet/Dry weight biomass per unit volume of the water body | *(to assign)* | Milligrams per cubic metre (P06 `UMMC`) |
| integer matrix of counts | Count of biological entity specified elsewhere | P01 `OCOUNT01` | Dimensionless (P06 `UUUU`) |

Where no safe vocabulary match exists, `measurementTypeID` is left empty for the data manager.
Abiotic columns (temperature, salinity, oxygen, chlorophyll-a, seafloor depth) are attached to
the event as event-level eMoF records.

### Identifiers
* `eventID` = `<zenodo id>_<source event id>`, or `<zenodo id>_ev_<hash>` of date + place + depth;
* `occurrenceID` = `<zenodo id>_<source occurrence id>`, or `<zenodo id>_occ_<hash>` of event +
  taxon + life stage + sex (stable between runs);
* `basisOfRecord` = `HumanObservation` (`MaterialSample` when the dataset mentions eDNA /
  metabarcoding).

### Robustness rules
* columns mapped to text terms (taxonomy, locality, basisOfRecord…) are dropped when they hold numbers;
* `basisOfRecord` values outside the DwC vocabulary are replaced by the default;
* start-of-haul positions (`StartLatitude`…) are used when the main coordinates are empty;
* taxon authorship is stripped (`Actinoptychus senarius (Ehr.) Ehr.` → `Actinoptychus senarius`);
* in species-as-columns matrices, headers are accepted as taxa only if WoRMS recognises them;
* a single Latin name in the Zenodo title is used for single-species tables without a taxon column;
* identical rows coming from repeated sheets/files are removed, conflicting occurrenceIDs made unique;
* taxon **reference lists** (code lists, checklists, taxonomies: one row per taxon, higher
  classification columns, no coordinates/stations and fewer than 5 numeric columns) are not
  turned into occurrences - e.g. the FAO ASFIS species code list shipped with model scripts.
  Taxon × sample matrices (ASV tables, abundance matrices) are kept.

### Station tables
If a dataset has a separate station table (station/locality + coordinates, no taxa), events
without coordinates receive them through the station code.

## 4. WoRMS matching

All distinct names are sent to the WoRMS REST service `AphiaRecordsByMatchNames` (≤ 50 names per
call, results cached in `<cache>/_worms_cache.json`). Accepted match types: exact, exact genus /
subgenus, phonetic, near_1, near_2. A match fills `scientificNameID`
(`urn:lsid:marinespecies.org:taxname:<AphiaID>`), `taxonRank`, `kingdom` … `genus` and
`scientificNameAuthorship` (source values are kept when present). Unmatched names are listed in
the issues.

## 5. Status of a dataset

| Status | Meaning |
|---|---|
| **Converted** | archive written; core + all required terms filled for ≥ 90 % of rows (`eventDate`, `decimalLatitude`, `decimalLongitude`, `occurrenceID`, `scientificName`, `occurrenceStatus`, `basisOfRecord`) and `scientificNameID` for ≥ 80 % |
| **Partial** | archive written but some required terms are missing — listed in `missing_required` (typically coordinates or dates kept in a paper / README) |
| **Failed** | files unreadable, or no taxon column / species matrix found in any table |
| **Not attempted** | embargoed / restricted, no tabular files, or all tabular files > 30 MB without `--large-files` |

A *Partial* archive is a working draft: add the missing information (or a mapping override)
and re-run.

## 6. Installation

```bash
git clone https://github.com/GeoTuxMan/zenodo_dwc_converter.git
cd zenodo_dwc_converter
pip install -r requirements.txt
```

Python ≥ 3.9. A Zenodo token is optional (`ZENODO_TOKEN` environment variable, never commit it).

## 7. Usage

```bash
# first 50 open "Included" datasets with tabular files (Black Sea first)
python dwc_converter.py --inventory ../zenodo_inventory/data.json --limit 50

# specific Zenodo records
python dwc_converter.py --ids 15119508 4665967

# only datasets that passed Task 2 QC (Pass / Warnings), best data score first
python dwc_converter.py --qc-results ../zenodo_qc/output/qc_results.json

# also read files > 30 MB remotely, and add the results to the existing report
python dwc_converter.py --ids 15396200 5004746 --large-files --merge

# all options
python dwc_converter.py --inventory data.json --cache ./zenodo_cache --out ./output \
                        --mapping-dir ./mappings --limit 100 --include-review \
                        --large-files --merge --retry-passes 2 --retry-wait 600
```

`--merge` updates `conversion_results.json` in `--out` instead of replacing it - use it when
re-running a subset with `--ids`.

### Mapping overrides
When the automatic mapping is wrong, put `<zenodo id>.json` (or `.yaml`) in the mapping folder.
Keys are table names as reported in the results (`file.xlsx::Sheet1`, or `*` for all tables),
values map source columns to DwC terms (`_depth` = single depth column):

```json
{
  "*": {
    "Stn": "locationID",
    "Sampling day": "eventDate",
    "Taxon_valid": "scientificName",
    "Prof": "_depth"
  }
}
```

### Tests
```bash
python -m pytest -q tests
```

## 8. Outputs (`--out`, default `./output`)

| File | Content |
|---|---|
| `dwca/<zenodo id>_dwca.zip` | one Darwin Core Archive per converted dataset (`meta.xml`, `eml.xml` built from the Zenodo metadata: title, creators, abstract, keywords, licence, DOI, temporal and geographic coverage) |
| `conversion_results.json` | summary + one entry per dataset (status, layout, core, counts, WoRMS match %, field mappings, unmapped columns, missing terms, issues) — consumed by the dashboard |
| `conversion_results.csv` | the same, flattened |
| `Conversion_Report_DwC.xlsx` | sheets *Summary*, *Datasets*, *Field mappings*, *Issues* |

Every archive is re-opened after writing: `meta.xml` must parse, all referenced files must
exist, column counts must match `meta.xml` and core ids must be unique.

## 9. Latest run (12 October 2026)

Full run over all 321 open *Included* records of the Task 1 inventory (Black Sea first, then
datasets whose title/notes mention occurrence, abundance, survey or monitoring); the 22 records
whose tables were all above 30 MB (106 GB in total) were then processed with `--large-files`.

| Status | Datasets |
|---|---|
| Converted | 22 |
| Partial | 145 |
| Failed | 154 |
| Not attempted | 0 |

167 archives written (106 Event core, 61 Occurrence core; 117 long, 50 wide layouts) —
557,713 events, 4,606,540 occurrences, 2,145,894 eMoF records; mean WoRMS match 76 %. All
archives pass the structural check. Most *Partial* datasets lack coordinates and/or dates in
the files themselves (typical of experiment / supplementary-material uploads). Of the *Failed*
ones, 120 hold no taxon-level data (traits, genetics, effort tables, model output, acoustic or
satellite products) and 34 have no readable table (images, audio, rasters, scripts).

Large-file pass: 22 records, 106 GB on Zenodo, about 0.4 GB actually transferred - 1 Converted
(benthic amphipods SDM, 15396200: 136,405 occurrences), 2 Partial (round goby impact 5779207,
seal predation 5004746), 19 Failed (image libraries, acoustic recordings, satellite rasters).

The archives were then checked by Task 2 stage 2 (`zenodo_qc.py --dwca … --run-biocheck`).

Converted examples: Romanian Black Sea Meiobenthos BRIDGE-BS 2022–2024 (15119508), Meiobenthos
GeoEcoMar DOORS 2023 (10863810), Seine nursery biodiversity 1995–2019 (3824354), gelatinous
macrozooplankton North Sea (13903034), Pacific Arctic benthos (12553904, wide matrix).

Zenodo was intermittently unavailable (HTTP 504) during the run; the tool opens a circuit
breaker after three server errors and retries the skipped datasets in extra passes
(`--retry-passes`, `--retry-wait`).

## 10. Known limitations

* Heuristic mapping: unusual headers or multi-row headers may need an override.
* Information held only in papers or README files (coordinates, dates, units) is not extracted.
* Units are taken from header brackets only; unit-less matrices get a generic measurement type.
* Joins between tables are limited to station → coordinates look-ups.
* NetCDF, shapefiles, R data and sequence files are out of scope.

Archives should be reviewed (e.g. with Task 2 QC / EMODnetBiocheck) before publication on the IPT.

## 11. Credits
Developed at NIMRD "Grigore Antipa", Constanța, Romania, for EMODnet Biology. Each converted
dataset keeps the licence and authorship of its Zenodo record.
