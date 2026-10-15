# Zenodo → Darwin Core Converter — Task 3

Converts the marine biology datasets discovered on Zenodo into **Darwin Core Archives (DwC-A)**
.

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

### Dataset-level terms
Filled in the core table from the Zenodo record (values already present in the data are kept),
because EMODnetBiocheck treats `datasetName` and `institutionCode` as required:

| Term | Source |
|---|---|
| `datasetName` | Zenodo title |
| `datasetID` | DOI (`https://doi.org/…`) |
| `institutionCode` | value of an `institutionCode` column in the data when present; else EDMO code when a creator (name or affiliation) belongs to an institution listed in `EDMO_CODES` — NIMRD "Grigore Antipa" = **697**, GeoEcoMar = **850**, CNR-ISMAR Venice = **227**; otherwise the affiliation of the first creator that has one - **to be checked by the data manager**; left empty when Zenodo gives no affiliation |
| `institutionID` | `https://edmo.seadatanet.org/report/<EDMO code>` when the EDMO code is known |
| `license` | Creative Commons licence URL |

`institutions.json` (next to the script) fixes `institutionCode` / `institutionID` for single
records when neither Zenodo nor the data name the institution; it overrides both:

```json
{ "18697857": {"institutionCode": "227", "institutionID": "https://edmo.seadatanet.org/report/227"},
  "18688776": {"institutionCode": "CAFF"},
  "15610736": {"institutionCode": "OBIS"} }
```

### Robustness rules
* columns mapped to text terms (taxonomy, locality, basisOfRecord…) are dropped when they hold numbers;
* `basisOfRecord` values outside the DwC vocabulary are replaced by the default;
* start-of-haul positions (`StartLatitude`…) are used when the main coordinates are empty;
* taxon authorship is stripped (`Actinoptychus senarius (Ehr.) Ehr.` → `Actinoptychus senarius`);
* in species-as-columns matrices, headers are accepted as taxa only if WoRMS recognises them;
* a single species in the Zenodo title - Latin name, or a vernacular name resolved by WoRMS
  (`harbour porpoise` → *Phocoena phocoena*) - is used for tables without a taxon column, but only
  for tables that record the organism: a tag / individual column (`Ptt`, `Bird_ID`), counts or
  abundance, or a name saying sighting, record, catch, sample… Noise, vessel, environment, model
  and raw acoustic-logger tables (CPOD click counts) never become occurrences;
* identical rows coming from repeated sheets/files are removed, conflicting occurrenceIDs made unique;
* taxon **reference lists** (code lists, checklists, taxonomies: one row per taxon, higher
  classification columns, no coordinates/stations and fewer than 5 numeric columns) are not
  turned into occurrences - e.g. the FAO ASFIS species code list shipped with model scripts.
  Taxon × sample matrices (ASV tables, abundance matrices) are kept;
* rows left **without a date** take the date of the same sample - same station code, coordinates
  and depth layer in the same table - when that sample has one date (or one date on ≥ 90 % of
  its rows), otherwise the date of the same station when the station has a single date. The
  number of rows filled this way is reported in `issues`. Example: 15119469 (GeoEcoMar
  zooplankton, July 2024 cruise), where the date was typed only on some rows of each sample.

* tables whose taxon column holds no taxa - metadata sheets (persons, glossary, project list),
  data dictionaries ("Column name / Description / Unit") - are skipped when fewer than 25 % of
  their distinct names (at least 5) match WoRMS;
* R-style taxon names coming from column headers are restored
  (`Cladophora.aff.coelothrix` → `Cladophora coelothrix` + `identificationQualifier = aff.`,
  `Balanus.sp` → `Balanus`);
* headers such as `Sampling_station_country_name`, `Project_name` or `ShipName` are never read as
  `scientificName`;
* tables with a Kingdom ... Genus, Species hierarchy where `Species` holds only the epithet
  (`pertusum`, `cf. oculata?`) get the name Genus + epithet, else the lowest filled rank
  (`Umbellula`, `Hexactinellida`);
* a `time` column holding date and time (`2021-07-06 14:31:39 UTC`) gives the `eventDate` when no
  date column exists; a time without a year is never read as a date; `20150622` is read as
  2015-06-22;
* text files written by R `write.table` (header one field shorter than the rows, because of the
  row names) are realigned, so columns are not shifted by one; quoted space-separated files are
  split on the quotes, not on every space;
* species-as-columns headers may carry a unit (`Acari [#/10 cm**2]`, PANGAEA style);
* in a species-as-columns table each row is one sample: a repeated event code (campaign `T0`, `T1`)
  is combined with the sample column so samples are not merged.

### Derived dates and coordinates
Values missing from the files are taken from outside the data only in the cases below, and are
always labelled:
* `eventDate` - the *collected* period of the Zenodo record (`2024-08-06/2024-08-23`) for events
  without a date; `eventRemarks` says *eventDate derived*;
* coordinates - the centroid of the `locality` / `waterBody` named in the data, from the
  [Marine Regions](https://www.marineregions.org) gazetteer, only when the place is small enough
  (radius ≤ 100 km from its bounding box, or a fjord, lagoon, bay, island… without one) and the name
  is not ambiguous; `coordinateUncertaintyInMeters` holds the radius and `georeferenceRemarks` the
  place and its MRGID. Sea-basin names (*North Atlantic*, *Kattegat*) are never used;
* values stated in the **paper linked to the record** (`paper_sources.py`, below): station
  coordinates matched by station code, the sampling period, and the centre of the study area.

### Dates and coordinates from the linked papers (`paper_sources.py`)
Most uploads are the data behind a paper ("Data from: …"); the sampling period and the station
positions are written in its Methods. `paper_sources.py`, run before the converter,

1. takes the DOIs of the record's related identifiers and description, keeps scholarly works only
   (OpenAlex type *article*, *review*, *preprint*…; data-repository DOIs are skipped);
2. reads an open-access full text: Europe PMC XML, otherwise the open-access PDF found by OpenAlex
   (directly or through the `citation_pdf_url` of a repository landing page - HAL, DSpace…);
   closed papers are reported and skipped;
3. in the Methods / Study area text extracts
   * **station tables** (station + latitude + longitude columns, or station rows of a PDF table),
   * **every coordinate pair** (DMS, decimal minutes, decimal degrees with hemisphere) → centre of
     the study area; the radius is the half-diagonal of their bounding box plus the precision of the
     coordinates as written (whole minutes ≈ 1.3 km, whole degrees ≈ 79 km),
   * the **sampling period** from the sentences about sampling / collection / surveys (author-year
     citations, measurements such as *2000 µm*, page headers and years from the publication year on
     are ignored);
4. writes `output/paper_values.json` and `output/Paper_Values_Review.xlsx` (one row per record:
   paper, source, values, evidence sentences, apply flags).

A value is applied (`date_apply` / `area_apply` = true) only when the period spans ≤ 10 years, the
study area has a radius ≤ 100 km and its centre lies in a region declared for the record in the
inventory; set the flag to `false` in `paper_values.json` to reject a value, then re-run the
converter. In the archive: `eventRemarks` = *eventDate derived: sampling period stated in doi:…*,
`georeferenceRemarks` = *coordinates derived: centre of the study area given in doi:…* (with
`coordinateUncertaintyInMeters`) or *coordinates of station … from the station table of doi:…*.

```bash
python paper_sources.py                    # all records with an archive (cached in zenodo_cache/_papers)
python dwc_converter.py --ids … --merge    # reads output/paper_values.json
```

### Relational datasets
Datasets exported from a database keep dates, coordinates and taxon names in lookup tables
(records → fieldwork → stations, records → taxa). A table with a key column `<x>_ID` / `<x> code`
holding unique values is a lookup table; every table where that column repeats receives the
lookup's other columns through a join, in up to three rounds. Example: MEDUSA (20713997) -
jellyfish records get the date from `Fieldwork_dataset`, coordinates from
`Sampling_stations_dataset` and the name from `WORMS_dataset`.

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

## 7b. Other data types (`data_types.py`)
The same records also hold chemistry, physical oceanography, bathymetry, sediment, tracking,
acoustic, genetic and other data. `data_types.py` lists **every file** of the records processed by
the converter (members of zip archives included - read remotely when the zip is not cached),
classifies it into a data type and saves every non-biological table as **CSV (UTF-8) and XLSX**:

| Data type | Recognised by |
|---|---|
| Biological occurrences / abundance | taxon column or species matrix (→ Darwin Core Archive) |
| Biological measurements without taxon column | biomass, abundance, density, CPUE … of groups (no species names) |
| Traits / morphometrics / diet | length, weight, age, diet, stomach, reproduction … |
| Chemistry / biogeochemistry | nutrients, oxygen, carbonate system, pigments, metals, contaminants, isotopes |
| Physical oceanography / meteorology | CTD, temperature, salinity, currents, waves, light, turbidity, `.cnv` |
| Bathymetry / topography | bathymetry, multibeam, DEM, `.xyz` / `.asc` / bathymetric rasters |
| Sediment / geology | grain size, mud / sand / gravel, organic matter, cores |
| Tracking / telemetry | GPS, Argos, PTT, tags, dives, `.gpx` / `.kml` |
| Acoustic | hydrophones, C-POD / F-POD, clicks, noise / SPL, `.wav` |
| Genetic / eDNA / sequences | `.fastq`, `.fasta`, `.vcf`, trees, OTU / ASV tables |
| Imagery / video, GIS / rasters / model output, R / MATLAB binaries, documentation / code | file extension |

The type of a table comes from its file / sheet name and column headers (the record title breaks
ties); one-row lookup sheets and data dictionaries count as documentation, and the station /
sample tables of a relational biological dataset count as biology.

Each exported table keeps all original columns (empty rows / columns dropped), with
`date_ISO8601`, `latitude_decimal`, `longitude_decimal` added when date or coordinate columns are
recognised and not already in that form. The XLSX has three sheets: *Data*, *Columns* (original
header, unit from the header, parameter group, Darwin Core term, numeric / text, fill rate,
example) and *Source* (Zenodo record, DOI, creators, license, original file, citation). Records
under a **NoDerivatives** license (CC BY-ND) get a change of format only - no columns are added.
Files over 30 MB are listed, not downloaded.

**NetCDF** files are converted to **JSON** first (CF-JSON-like: source, global attributes,
dimensions, and every variable with its dimensions, attributes, shape and data as nested arrays;
times in ISO 8601, missing values as `null`, float32 to 7 significant digits), then **from the
JSON** to CSV / XLSX: one long table per group of variables sharing the same dimensions (columns =
dimension coordinates + variables; XLSX sheets *Data*, *Variables*, *Global attributes*, *Source*).
Gridded fields with more than 100,000 cells stay in the JSON only (as a table they would be tens of
millions of rows); NetCDF files inside zip archives larger than 30 MB are listed, not converted
(`--max-nc-mb`, `--max-nc-per-record` change these limits).

```bash
python data_types.py                 # all records of output/conversion_results.json
python data_types.py --ids 3998049   # some records (merged into the existing results)
```

Outputs: `output/other_data/<type>/<zenodo id>/<table>.csv|.xlsx|.json`, `output/data_types.json`
(per file: type, size, table, rows, columns, parameter groups, processed / saved, download paths;
summary per type: records, files found, tables processed, files saved) and
`output/Data_Types_Report.xlsx`.

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
| `paper_values.json`, `Paper_Values_Review.xlsx` | dates / coordinates found in the linked papers (`paper_sources.py`) |
| `other_data/`, `data_types.json`, `Data_Types_Report.xlsx` | data types of all files and the exported non-biological tables (`data_types.py`) |

Every archive is re-opened after writing: `meta.xml` must parse, all referenced files must
exist, column counts must match `meta.xml` and core ids must be unique.

## 9. Latest run (15 October 2026)

Full run over all 321 open *Included* records of the Task 1 inventory (Black Sea first, then
datasets whose title/notes mention occurrence, abundance, survey or monitoring); the 22 records
whose tables were all above 30 MB (106 GB in total) were then processed with `--large-files`,
and all 321 records were reconverted with `--large-files` after the fixes of 13 October
(see *Relational datasets* and the taxon-name rules above), then again after the hierarchy-name,
observation-table and derived date/coordinate rules (fewer false occurrences from noise, vessel and
logger tables; *Derived dates and coordinates* filled 2 datasets with dates and 5 with coordinates).

**Step 3 - linked papers** (`paper_sources.py`): 68 of the 153 records with an archive link to a
paper; 43 full texts were open access (Europe PMC or open-access PDF), 25 were closed. Dates or
coordinates were found for 38; after the checks (period ≤ 10 years, study area ≤ 100 km, inside
the declared region and not on land, not a compilation of records) 24 periods and 10 study areas
were applied: eventDate of 1,153 events in 14 datasets and coordinates of 618 events in 7
datasets. Four datasets became *Converted* (macrofouling off Elba 5014110, kittiwakes at Hornøya
4971554, mudflat meiofauna 4943866, round goby 5779207). No station table could be matched to the
data (station positions are mostly in figures or supplements).

| Status | Datasets |
|---|---|
| Converted | 36 |
| Partial | 117 |
| Failed | 168 |
| Not attempted | 0 |

153 archives written (112 Event core, 41 Occurrence core; 103 long, 50 wide layouts) —
521,331 events, 5,326,879 occurrences, 2,575,204 eMoF records; mean WoRMS match 88 %
(76 % before the fixes, which removed metadata sheets and data dictionaries read as taxa). All
archives pass the structural check. Most *Partial* datasets lack coordinates and/or dates in
the files themselves (typical of experiment / supplementary-material uploads). Of the *Failed*
ones, 141 hold no taxon-level data (traits, genetics, effort tables, model output, acoustic or
satellite products, species codes without a lookup table) and 27 have no readable table (images, audio, rasters, scripts).

Large-file pass: 22 records, 106 GB on Zenodo, about 0.4 GB actually transferred - 1 Converted
(benthic amphipods SDM, 15396200: 136,405 occurrences), 2 Partial (round goby impact 5779207,
seal predation 5004746), 19 Failed (image libraries, acoustic recordings, satellite rasters).
In the rebuild, three datasets that had skipped some large files gained data: MEDUSA jellyfish
records (20713997), SAMBAH harbour porpoise detections (6639075, now Converted) and BENTH-MARC
zoobenthos (20782806).

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
* Joins between tables need a shared key column (`<x>_ID`, `<x> code`) or station → coordinates look-ups.
* Shapefiles, rasters, R data and sequence files are only listed and classified (`data_types.py`), not converted.

Archives should be reviewed (e.g. with Task 2 QC / EMODnetBiocheck).

## 11. Credits
Developed at NIMRD "Grigore Antipa", Constanța, Romania, for EMODnet Biology. Each converted
dataset keeps the licence and authorship of its Zenodo record.
