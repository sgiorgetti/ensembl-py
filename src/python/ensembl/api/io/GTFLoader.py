# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Ensembl GTF files parser/loader"""

import contextlib
from collections import namedtuple
from dataclasses import dataclass, field
from pathlib import Path
import re
import warnings
from ensembl.api.core import Feature, Transcript, Exon, Slice

GTF_DELIMITERS = ["\t", ";"]
PARENTS = {"gene": "", "transcript": "gene", "exon": "transcript"}
STRAND_MAP = {"+": 1, "-": -1, ".": 0, "?": 0}

# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
GTFRawRecord = namedtuple('RawRecord', 'seqname, source, feature, start, end, score, strand, frame, attributes')


def _get_attributes(gtf_row: GTFRawRecord) -> dict[str,str]:
    if not gtf_row:
        return None
    return dict(re.findall(r'(\S+)\s+"?(.*?)"?;', str(gtf_row.attributes)))

@dataclass
class GTFRecord():
    seq_region_name: str
    source: str = field(repr=False, compare=False)
    feature_type: str
    start: int
    end: int
    score: int = field(repr=False, compare=False)
    strand: int
    phase: int
    attributes: dict = field(default_factory=dict, repr=False, compare=False)

    @classmethod
    def fromGTFRow(cls, row: str):
        record = GTFRawRecord._make(row.rstrip().split(GTF_DELIMITERS[0]))
        return cls(seq_region_name = record.seqname,
            source = record.source,
            feature_type = str(record.feature).lower(),
            start = int(record.start),
            end = int(record.end),
            score = record.score,
            strand = 1 if STRAND_MAP.get(record.strand) == 0 else STRAND_MAP.get(record.strand),
            phase = -1 if record.frame == '.' else record.frame,
            attributes = _get_attributes(record)
        )

def _get_feature_parent_id(gtf_row: GTFRecord) -> str:
    if not PARENTS.get(gtf_row.feature_type):
        return None
    return _get_attributes(gtf_row).get(f"{PARENTS.get(gtf_row.feature_type)}_id")

def _get_feature_id(gtf_row: GTFRecord) -> str:
    return _get_attributes(gtf_row).get(f"{gtf_row.feature_type}_id")


@contextlib.contextmanager
def _open_gtf_fileh(filename: Path):
    if not filename.exists():
        raise ValueError(f"Cannot open file {filename}.")
    if not filename.is_file():
        raise ValueError(f"Provided {filename} is not a file.")
    try:
        fh = open(filename, 'rt', 1, encoding="utf-8")
        yield fh
    except OSError as exc:
        raise exc
    finally:
        fh.close()

def _build_transcript(reg_slice: Slice, gtf_tr_row: GTFRecord, exons: list[Exon], analysis_name: str = "gb") -> Transcript:
    tr = Transcript.fastinit(start=gtf_tr_row.start,
                                end=gtf_tr_row.end,
                                length=None,
                                analysis_name=analysis_name,
                                strand=gtf_tr_row.strand,
                                reg_slice=reg_slice,
                                biotype=gtf_tr_row.attributes.get('biotype'))
    tr.stable_id = gtf_tr_row.attributes.get('transcript_id')
    if gtf_tr_row.strand == -1:
        tr.set_exons(sorted(exons, key=lambda e: e.start, reverse=True))
    else:
        tr.set_exons(sorted(exons, key=lambda e: e.start))
    return tr

def parse_gtf(filename: Path) -> list[Feature]:
    reg_slice = Slice.fastinit('chromosome', "1", 20000000, 1, 20000000, 1)
    exons = []
    transcripts = []
    genes = []
    tr_row = None
    first_loop = True
    with _open_gtf_fileh(filename) as gtf_h:
        ii = 0
        for line in gtf_h:
            exon = None
            if line.startswith("#"):
                continue
            gtf_row = GTFRecord.fromGTFRow(line.rstrip())
            if gtf_row.seq_region_name != reg_slice.region.name:
                warnings.warn(f"Found new region {gtf_row.seq_region_name} in file. Stopping parsing.")
                break
            if gtf_row.feature_type not in ('transcript', 'exon'):
                warnings.warn(f"Found unmanaged feature type {gtf_row.feature_type}. Skipping ...")
                continue
            if gtf_row.feature_type == 'transcript':
                if not first_loop:
                    tr = _build_transcript(reg_slice, tr_row, exons)
                    transcripts.append(tr)
                    exons = []
                    tr_row = None
                    ii += 1
                tr_row = gtf_row
                continue
            if gtf_row.feature_type != 'exon':
                continue
            print(f"Exon no: {gtf_row.attributes.get('exon_number')} - Transcript ID: {gtf_row.attributes.get('transcript_id')}")
            exon = Exon.fastinit(start=gtf_row.start,
                                    end=gtf_row.end,
                                    strand=gtf_row.strand,
                                    gen_slice=reg_slice,
                                    phase=gtf_row.phase)
            exons.append(exon)

            first_loop = False

        if len(exons) > 0:
            tr = _build_transcript(reg_slice, tr_row, exons)
            transcripts.append(tr)
    return transcripts

if __name__ == "__main__":
    file = Path("./data/transcriptomic_raw.gtf")
    trs = parse_gtf(file)
    for t in trs:
        print(t)
