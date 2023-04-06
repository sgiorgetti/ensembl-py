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

from sqlalchemy import select, func, bindparam
from sqlalchemy.orm import Session, Bundle, aliased
from Bio.Seq import Seq, MutableSeq, reverse_complement

from ensembl.core.models import Dna as DnaORM, SeqRegionAttrib as SeqRegionAttribORM, SeqRegion as SeqRegionORM, AttribType as AttribTypeORM, CoordSystem as CoordSystemORM

from ensembl.api.dbsql.AssemblyAdaptor import AssemblyAdaptor, CoordSystemAdaptor
from ensembl.api.dbsql.MetaAdaptor import MetaAdaptor

from ensembl.api.core.Slice import Slice, Strand

from dataclasses import dataclass
from functools import lru_cache
import warnings

__all__ = [ 'SequenceAdaptor', 'SeqEdit' ]


class SequenceAdaptor():
    def __init__(self, session: Session, species_id: int = None) -> None:
        if not session or not isinstance(session, Session):
            raise ValueError('Need a session to a DB')
        if species_id and not isinstance(species_id, int):
            raise ValueError(f'Species ID must be None or a single int.')
        self._session = session
        self._species_id = species_id
        self._meta = MetaAdaptor(session)

        if not self._meta.is_multispecies():
            self._species_id = 1
            self._species_prod_name = self._meta.get_species_prod_name_by_id(self._species_id)
        # the multispecies db mgmt is missing here!!!
        self._rna_edits_cache()

    def __repr__(self) -> str:
        if self._species_id:
            return f'{self.__class__.__name__}({self._session.connection().engine.url}:{self._species_id})'
        return f'{self.__class__.__name__}({self._session.connection().engine.url}'
    
    @property
    def species_id(self) -> int:
        return self._species_id
    

    @lru_cache(maxsize=25)
    def _fetch_seq(self, seq_region_id: int) -> str:
        stmt = (
            select(func.upper(DnaORM.sequence).label("seq"))
            .where(DnaORM.seq_region_id == seq_region_id)
        )
        seq = self._session.scalars(stmt).first()
        return seq
    
    def _fetch_raw_seq(self, seq_region_id: int, start: int, length: int) -> Seq:
        seq = self._fetch_seq(seq_region_id)[start:start+length]
        return Seq(seq)

    def fetch_by_Slice_start_end_strand(self, slice: Slice, start: int, end: int, strand: Strand) -> Seq:
        if not slice or not isinstance(slice, Slice):
            raise ValueError(f'Slice argument is required.')

        start = start if start else 1
        strand = strand if strand else Strand.FORWARD

        # This is both verbatim Perl and horrible
        if not end or start > end or start < 0 or end < 0 or slice.start > slice.end or slice.is_circular():
            if not end or start > end:
                return self._fetch_by_Slice_start_end_strand_circular(slice, start, end, strand)
            if end and end < 0:
                end += slice.seq_region_length
            if start < 0:
                start += slice.seq_region_length
            if slice.start > slice.end:
                return self._fetch_by_Slice_start_end_strand_circular(slice, start, end, strand)
        if not end and not slice.is_circular():
            end = slice.end - slice.start + 1
        if start > end:
            raise ValueError('Start must be less than or equal to end.')
        
        
        #get a new slice that spans the exact region to retrieve dna from
        right_expand  = end - slice.length() #negative is fine
        left_expand   = 1 - start #negative is fine
        if right_expand or left_expand:
            slice = slice.expand(left_expand, right_expand)

        #### Projecting and trimming
        #haplotypes and PARs are not supported
        if AssemblyAdaptor.check_assembly_exceptions(self._session):
            raise NotImplementedError('Haplotypes and PARs are not supported.')
        
        # we need to project this slice onto the sequence coordinate system
        # even if the slice is in the same coord system, we want to trim out
        # flanking gaps (if the slice is past the edges of the seqregion)
        cs_seq = CoordSystemAdaptor.fetch_sequence_level(self._session)
        projections = SliceAdaptor.project(self._session, slice, cs_seq)

        seq = ''
        total = 0
        tmp_seq = ''
        for segment in projections:
            (p_start, p_end, seq_slice) = segment

            gap = p_start - total - 1
            seq += 'N' * gap if gap else ''

            tmp_seq = self._fetch_raw_seq(seq_slice.seq_region_id,
                                          seq_slice.start,
                                          seq_slice.length())
            if not tmp_seq:
                raise Exception(f'No sequence found for seq_region {seq_slice.seq_region_id}:{seq_slice.start}')
            
            if seq_slice.strand == Strand.REVERSE:
                tmp_seq = reverse_complement(tmp_seq)

            seq += tmp_seq
            total = p_end
        
        # check for any remaining gaps at the end
        gap = slice.length() - len(seq)
        seq += 'N' * gap if gap else ''
        
        if self._rna_edits_cache and self._rna_edits_cache.get(slice.seq_region_id):
            # $self->_rna_edit($slice,\$seq);
            seq = self._rna_edit(slice, seq)

        if strand == Strand.REVERSE:
                seq = reverse_complement(seq)
        
        return seq
    

    def _fetch_by_Slice_start_end_strand_circular(self, slice: Slice, start: int, end: int, strand: Strand) -> Seq:
        if not slice or not isinstance(slice, Slice):
            raise ValueError(f'Slice argument is required.')
        if not start:
            start = 1
        if not strand:
            strand = Strand.FORWARD
        if not end:
            end = slice.end = slice.start + 1

        # Similar to Perl - weak Kung Fu - also lacking modular arithmetic (like always)
        # Do not like it at all - to redo later on
        if start > end and slice.is_circular():
            zeropoint = slice.seq_region_length - slice.start + 1
            seq1 = self._fetch_by_Slice_start_end_strand_circular(slice, 1, zeropoint, Strand.FORWARD)
            seq2 = self._fetch_by_Slice_start_end_strand_circular(slice, zeropoint+1, slice.length(), Strand.FORWARD)
            if strand == Strand.FORWARD:
                return seq1 + seq2
            else:
                return reverse_complement(seq2 + seq1) 

        right_expand  = end - slice.length() #negative is fine
        left_expand   = 1 - start #negative is fine
        if right_expand or left_expand:
            if strand == Strand.FORWARD:
                slice = slice.expand(left_expand, right_expand)
            else:
                slice = slice.expand(right_expand, left_expand)

        #### Projecting and trimming
        #haplotypes and PARs are not supported
        if AssemblyAdaptor.check_assembly_exceptions(self._session):
            raise NotImplementedError('Haplotypes and PARs are not supported.')
        
        # we need to project this slice onto the sequence coordinate system
        # even if the slice is in the same coord system, we want to trim out
        # flanking gaps (if the slice is past the edges of the seqregion)
        cs_seq = CoordSystemAdaptor.fetch_sequence_level(self._session)
        projections = SliceAdaptor.project(self._session, slice, cs_seq)

        seq = ''
        total = 0
        tmp_seq = ''
        for segment in projections:
            (p_start, p_end, seq_slice) = segment

            gap = p_start - total - 1
            seq += 'N' * gap if gap else ''

            tmp_seq = self._fetch_raw_seq(seq_slice.seq_region_id,
                                          seq_slice.start,
                                          seq_slice.length())
            if not tmp_seq:
                raise Exception(f'No sequence found for seq_region {seq_slice.seq_region_id}:{seq_slice.start}')
            
            if seq_slice.strand == Strand.REVERSE:
                tmp_seq = reverse_complement(tmp_seq)

            seq += tmp_seq
            total = p_end
        
        # check for any remaining gaps at the end
        gap = slice.length() - len(seq)
        seq += 'N' * gap if gap else ''
        
        # self._rna_edits_cache ????
        if self._rna_edits and self._rna_edits.get(slice.seq_region_id):
            seq = self._rna_edit(slice, seq)

        if strand == Strand.REVERSE:
                seq = reverse_complement(seq)
        
        return seq
    
    def _rna_edits_cache(self) -> None:
        if not self._species_id and self._meta.is_multispecies():
            warnings.warn(f'DB is multispecies: consider to set species_id for better performance', UserWarning)
        if self._meta.is_multispecies():
            subq = (select(SeqRegionAttribORM.seq_region_id, SeqRegionAttribORM.value)
                .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
                .where(AttribTypeORM.code == '_rna_edit')
                .subquery()
                )
            sr_attr_subq = aliased(SeqRegionAttribORM, subq, name="seq_region_attrib")
            stmt = (select(
                    Bundle("SRAttrRNAEdits", sr_attr_subq.seq_region_id, sr_attr_subq.value),
                    Bundle("CoordSystem", CoordSystemORM.species_id)
                    )
            .join(SeqRegionORM.coord_system)
            .join(sr_attr_subq, SeqRegionORM.seq_region_id == sr_attr_subq.seq_region_id, isouter=True)
            )
            if self._species_id:
                stmt = stmt.where(CoordSystemORM.species_id == self._species_id)
            stmt = stmt.order_by(CoordSystemORM.species_id)
        else:
            stmt = (select(
                    bindparam("species_id", self._species_id),
                    SeqRegionAttribORM.seq_region_id,
                    SeqRegionAttribORM.value.label("seq_edit")
                )
                .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
                .where(AttribTypeORM.code == '_rna_edit')
                )
        
        #  Very un-pythonic - bleah!
        d = {}
        for r in self._session.execute(stmt).all():
            if d.get(r.species_id):
                if d.get(r.species_id).get(r.seq_region_id):
                    d[r.species_id][r.seq_region_id].append(SeqEdit.fromString(r.seq_edit))
                else:
                    d[r.species_id][r.seq_region_id] = [SeqEdit.fromString(r.seq_edit),]
            else:
                d[r.species_id] = { r.seq_region_id: [SeqEdit.fromString(r.seq_edit),] }
        
        self._rna_edits_cache = d
    

    def _rna_edit(self, slice: Slice, seq: Seq) -> Seq:
        s_start = slice.start
        s_end = s_start + len(seq) - 1
        my_seq = MutableSeq(seq)

        for seq_edit in self._rna_edits_cache[slice.coord_system.species_id][slice.seq_region_id]:
            # check that RNA edit is not outside the requested region : happens quite often with LRG regions
            if seq_edit.end < s_start or s_end < seq_edit.start:
                continue

            edit_offset = max(s_start - seq_edit.start, 0)
            edit_len = seq_edit.length()

            # If the edit isn't fully encompassed by the slice, we need to extract the
            # edit's sub-sequence that we're patching on to the sequence
            if seq_edit.start < s_start or seq_edit.end > s_end:
                edit_len = seq_edit.length() - edit_offset - max(seq_edit.end - s_end, 0)

            my_seq[edit_offset:edit_len] = seq_edit.edit[edit_offset:edit_len]

        return Seq(my_seq)



@dataclass
class SeqEdit():
    start: int
    end: int
    edit: str

    @classmethod
    def fromString(cls, seq_edit: str):
        if not seq_edit:
            raise ValueError(f'Missing argument seq edit')
        (s,e,edit) = seq_edit.split()
        return cls(int(s),int(e),edit)

    def length(self) -> int:
        return self.end - self.start + 1




from ensembl.database.dbconnection import DBConnection
from ensembl.api.dbsql import SliceAdaptor
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import time
# st = time.process_time()

def main():
    dbc = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
    # dbc = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/bos_taurus_core_110_12')
    with dbc.session_scope() as session:
        slice_name = 'chromosome:GRCh38:13:32315086:32400268:1'
        # slice_name = 'primary_assembly:ARS-UCD1.2:7:16782097:16851987:-1'
        slice = SliceAdaptor.fetch_by_name(session, slice_name)
        sa = SequenceAdaptor(session)
        bioseq = sa.fetch_by_Slice_start_end_strand(slice, 
                                                    slice.start-slice.seq_region_start,
                                                    slice.length(), slice.strand)
        bioseq_rev = reverse_complement(bioseq)
        print(f'Found sequence {len(bioseq)}bp long')
        # with open('my_sequence.fa', 'w') as fh:
        #     if slice.strand.value == -1:
        #         record = SeqRecord(bioseq_rev, id=f'{slice.seq_region_name}', name='dna:primary_assembly',
        #                            description=f'{slice.name()}')
        #     else:
        #         record = SeqRecord(bioseq, id=f'{slice.seq_region_name}', name='dna:primary_assembly',
        #                            description=f'{slice.name()}')
        #     SeqIO.write(record, fh, "fasta")


if __name__ == '__main__':
    main()