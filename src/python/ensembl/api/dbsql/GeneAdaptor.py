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

from sqlalchemy import and_, select
from sqlalchemy.orm import Session, Bundle
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.engine.row import Row

from ensembl.core.models import Gene as GeneORM, GeneAttrib as GeneAttribORM, AttribType as AttribTypeORM, SeqRegion as SeqRegionORM, CoordSystem as CoordSystemORM
from ensembl.core.models import Xref, Analysis

from typing import Optional

from ensembl.api.core import Gene, Strand
from ensembl.api.core.Slice import Slice
from ensembl.api.dbsql.SliceAdaptor import SliceAdaptor
from ensembl.api.dbsql.TranscriptAdaptor import TranscriptAdaptor
from ensembl.api.dbsql.BiotypeAdaptor import BiotypeAdaptor

import warnings
# from ensembl.api.dbsql.Utils import timeme

__all__ = [ 'GeneAdaptor' ]

class GeneAdaptor():
    """Contains all the Gene related functions over Gene ORM
    """

    @classmethod
    def fetch_by_stable_id(cls, session: Session, stable_id: str, load_transcripts: bool = False,
                           load_exons: bool = False, load_translations: bool = False) -> Gene:
        """
        Arg [1]    : Session session - the ORM session object to connect to the DB to
        Arg [2]    : Str stable_id 
                     The stable ID of the gene to retrieve
        Arg [3]    : bool load_transcripts (default False)
                     flag to load gene's transcripts
        Arg [4]    : bool load_exons (default False)
                     flag to load transcript's exons. Ignored unless load_transcripts=True.
        Arg [5]    : bool load_translations (default False)
                     flag to load transcript's translations. Ignored unless load_transcripts=True.
        Example    : gene = GeneAdaptor.fetch_by_stable_id('ENSG00000148944');
        Description: Retrieves a gene object from the database via its stable id.
                     The gene will be retrieved in its native coordinate system (i.e.
                     in the coordinate system it is stored in the database). It may
                     be converted to a different coordinate system through a call to
                     transform() or transfer(). If the gene or exon is not found
                     undef is returned instead.
        Returntype : ensembl.api.core.Gene
        Exceptions : if we cant get the gene in given coord system
        Caller     : general
        Status     : At Risk
                   : under development
        """
        # my $constraint = "g.stable_id = ? AND g.is_current = 1";
        # $self->bind_param_generic_fetch($stable_id, SQL_VARCHAR);
        # my ($gene) = @{$self->generic_fetch($constraint)};
        if stable_id.find('.') > 0:
            (unversioned_stable_id, version) = stable_id.split('.')
        else:
            unversioned_stable_id = stable_id
            version = 0

        if not isinstance(version, int):
            raise ValueError('Stable ID version must be an int')
        
        res = (
            session.query(
                GeneORM,
                Bundle("Xref",
                       Xref.dbprimary_acc,
                       Xref.description
                ),
                Bundle("Analysis",
                       Analysis.logic_name
                )
            )
            .join(GeneORM.analysis)
            .join(GeneORM.display_xref, isouter=True)
            .filter(GeneORM.stable_id == unversioned_stable_id)
            .filter(GeneORM.is_current == '1')
            .first()
        )

        if not res:
                raise NoResultFound(f'Could not find {unversioned_stable_id} as current stable_id')
        if version > 0 and res.Gene.version != version:
            warnings.warn(f"Could not find {stable_id}. Found {unversioned_stable_id}.{res.Gene.version} instead!", UserWarning)

        gene_attribs = (
            session.query(GeneAttribORM, AttribTypeORM)
            .join(GeneAttribORM.attribs)
            .filter(GeneAttribORM.gene_id == res.Gene.gene_id)
            .all()
        )

        gene = GeneAdaptor._generow_to_gene(session, res, gene_attribs)
        if load_transcripts:
                tr = TranscriptAdaptor.fetch_all_by_gene_id(session, res.Gene.gene_id, load_exons, load_translations)
                gene.set_transcripts(tr)

        return gene
    
    @classmethod
    def fetch_all_by_Slice(cls,
                           session: Session,
                           slice: Slice,
                           logic_name: Optional[str] = None,
                           load_transcripts: Optional[bool] = False,
                           source: Optional[str] = None,
                           biotype: Optional[str] = None,
                           load_exons: bool = False,
                           load_translations: bool = False
                           ) -> list[Gene]:
        """
        Arg [1]    : Session session - the ORM session object to connect to the DB to
        Arg [2]    : ensembl.api.core.Slice slice
                     The slice to fetch genes on.
        Arg [3]    : (optional) str logic_name
                     the logic name of the type of features to obtain
        Arg [4]    : (optional) bool load_transcripts
                     if true, transcripts will be loaded immediately rather than
                     lazy loaded later.
        Arg [5]    : (optional) str source
                     the source name of the features to obtain.
        Arg [6]    : (optional) string biotype
                     the biotype of the features to obtain.
        Arg [7]    : bool load_exons (default False)
                     flag to load transcript's exons. Ignored unless load_transcripts=True.
        Arg [8]    : bool load_translations (default False)
                     flag to load transcript's translations. Ignored unless load_transcripts=True.
        Example    : genes = GeneAdaptor.fetch_all_by_Slice()
        Description: Overrides superclass method to optionally load transcripts
                     immediately rather than lazy-loading them later.  This
                     is more efficient when there are a lot of genes whose
                     transcripts are going to be used.
        Returntype : reference to list of genes 
        Exceptions : thrown if exon cannot be placed on transcript slice
        Caller     : aSlice.get_all_Genes
        Status     : At Risk
                   : under development
        """
        stmt = (select(
                GeneORM,
                Bundle("Xref",
                       Xref.dbprimary_acc,
                       Xref.description
                ),
                Bundle("Analysis",
                       Analysis.logic_name
                ),
                Bundle("SeqRegion",
                        SeqRegionORM.seq_region_id,
                        SeqRegionORM.name,
                        SeqRegionORM.length
                ),
                Bundle("CoordSystem",
                        CoordSystemORM.coord_system_id,
                        CoordSystemORM.name,
                        CoordSystemORM.version,
                        CoordSystemORM.rank
                )
            )
            .join(GeneORM.seq_region)
            .join(SeqRegionORM.coord_system)
            .join(GeneORM.analysis)
            .join(GeneORM.display_xref, isouter=True)
            .where(and_(GeneORM.is_current == '1',
                        SeqRegionORM.name == slice.seq_region_name,
                        CoordSystemORM.name == slice.coord_system.name,
                        CoordSystemORM.version == slice.coord_system.version,
                        GeneORM.seq_region_start >= slice.seq_region_start,
                        GeneORM.seq_region_end <= slice.seq_region_end,
                        GeneORM.seq_region_strand == slice.seq_region_strand.value
                        ))
        )
        if logic_name:
            stmt = stmt.filter(Analysis.logic_name == logic_name)
        if source:
            stmt = stmt.filter(GeneORM.source == source)
        if biotype:
            stmt = stmt.filter(GeneORM.biotype == biotype)
        
        rows = session.execute(stmt).all()

        if not rows:
                return []
        
        genes = []
        for row in rows:
            gene_attribs = (
                session.query(GeneAttribORM, AttribTypeORM)
                .join(GeneAttribORM.attribs)
                .filter(GeneAttribORM.gene_id == row.Gene.gene_id)
                .all()
            )
            gene = GeneAdaptor._generow_to_gene(session, row, gene_attribs, slice)
            if load_transcripts:
                tr = TranscriptAdaptor.fetch_all_by_gene_id(session, row.Gene.gene_id, load_exons, load_translations)
                gene.set_transcripts(tr)
            
            genes.append(gene)

        return genes
        

    @classmethod
    def fetch_all(cls, session: Session,
                  species_id: int = 1,
                  logic_name: str = None,
                  source: str = None,
                  biotype: str = None,
                  load_transcripts: bool = False,
                  load_exons: bool = False,
                  load_translations: bool = False) -> tuple:
        # my $constraint = 'g.biotype != "LRG_gene" and g.is_current = 1';
        stmt = (select(
                GeneORM,
                Bundle("Xref",
                       Xref.dbprimary_acc,
                       Xref.description
                ),
                Bundle("Analysis",
                       Analysis.logic_name
                ),
                Bundle("SeqRegion",
                        SeqRegionORM.seq_region_id,
                        SeqRegionORM.name,
                        SeqRegionORM.length
                ),
                Bundle("CoordSystem",
                        CoordSystemORM.coord_system_id,
                        CoordSystemORM.name,
                        CoordSystemORM.version,
                        CoordSystemORM.rank
                )
            )
            .join(GeneORM.seq_region)
            .join(SeqRegionORM.coord_system)
            .join(GeneORM.analysis)
            .join(GeneORM.display_xref, isouter=True)
            .where(GeneORM.is_current == '1')
            .where(CoordSystemORM.species_id == species_id)
        )
        if logic_name:
            stmt = stmt.filter(Analysis.logic_name == logic_name)
        if source:
            stmt = stmt.filter(GeneORM.source == source)
        if biotype:
            stmt = stmt.filter(GeneORM.biotype == biotype)
        
        rows = session.execute(stmt).all()

        if not rows:
                return []
        
        genes = []
        for row in rows:
            gene_attribs = (
                session.query(GeneAttribORM, AttribTypeORM)
                .join(GeneAttribORM.attribs)
                .filter(GeneAttribORM.gene_id == row.Gene.gene_id)
                .all()
            )
            gene = GeneAdaptor._generow_to_gene(session, row, gene_attribs, slice)
            if load_transcripts:
                tr = TranscriptAdaptor.fetch_all_by_gene_id(session, row.Gene.gene_id, load_exons, load_translations)
                gene.set_transcripts(tr)
            
            genes.append(gene)

        return genes
    

    @classmethod
    def list_seq_region_ids(cls, session: Session, species_id: int) -> list[int]:
        
        stmt = (select(GeneORM.seq_region_id)
                .join(SeqRegionORM.coord_system)
                .join(GeneORM.seq_region)
                .where(CoordSystemORM.species_id == species_id)
        )
        return session.scalars(stmt).all()
    

    @classmethod
    def _generow_to_gene(cls, session: Session, g_row: Row, gene_attribs: list[Row], slice: Slice = None) -> Gene:
        
        if not slice:
            slice = SliceAdaptor.fetch_by_seq_region_id(session, g_row.Gene.seq_region_id)
        biotype = BiotypeAdaptor.fetch_by_name_object_type(session, g_row.Gene.biotype, 'gene')

        #if slice != g_slice: throw or extend?

        g = Gene(
            g_row.Gene.stable_id,
            g_row.Gene.version,
            g_row.Gene.gene_id,
            slice,
            g_row.Gene.seq_region_start,
            g_row.Gene.seq_region_end,
            Strand(g_row.Gene.seq_region_strand),
            g_row.Analysis.logic_name,
            transcripts=None,
            biotype=biotype,
            source=g_row.Gene.source,
            external_name=g_row.Xref.dbprimary_acc,
            description=g_row.Xref.description,
            created_date=g_row.Gene.created_date,
            modified_date=g_row.Gene.modified_date
        )        

        for ga in gene_attribs:
            g.add_attrib(ga.AttribType.code, ga.GeneAttrib.value)
                
        return g