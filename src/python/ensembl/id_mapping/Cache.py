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

from sqlalchemy.orm import Session

from ensembl.api.dbsql.SliceAdaptor import SliceAdaptor
from ensembl.api.dbsql.GeneAdaptor import GeneAdaptor
from ensembl.api.dbsql.CoordSystemAdaptor import CoordSystemAdaptor

from ensembl.api.core.Assembly import CoordSystem
from ensembl.api.core.Gene import Gene
from ensembl.id_mapping.TinyGene import TinyGene

from enum import Enum

__all__ = [ 'Cache' ]

class DBType(Enum):
    SOURCE = 1
    TARGET = 2

class Cache():
    """
    ensembl.id_mapping.Cache - a cache to hold data objects used by the 
    id_mapping application
    """
    _cache_method: str = ''
    _common_cs: list[CoordSystem] = []
    _conf = {}

    @classmethod
    def build_cache_by_slice(cls, src_session: Session, tgt_session: Session, dbtype: DBType, slice_name: str):
        """
        Arg[1]      : String dbtype - db type (source|target)
        Arg[2]      : String slice_name - the name of a slice (format as returned by
                      the name property of ensembl.api.core.Slice)
        Example     : (num_genes, filesize) = Cache.build_cache_by_slice(
                      'source', 'chromosome:NCBI36:X:1:1000000:-1');
        Description : Builds a cache of genes, transcripts, translations and exons
                      needed by the IdMapping application and serialises the resulting
                      cache object to a file, one slice at a time.
        Return type : list of the number of genes processed and the size of the
                      serialised cache file
        Exceptions  : thrown on invalid slice name
        Caller      : general
        Status      : At Risk
                    : under development
        """
        # # set cache method (required for loading cache later)
        cls._cache_method('BY_SEQ_REGION')

        if dbtype.name == 'SOURCE':
            session = src_session
        elif dbtype == 'TARGET':
            session = tgt_session
        else:
            raise Exception('Undefined dbtype: must be (SOURCE|TARGET)')

        slice = SliceAdaptor.fetch_by_name(session, slice_name)
        genes = GeneAdaptor.fetch_all_by_Slice(session, slice, load_transcripts=True)

        # find common coord_system
        common_cs_found = Cache.find_common_coord_systems(src_session, tgt_session)
        
        # find out whether native coord_system is a common coord_system.
        # if so, you don't need to project.
        # also don't project if no common coord_system present
        need_project = True
        if cls.is_common_cs(slice.region.coord_system) or not cls.highest_common_cs:
            need_project = False

        # build cache
        mytype = f'{dbtype.name}.{slice.name}'
        num_genes = cls.build_cache_from_genes( mytype, genes, need_project )
        
        # write cache to file, then flush cache to reclaim memory
        size = cls.write_all_to_file(mytype)

        return (num_genes, size)
    
  
    @classmethod
    def find_common_coord_systems(cls, src_session: Session, tgt_session: Session) -> list[CoordSystem]:
        with src_session, tgt_session:
            source_cs = CoordSystemAdaptor.fetch_all_default(src_session)
            target_cs = CoordSystemAdaptor.fetch_all_default(tgt_session)
        
        for s_cs in source_cs:
            for t_cs in target_cs:
                if s_cs.name == t_cs.name:
                    if s_cs.version and s_cs.version != t_cs.version:
                        continue

                    if Cache.seq_regions_compatible(src_session, tgt_session, s_cs):
                        cls._common_cs.append(s_cs)
                        break
        
        return cls._common_cs


    @staticmethod
    def seq_regions_compatible(src_session: Session, tgt_session: Session, source_cs: CoordSystem) -> bool:
        if not isinstance(source_cs, CoordSystem):
           raise Exception(f'You must provide a CoordSystem')
        
        src_seq_regions = SliceAdaptor.fetch_all(src_session, source_cs.name, source_cs.version)
        tgt_seq_regions = SliceAdaptor.fetch_all(tgt_session, source_cs.name, source_cs.version)
        # sanity check to prevent divison by zero
        s_sr_count = len(src_seq_regions)
        t_sr_count = len(tgt_seq_regions)
        if s_sr_count == 0 or t_sr_count == 0:
            return False
        
        sr_match = { sr.name: sr.region.length for sr in src_seq_regions }

        equal = 0.0

        for t_sr in tgt_seq_regions:
            if sr_match.get(t_sr.name):
                equal += 1
                if sr_match.get(t_sr.name) != t_sr.region.length:
                    return False
        
        if equal/s_sr_count > 0.5 and equal/t_sr_count > 0.5:
            return True
        else:
            # logger
            return False
        
    
    @classmethod
    def is_common_cs(cls, cs: CoordSystem) -> bool:
        # We defined two CS being equal, when they have both name and version in common
        # See ensembl.api.core.Assembly.CoordSystem class
        if cs == cls._common_cs:
            return True
        return False
    
    
    @classmethod
    def highest_common_cs(cls) -> CoordSystem:
        # We defined two CS being equal, when they have both name and version in common
        # See ensembl.api.core.Assembly.CoordSystem class
        if cls._common_cs:
            return cls._common_cs[0]
        return None
    

    @classmethod
    def cache_method(cls, cm: str = None) -> str:
        if cm:
            cls._cache_method = cm
            return None
        return cls._cache_method
    
    @classmethod
    def conf(cls, conf_values: dict) -> dict:
        if conf_values:
            cls._conf = conf_values
            return None
        return cls._conf
    

    @classmethod
    def build_cache_from_genes(cls, mytype: str, genes: list[Gene], need_project: bool):
        if not mytype:
            raise ValueError('You must provide a type')
        if not genes or not isinstance(genes, list) or not isinstance(genes[0], Gene):
            raise ValueError('You must provide a non-empty list of genes')
        
        if cls._conf.get('biotypes') or cls._conf.get('biotypes_include') or cls._conf.get('biotypes_exclude'):
            genes = cls.filter_biotypes(genes)
        num_genes = len(genes)

        # initialise cache for the given type.
        # $self->{'cache'}->{$type} = {};

        # loop over genes sorted by gene location.
        # the sort will hopefully improve assembly mapper cache performance and
        # therefore speed up exon sequence retrieval
        genes.sort(key=lambda x: x.get_slice().location.start)

        for gene in genes:
            if need_project == 'CHECK':
                # find out whether native coord_system is a common coord_system.
                # if so, you don't need to project.
                # also don't project if no common coord_system present
                if cls.highest_common_cs:
                    if cls.is_common_cs(gene.get_slice().region.coord_system):
                        need_project = False
                else:
                    need_project = True
        
            # create lightweight gene
            lgene = TinyGene(gene.internal_id,
                             gene.unversioned_stable_id,
                             gene.version,
                             gene.get_metadata('created_date'),
                             gene.get_metadata('modified_date'),
                             gene.get_slice().location.start,
                             gene.get_slice().location.end,
                             gene.get_slice().strand.value,
                             gene.get_slice().name,
                             gene.biotype)

            # build gene caches
            # cls.add( 'genes_by_id', $type, $gene->dbID, $lgene );

            # create lightweight transcripts
            lgene.add_fat_transcripts(gene.get_transcripts())

            # build transcript caches
            # $self->add( 'transcripts_by_id',      $type, $tr->dbID, $ltr );
            # $self->add( 'genes_by_transcript_id', $type, $tr->dbID, $lgene );
