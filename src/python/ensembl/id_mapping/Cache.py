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
from ensembl.api.core.Strand import Strand
from ensembl.id_mapping.TinyFeature import TinyFeature
from ensembl.id_mapping.TinyGene import TinyGene
from ensembl.id_mapping.TinyTranslation import TinyTranslation
from ensembl.id_mapping.Utils import format_bytes, timeme

from typing import Union
from pathlib import Path, PurePath
import warnings
import pickle

__all__ = [ 'Cache' ]

class Cache():
    """
    ensembl.id_mapping.Cache - a cache to hold data objects used by the 
    id_mapping application
    """
    _conf = { 'basedir' : '/home/stefano/ensembl-py/scripts/id_mapping' }
    # _cache has the following format
    # c_item = { cache_type: { cache_name: { cache_key: [value,] } } }
    _cache = {}
    _instance = {'cache_method': '',
                 'loaded': [],
                 'hccs': '',
                 'ccs': []}
    _dump_path = None


    @classmethod
    def build_cache_all(cls, src_session: Session, tgt_session: Session, dbtype: str) -> tuple[int,str]:
        # set cache method (required for loading cache later)
        cls.cache_method('ALL')
        
        dbtype = dbtype.lower()

        if dbtype == 'source':
            session = src_session
        elif dbtype == 'target':
            session = tgt_session
        else:
            raise Exception('Undefined dbtype: must be (source|target)')
        
        genes = GeneAdaptor.fetch_all(session, load_transcripts=True,
                                      load_exons=True, load_translations=True)

        # find common coord_system
        common_cs_found = Cache.find_common_coord_systems(src_session, tgt_session)

        # build cache
        mytype = f'{dbtype}.ALL'
        num_genes = cls.build_cache_from_genes( mytype, genes, 'CHECK' )
        
        # write cache to file, then flush cache to reclaim memory
        size_str = cls.write_all_to_file(mytype)

        return (num_genes, size_str)
        


    @classmethod
    @timeme
    def build_cache_by_slice(cls, src_session: Session, tgt_session: Session, dbtype: str, slice_name: str) -> tuple[int,str]:
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
        # set cache method (required for loading cache later)
        cls.cache_method('BY_SEQ_REGION')

        mytype = dbtype.lower()

        if mytype == 'source':
            session = src_session
        elif mytype == 'target':
            session = tgt_session
        else:
            raise Exception('Undefined dbtype: must be (source|target)')

        slice = SliceAdaptor.fetch_by_name(session, slice_name)
        if not slice:
            raise Exception(f'Could not retrieve slice {slice_name}.')
        
        print(f'Got slice {slice.name()}')  
        genes = GeneAdaptor.fetch_all_by_Slice(session, slice, load_transcripts=True,
                                               load_exons=True, load_translations=True)
        print(f'Found {len(genes)} genes')

        # find common coord_system
        common_cs_found = Cache.find_common_coord_systems(src_session, tgt_session)
        
        # find out whether native coord_system is a common coord_system.
        # if so, you don't need to project.
        # also don't project if no common coord_system present
        need_project = True
        if cls.is_common_cs(slice.coord_system) or not cls.highest_common_cs:
            need_project = False

        # build cache
        mytype = f'{mytype}.{slice.name()}'
        num_genes = cls.build_cache_from_genes( mytype, genes, need_project )
        
        # write cache to file, then flush cache to reclaim memory
        size_str = cls.write_all_to_file(mytype)

        return (num_genes, size_str)
    
  
    @classmethod
    @timeme
    def find_common_coord_systems(cls, src_session: Session, tgt_session: Session) -> list[CoordSystem]:
        source_cs = CoordSystemAdaptor.fetch_all_default(src_session)
        target_cs = CoordSystemAdaptor.fetch_all_default(tgt_session)
        
        for s_cs in source_cs:
            for t_cs in target_cs:
                if s_cs.name == t_cs.name:
                    if s_cs.version and s_cs.version != t_cs.version:
                        continue
                    
                    if Cache.seq_regions_compatible(src_session, tgt_session, s_cs):
                        cls._instance['ccs'].append(s_cs)
                        if len(cls._instance['ccs']) == 1:
                            cls._instance['hccs'] = s_cs
                        break
        
        return cls._instance['ccs']


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
        
        sr_match = { sr.seq_region_name: sr.length() for sr in src_seq_regions }

        equal = 0.0

        for t_sr in tgt_seq_regions:
            if sr_match.get(t_sr.seq_region_name):
                equal += 1
                if sr_match.get(t_sr.seq_region_name) != t_sr.length():
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
        if cs in cls._instance.get('ccs'):
            return True
        return False
    
    
    @classmethod
    def highest_common_cs(cls) -> CoordSystem:
        return cls._instance.get('hccs')
    
    @classmethod
    def highest_common_cs_name(cls) -> str:
        return cls._instance.get('hccs').name
    
    @classmethod
    def highest_common_cs_version(cls) -> str:
        return cls._instance.get('hccs').version
    

    @classmethod
    def cache_method(cls, cm: str = None) -> str:
        if cm:
            cls._instance['cache_method'] = cm
            return None
        return cls._instance.get('cache_method')
    
    @classmethod
    def conf(cls, conf_values: dict) -> dict:
        if conf_values:
            cls._conf = conf_values
            return None
        return cls._conf
    

    @classmethod
    @timeme
    def build_cache_from_genes(cls, mytype: str, genes: list[Gene], need_project: bool) -> int:
        """
        Arg[1]      : mytype: str - cache type e.g. (source|target)
        Arg[2]      : genes: list[ensembl.api.core.Gene] - genes to build cache from
        Arg[3]      : need_project: bool - indicate if we need to project exons to
                      common coordinate system
        Example     : Cache.build_cache_from_genes(
                      'source.chromosome:NCBI36:X:1:100000:1', genes);
        Description : Builds the cache by fetching transcripts, translations and exons
                      for a list of genes from the database, and creating lightweight
                      ensembl.id_mapping.TinyFeature objects containing only the
                      data needed by the IdMapping application. These objects are
                      attached to a name cache in this cache object. Exons only need
                      to be projected to a commond coordinate system if their native
                      coordinate system isn't common to source and target assembly
                      itself.
        Return type : int - number of genes after filtering
        Exceptions  : thrown on wrong or missing arguments
        Caller      : internal
        Status      : At Risk
                    : under development
        """
        if not mytype:
            raise ValueError('You must provide a type')
        if not genes or not (isinstance(genes, list) or isinstance(genes, tuple)) or not isinstance(genes[0], Gene):
            raise ValueError('You must provide a non-empty list of genes')
        
        genes = list(genes)
        if cls._conf.get('biotypes') or cls._conf.get('biotypes_include') or cls._conf.get('biotypes_exclude'):
            genes = cls.filter_biotypes(genes)
        num_genes = len(genes)

        # initialise cache for the given type.
        # mytype = '(source|target).{slice.seq_region_name}'
        cls._cache[mytype.lower()] = {}

        # loop over genes sorted by gene location.
        # the sort will hopefully improve assembly mapper cache performance and
        # therefore speed up exon sequence retrieval
        genes.sort(key=lambda x: x.get_slice().start)

        for gene in genes:
            if need_project == 'CHECK':
                # find out whether native coord_system is a common coord_system.
                # if so, you don't need to project.
                # also don't project if no common coord_system present
                if cls.highest_common_cs:
                    if cls.is_common_cs(gene.get_slice().coord_system):
                        need_project = False
                else:
                    need_project = True
        
            # create lightweight gene
            lgene = TinyGene(gene.internal_id,
                             gene.stable_id,
                             gene.version,
                             gene.created_date,
                             gene.modified_date,
                             gene.start,
                             gene.end,
                             gene.strand,
                             gene.get_slice().seq_region_name,
                             gene.biotype.name)

            # build gene caches
            cls.add( 'genes_by_id', mytype, gene.internal_id, lgene )

            # create lightweight transcripts
            for tr in gene.get_transcripts():
                ltr = lgene.add_fat_transcript(tr)

                # build transcript caches
                cls.add( 'transcripts_by_id', mytype, tr.internal_id, ltr )
                cls.add( 'genes_by_transcript_id', mytype, tr.internal_id, lgene )

                # translation (if there is one)
                if tr.translation:
                    tl = tr.translation
                    ltl = TinyTranslation(
                        tl.internal_id,
                        tl.stable_id,
                        tl.version,
                        tl.created_date,
                        tl.modified_date,
                        tr.internal_id,
                        tl.seq
                    )

                    ltr.add_translation(ltl)

                    cls.add('translations_by_id', mytype, tl.internal_id, ltl)

                #exons
                for exon in tr.get_exons():
                    lex = ltr.add_fat_exon(exon, need_project)
                    if need_project:
                        #  exon.project (highest CCS)
                        # TO DO
                        raise NotImplementedError()

                    cls.add( 'exons_by_id', mytype, exon.internal_id, lex )
                    cls.add_list( 'transcripts_by_exon_id', mytype, exon.internal_id, ltr )

        return num_genes


    @classmethod
    def add(cls, cache_name: str, cache_type: str, cache_key: str, value: TinyFeature):
        """
        Arg[1]      : cache_name: str - a cache name (e.g. 'genes_by_id')
        Arg[2]      : cache_type: str - a cache type (e.g. "source.$slice_name")
        Arg[3]      : cache_key: str - key of this entry (e.g. a gene dbID)
        Arg[4]      : value: ensembl.id_mapping.TinyFeature - value to cache
        Example     : Cache.add('genes_by_id',
                        'source.chromosome:NCBI36:X:1:1000000:1', '1234', tiny_gene);
        Description : Adds a TinyFeature object to a named cache.
        Return type : ensembl.id_mapping.TinyFeature
        Exceptions  : thrown on wrong or missing arguments
        Caller      : internal
        Status      : At Risk
                    : under development
        """
        if not cache_name:
            raise ValueError(f'You must provide a cache name (e.g. genes_by_id).')
        if not cache_type:
            raise ValueError(f'You must provide a cache type.')
        if not cache_key:
            raise ValueError(f'You must provide a cache key (e.g. a gene dbID).')
        if not value:
            raise ValueError(f'You must provide a TinyFeature object as value')
        if (isinstance(value, list) or isinstance(value, tuple)) and not isinstance(value[0], TinyFeature):
            raise ValueError(f'You must provide a TinyFeature object/list as value')
        elif not isinstance(value, TinyFeature):
            raise ValueError(f'You must provide a TinyFeature object as value, instead of a {type(value)}.')
        
        c_item = { cache_type: { cache_name: { cache_key: [value,]}}}
        if not cls._cache.get(cache_type):
            cls._cache.update(c_item)
        elif not cls._cache.get(cache_type).get(cache_name):
            cls._cache[cache_type].update(c_item[cache_type])
        elif not cls._cache.get(cache_type).get(cache_name).get(cache_key):
            cls._cache[cache_type][cache_name].update(c_item[cache_type][cache_name])
        else:
            cls._cache.get(cache_type).get(cache_name).get(cache_key).append(value)
         
        return cls._cache.get(cache_type).get(cache_name).get(cache_key)
    
    
    @classmethod
    def add_list(cls, cache_name: str, cache_type: str, cache_key: str, value_list: Union[list,tuple]):
        """
        Synonym to Cache.add provided for backward compatibility only
        Caller      : internal
        Status      : At Risk
                    : under development
        """
        return cls.add(cache_name, cache_type, cache_key, value_list)
    

    @classmethod
    def filter_biotypes(cls, genes: list[Gene]) -> list[Gene]:
        """
        Arg[1]      : genes: list[ensembl.api.core.Gene] - the genes to filter
        Example     : filtered = Cache.filter_biotypes(genes)
        Description : Filters a list of genes by biotype.  Biotypes are
                      taken from the IdMapping configuration parameter
                      'biotypes_include' or 'biotypes_exclude'.
                      If the configuration parameter 'biotypes_exclude' is
                      defined, then rather than returning the genes whose
                      biotype is listed in the configuration parameter
                      'biotypes_include' the method will return the genes
                      whose biotype is *not* listed in the 'biotypes_exclude'
                      configuration parameter.

                      It is an error to define both these configuration
                      parameters.

                      The old parameter 'biotypes' is equivalent to
                      'biotypes_include'.
        Return type : list[ensembl.api.core.Gene] (or empty list)
        Exceptions  : none
        Caller      : internal
        Status      : At Risk
                    : under development
        """
        if (cls._conf['biotypes_include'] or cls._conf['biotypes']) and cls._conf()['biotypes_exclude']:
            # logger->error("...")
            warnings.warn("You may not use both 'biotypes_include' and 'biotypes_exclude' in the configuration",
                            UserWarning)
        
        filtered = []
        opt_reverse = False
        if cls._conf['biotypes_include']:
            biotypes = cls._conf['biotypes_include']
        elif cls._conf['biotypes']:
            biotypes = cls._conf['biotypes']
        else:
            biotypes = cls._conf['biotypes_exclude']
            opt_reverse = True

        for gene in genes:
            if (not opt_reverse and gene.biotype in biotypes or
                opt_reverse and gene.biotype not in biotypes):
                filtered.append(gene)
        
        return filtered
    

    @classmethod
    def get_by_key(cls, cache_name: str, cache_type: str, cache_key: str) -> TinyFeature:
        """
        Arg[1]      : cache_name: str - a cache name (e.g. 'genes_by_id')
        Arg[2]      : cache_type: str - a cache type (e.g. f"source.{slice_name}")
        Arg[3]      : cache_key: str - key of this entry (e.g. a gene dbID)
        Example     : tf = Cache.get_by_key('genes_by_id',
                        'source.chromosome:NCBI36:X:1:1000000:1', '1234');
        Description : Gets a TinyFeature object from the cache.
        Return type : ensembl.id_mapping.TinyFeature
        Exceptions  : thrown on wrong or missing arguments
        Caller      : internal
        Status      : At Risk
                    : under development
        """
        if not cache_name:
            raise ValueError(f'You must provide a cache name (e.g. genes_by_id).')
        if not cache_type:
            raise ValueError(f'You must provide a cache type.')
        if not cache_key:
            raise ValueError(f'You must provide a cache key (e.g. a gene dbID).')
        
        if cache_type not in cls._instance.get('loaded'):
            cls.read_and_merge(cache_type)
        
        return cls._cache[cache_type][cache_name][cache_key]
    

    @classmethod
    def get_by_name(cls, cache_name: str, cache_type: str) -> dict[str, TinyFeature]:
        """
        Arg[1]      : cache_name: str - a cache name (e.g. 'genes_by_id')
        Arg[2]      : cache_type: str - a cache type (e.g. f"source.{slice_name}")
        Example     : tf = Cache.get_by_name('genes_by_id',
                        'source.chromosome:NCBI36:X:1:1000000:1');
        Description : Gets a dict[str,TinyFeature] object from the cache.
        Return type : dict[str,TinyFeature]
        Exceptions  : thrown on wrong or missing arguments
        Caller      : internal
        Status      : At Risk
                    : under development
        """
        if not cache_name:
            raise ValueError(f'You must provide a cache name (e.g. genes_by_id).')
        if not cache_type:
            raise ValueError(f'You must provide a cache type.')
        
        if cache_type not in cls._instance.get('loaded'):
            cls.read_and_merge(cache_type)
        
        if cls._cache[cache_type][cache_name]:
            return cls._cache[cache_type][cache_name]
        return {}
    

    @classmethod
    def get_count_by_name(cls, cache_name: str, cache_type: str) -> int:
        """
        Arg[1]      : cache_name: str - a cache name (e.g. 'genes_by_id')
        Arg[2]      : cache_type: str - a cache type (e.g. f"source.{slice_name}")
        Example     : tf = Cache.get_by_key('genes_by_id',
                        'source.chromosome:NCBI36:X:1:1000000:1', '1234');
        Description : get the count of TinyFeature objects from the cache.
        Return type : int
        Exceptions  : thrown on wrong or missing arguments
        Caller      : internal
        Status      : At Risk
                    : under development
        """
        if not cache_name:
            raise ValueError(f'You must provide a cache name (e.g. genes_by_id).')
        if not cache_type:
            raise ValueError(f'You must provide a cache type.')
        
        if cache_type not in cls._instance.get('loaded'):
            cls.read_and_merge(cache_type)
        
        return len(cls._cache[cache_type][cache_name])
        

    @classmethod
    def dump_path(cls):
        if not cls._dump_path:
            cls._dump_path = PurePath(cls._conf['basedir'], 'cache')
        return Cache._dump_path
    
    
    @classmethod
    def cache_file(cls, cache_type: str) -> Path:
        if not cache_type:
            raise ValueError(f'You must provide a cache type.')
        return Path(cls.dump_path(), f'{cache_type}.object_cache.ser')
    

    @classmethod
    def instance_file(cls) -> Path:
        return Path(cls.dump_path(), f'cache_instance.ser')
    

    @classmethod
    def write_instance_to_file(cls) -> int:
        filepath = cls.instance_file()
        filepath.parent.mkdir(parents=True, exist_ok=True)
        return cls._serialize_dump(cls._instance, filepath)
    

    @classmethod
    def write_to_file(cls, cache_type: str) -> int:
        if not cache_type:
            raise ValueError(f'You must provide a cache type.')
        
        if not cls._cache[cache_type]:
            # logger->error("...")
            warnings.warn(f"No features found in {cache_type}. Won't write cache file.",
                            UserWarning)
            return 0
        filepath = cls.cache_file(cache_type)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        return cls._serialize_dump(cls._cache[cache_type], filepath)
    

    @classmethod
    def write_all_to_file(cls, cache_type: str) -> int:
        if not cache_type:
            raise ValueError(f'You must provide a cache type.')

        size = 0
        size += cls.write_to_file(cache_type)
        size += cls.write_instance_to_file()

        return format_bytes(size)
    

    @classmethod
    def _serialize_dump(cls, cache_obj, file_path: Path) -> int:
        try:
            with open(file_path, 'wb') as fh:
                pickle.dump(cache_obj, fh, pickle.HIGHEST_PROTOCOL)
        except OSError:
            raise OSError(f'Unable to store {file_path}')
        
        return Path(file_path).stat().st_size
    

    @classmethod
    def read_from_file(cls, cache_type: str):
        if not cache_type:
            raise ValueError(f"You must provide a cache type.")

        cache_file = cls.cache_file(cache_type)

        if cache_file.exists():     
            #$self->logger->info("Reading cache from file...\n", 0, 'stamped');
            #$self->logger->info("Cache file $cache_file.\n", 1);
            try:
                with open(cache_file, 'rb') as fh:
                    # The protocol version used is detected automatically, so we do not
                    # have to specify it.
                    cls._cache[cache_type] = pickle.load(fh)
            except OSError:
                    raise OSError(f'Unable to retrieve cache {cache_type}')
      
            #$self->logger->info("Done.\n", 0, 'stamped');
        else:
            # $self->logger->warning("Cache file $cache_file not found or empty.\n");
            warnings.warn(f"Cache file {cache_file} not found or empty.",
                            UserWarning)

        return cls._cache[cache_type]
    

    @classmethod
    def read_and_merge(cls, dbtype: str) -> None:
        if not dbtype:
            raise ValueError(f'Db type (source|target) must be provided.')
        if dbtype not in ('source', 'target'):
            raise ValueError("Db type must be 'source' or 'target'.")

        # read cache from single or multiple files, depending on caching strategy
        cache_method = cls.cache_method()
        if cache_method == 'ALL':
            cls.read_from_file(f"{dbtype}.ALL")
        elif cache_method == 'BY_SEQ_REGION':
            for slice_name in cls.slice_names(dbtype):
                cls.read_from_file(f"{dbtype}.{slice_name}");
        else:
            raise ValueError(f"Unknown cache method: {cache_method}.")

        cls.merge(dbtype)

        # flag as being loaded
        cls._instance['loaded'].append(dbtype)

    
    @classmethod
    def merge(cls, dbtype: str) -> None:
        if not dbtype:
            raise ValueError(f'Db type (source|target) must be provided.')
        if dbtype not in ('source', 'target'):
            raise ValueError("Db type must be 'source' or 'target'.")

        for mytype in cls._cache.keys():
            if not dbtype in mytype:
                continue

            for name in cls._cache[dbtype].keys():
                for key in cls._cache[dbtype][name]:
                    if cls._cache[dbtype][name][key]:
                        # warning("Duplicate key in cache: $name|$dbtype|$key. Skipping.\n");
                        pass
                    else:
                        cls._cache[dbtype][name][key] = cls._cache[mytype][name][key]

                    del(cls._cache[mytype][name][key])
            
                del(cls._cache[mytype][name])
            
            del(cls._cache[mytype])

    
    @classmethod
    def slice_names(cls, session: Session, dbtype: str):
        if not dbtype:
            raise ValueError(f'You must provide a db type (source|target).')
        if dbtype not in ('source', 'target'):
            raise ValueError("Db type must be 'source' or 'target'.")
        if not session:
            raise ValueError(f'You must provide a session to connect to a core DB.')

        slice_names = []

        if cls._conf['chromosomes']:
        # if ( $self->conf->param('chromosomes') ) {
            # Fetch the specified chromosomes.
            for chr in cls._conf['chromosomes']:
                slice = SliceAdaptor.fetch_by_seq_region(session, chr, 'chromosome')
                if len(slice) > 1:
                    raise Exception()
                slice_names.append(slice[0].name())
        elif cls._conf['region']:
            # Fetch the slices on the specified regions.  Don't use
            # SliceAdaptor->fetch_by_name() since this will fail if assembly
            # versions are different for source and target db.
            (cs_name, cs_version, sr_name, start, end, strand) = ':'.split(cls._conf['region'])
            strand = Strand(strand)

            slice = SliceAdaptor.fetch_by_seq_region(session, 
                                                     sr_name,
                                                     cs_name,
                                                     cs_version, 
                                                     start, 
                                                     end)
            if len(slice) > 1:
                    raise Exception()
            slice_names.append(slice[0].name())
        else:
            # Fetch all slices that have genes on them.
            for srid in GeneAdaptor.list_seq_region_ids(session, 'species_id goes here'):
                slice = SliceAdaptor.fetch_by_seq_region_id(session, srid)
                # outdated - there are no more assembly exceptions
                # slices = SliceAdaptor.fetch_by_seq_region(session, slice.seq_region_name, slice.coord_system_name, cls._instance.get('hccsv'))
                slice_names.append(slice.name())

        return slice_names


from ensembl.database.dbconnection import DBConnection
def main():
    dbc1 = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
    dbc2 = DBConnection('mysql://ensro@mysql-ens-mirror-1.ebi.ac.uk:4240/homo_sapiens_core_109_38')
    # dbc1 = DBConnection('mysql://ensro@mysql-ens-mirror-1.ebi.ac.uk:4240/esox_lucius_core_109_4')
    # dbc2 = DBConnection('mysql://ensro@mysql-ens-mirror-1.ebi.ac.uk:4240/esox_lucius_core_108_4')
    with dbc1.session_scope() as tgt_session, dbc2.session_scope() as src_session:
        # aa = Cache.find_common_coord_systems(src_session, tgt_session)
        # print(f'Found {len(aa)} records')
        # for a in aa:
        #     print(f'\t{a}')
        # g = GeneAdaptor.fetch_by_stable_id(src_session, 'ENSG00000186092', load_transcripts=True,
        #                                    load_exons=True, load_translations=True)
        # genes = (g,)
        # num = Cache.build_cache_from_genes('source', genes, False)
        # print(num)
        slice_name = 'chromosome:GRCh38:13:32315086:32400268:1'
        slice_name = 'chromosome:GRCh38:20:1:64444167:1'
        # slice_name = 'chromosome:GRCh38:MT:1:16569'
        (num_genes, bytesize) = Cache.build_cache_by_slice(src_session, tgt_session, 'source', slice_name)
        print(f'Dumped {num_genes} genes in a file {bytesize} in size')
        
if __name__ == '__main__':
    main()