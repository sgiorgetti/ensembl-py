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

__all__ = [ 'Cache' ]

class Cache():
    """
    ensembl.id_mapping.Cache - a cache to hold data objects used by the 
    id_mapping application
    """

    @staticmethod
    def build_cache_by_slice(session: Session, target, dbtype, slice_name):
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
        # $self->cache_method('BY_SEQ_REGION');
        slice = SliceAdaptor.fetch_by_name(session, slice_name)
        genes = GeneAdaptor.fetch_all_by_Slice(session, slice, load_transcripts=True)

        # find common coord_system
        # $common_cs_found = $self->find_common_coord_systems;
        pass
    
  

# sub build_cache_by_slice {
#   my $self       = shift;
#   my $dbtype     = shift;
#   my $slice_name = shift;

#   # set cache method (required for loading cache later)
#   $self->cache_method('BY_SEQ_REGION');

#   my $dba = $self->get_DBAdaptor($dbtype);
#   my $sa  = $dba->get_SliceAdaptor;

#   my $slice = $sa->fetch_by_name($slice_name);
#   unless ($slice) {
#     throw("Could not retrieve slice $slice_name.");
#   }

#   my $genes = $slice->get_all_Genes( undef, undef, 1 );

#   # find common coord_system
#   my $common_cs_found = $self->find_common_coord_systems;

#   # find out whether native coord_system is a common coord_system.
#   # if so, you don't need to project.
#   # also don't project if no common coord_system present
#   my $need_project = 1;

#   my $csid = join( ':',
#                    $slice->coord_system_name,
#                    $slice->coord_system->version );

#   if ( $self->is_common_cs($csid) or !$self->highest_common_cs ) {
#     $need_project = 0;
#   }

#   # build cache
#   my $type = "$dbtype.$slice_name";
#   my $num_genes =
#     $self->build_cache_from_genes( $type, $genes, $need_project );
#   undef $genes;

#   # write cache to file, then flush cache to reclaim memory
#   my $size = $self->write_all_to_file($type);

#   return $num_genes, $size;
# } ## end sub build_cache_by_slice
