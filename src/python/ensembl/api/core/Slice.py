from ensembl.api.core.Location import Location
from ensembl.api.core.Region import Region
from ensembl.api.core.Strand import Strand

from ensembl.database.dbconnection import DBConnection

__all__ = [ 'Slice' ]

class Slice():

    def __init__(self, region: Region = None, location: Location = None, strand: Strand = None) -> None:
        if not region:
            raise ValueError('Region object is required to instantiate a Slice')
        self._location = location
        self._region = region
        self._strand = strand

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._region._coord_system._version}:{self._region.name}:{self._location.start}-{self._location.end}:{self._strand.name})'

    @property
    def location(self) -> Location:
        return self._location

    @location.setter
    def location(self, value: Location) -> None:
        self._location = value

    @property
    def region(self) -> Region:
        return self._region

    @property
    def strand(self) -> Strand:
        return self._strand

    @strand.setter
    def strand(self, value: Strand) -> None:
        self._strand = value

    @property
    def name(self) -> str:
        # 'chromosome:NCBI36:X:1:1000000:-1'
        vals = (
            self._region.code.name.lower(),
            self._region.coord_system.version,
            self._region.name,
            self._location.start,
            self._location.end,
            self._strand.value
        )
        return ':'.join(vals)
    
    def get_all_genes(self):
        # Arg [1]    : (optional) string $logic_name
        #             The name of the analysis used to generate the genes to retrieve
        # Arg [2]    : (optional) string $dbtype
        #             The dbtype of genes to obtain.  This assumes that the db has
        #             been added to the DBAdaptor under this name (using the
        #             DBConnection::add_db_adaptor method).
        # Arg [3]    : (optional) boolean $load_transcripts
        #             If set to true, transcripts will be loaded immediately rather
        #             than being lazy-loaded on request.  This will result in a
        #             significant speed up if the Transcripts and Exons are going to
        #             be used (but a slow down if they are not).
        # Arg [4]    : (optional) string $source
        #             The source of the genes to retrieve.
        # Arg [5]    : (optional) string $biotype
        #             The biotype of the genes to retrieve.
        # Example    : @genes = @{$slice->get_all_Genes};
        # Description: Retrieves all genes that overlap this slice, including those on
        #             the reverse strand.
        # Returntype : listref of Bio::EnsEMBL::Genes
        # Exceptions : none
        # Caller     : none
        # Status     : Stable
        # my ($self, $logic_name, $dbtype, $load_transcripts, $source, $biotype) = @_;
        # if(my $adaptor = $self->_get_CoreAdaptor('Gene', $dbtype)) {
        #     return $adaptor->fetch_all_by_Slice( $self, $logic_name, $load_transcripts, $source, $biotype);
        # }
        # return [];
        with self._get_coreadaptor_conn('Gene').session_scope() as session:
            pass

    #     SliceAdaptor.fetch_by_name(session, self.name)
        

    def _get_coreadaptor_conn(self, obj_type: str, dbtype: str = 'core') -> DBConnection:
        """
        Arg  [1]    : Str object_type to retrieve an adaptor for
        Arg  [2]    : Str dbtype to search for the given adaptor in. Defaults to core
        Description : Searches for the specified adaptor in the Registry and returns it. Otherwise
                      it will return nothing if the adaptor was not found
        ReturnType  : Bio::EnsEMBL::DBSQL::BaseAdaptor derived instance (specific to core-like dbs)
        Exceptions  : missing object_type
        """
        if not obj_type:
            raise ValueError('Object type is a required parameter')
        if not dbtype:
            dbtype = 'core'
        return self._get_Adaptor(obj_type, dbtype)

    def _get_Adaptor(self, obj_type: str, dbtype: str, check_db: bool = False) -> DBConnection:
        """
        Arg  [1]    : Str object_type to retrieve an adaptor for
        Arg  [2]    : Str dbtype to search for the given adaptor in
        Arg  [3]    : Boolean Turn off the checking of Registry->get_db() for your 
                      adaptor.
        Description : Searches for the specified adaptor in the Registry and returns it. Otherwise
                      it will return nothing if the adaptor was not found. We consult the 
                      "special" adaptors held by Bio::EnsEMBL::Registry::get_db() method and then
                      fall back to the normal methods of finding an adaptor.
                      This method will warn when adaptors are missing but will never through an
                      exception. It is up to the calling code to decide how to handle the unavailablity
                      of an adaptor.
        ReturnType  : Bio::EnsEMBL::DBSQL::BaseAdaptor derrived instance. Otherwise it returns nothing
        Exceptions  : none
        """
        ######Example:   Gene      core       ??
        # my ($self, $object_type, $dbtype, $do_not_check_db) = @_;

        if not obj_type:
            raise ValueError('Object type is a required parameter')

        species = self._region.coord_system.species
        assembly_name = self._region._coord_system.version

        # #First we query for the DBAdaptor using get_db(). This is a deprecated method
        # #call but "special" adaptors can be registered via this method. We must
        # #consult here 1st to find the possible special adaptor
        # if(!$do_not_check_db && $dbtype) {
        #     my $db = $registry->get_db($local_db, $dbtype);
        #     if($db) {
        #     # If we got a return then use this DBAdaptor's species name, group and the given object type.
        #     # Special adaptors can have different species names
        #     $adaptor = $registry->get_adaptor($db->species(), $db->group(), $object_type);
        #     }
        #     else {
        #     #Otherwise just use the current species, dbtype and object type
        #     $adaptor = $registry->get_adaptor($species, $dbtype, $object_type);
        #     }
        # }
        # # Otherwise our query group is the one attached to the current adaptor
        # else {
        #     #If not set use the group attached to the local adaptor 
        #     $dbtype ||= $local_db->group();
        #     $adaptor = $registry->get_adaptor($species, $dbtype, $object_type);
        # }
        # return $adaptor if $adaptor;

        #A string like 'anonymous@ensembldb.ensembl.org:3306/homo_sapiens_core_109_38'
        #is returned by the (new)Registry

        return DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
        
        
