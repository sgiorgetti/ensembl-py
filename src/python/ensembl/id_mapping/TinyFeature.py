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

# =head1 CONTACT

#   Please email comments or questions to the public Ensembl
#   developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

#   Questions may also be sent to the Ensembl help desk at
#   <http://www.ensembl.org/Help/Contact>.

# =cut

# =head1 NAME

# Bio::EnsEMBL::IdMapping::TinyFeature - lightweight feature object

# =head1 SYNOPSIS

# This object isn't instantiated. See objects which inherit from it
# (TinyGene, TinyTranscript, etc.) for examples.

# =head1 DESCRIPTION

# This is the base class for the lightweight feature objects used by the
# stable Id maping application. For performance reasons, these objects
# are instantiated using a new_fast() method. The internal implementation
# is an arrayref (rather than the more common hashref), which optimises
# memory usage.

# There are no adaptors to fetch TinyFeatures from the database. You
# rather use the normal feature adaptors and then create the TinyFeatures
# from the heavy objects you get. The memory saving will therefore mainly
# take effect when serialising and reloading these objects.

# Also note that TinyFeatures don't have a slice attached to them - all
# location information (where required) is stored on the feature object
# directly.

# =head1 METHODS

#   new_fast
#   id
#   stable_id
#   version
#   created_date
#   modified_date
#   to_string

# =cut

__all__ = [ 'TinyFeature' ]

# import abc
# class TinyFeature(metaclass=abc.ABCMeta):
class TinyFeature():

    def __init__(self,
                 internal_id: int,
                 stable_id: str,
                 version: int,
                 created_date: int,
                 modified_date: int) -> None:
        """
        Description : Constructor.
        Return type : None
        Exceptions  : none
        Caller      : Bio.EnsEMBL.Id_Mapping.Cache
        Status      : At Risk
                    : under development
        """
        self._internal_id = internal_id
        self._stable_id = stable_id
        self._version = version
        self._created_date = created_date
        self._modified_date = modified_date

    def __repr__(self) -> str:
        """
        Example     : print(f"Created {feature}");
        Description : Prints a string representation of the feature for debug
                      purposes.
        Return type : String
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return f'{self.__class__.__name__}({self._internal_id}:{self._stable_id}.{self._version})'
    
    def __str__(self) -> str:
        return f"{self._internal_id}:{self._stable_id}.{self._version}"
    

    @classmethod
    def fast_from_tuple(cls, data: tuple) -> None:
        if len(data) != 5:
            raise ValueError(f"Argument tuple length must be 5: (id, stable_id, version, created_date, modified_date)")
        cls.__init__(data)

    @classmethod
    def fast_from_dict(cls, data: dict) -> None:
        version = data.get('version') if data.get('version') else 1
        cls.__init__(
            data.get('id'),
            data.get('stable_id'),
            version,
            data.get('created_date'),
            data.get('modified_date')
        )

    @property
    def internal_id(self) -> int:
        """
        Description : Getter for the feature's internal Id.
        Return type : Int
        Exceptions  : none
        Caller      : Bio::EnsEMBL::IdMapping::Cache
        Status      : At Risk
                    : under development
        """
        return self._internal_id

    @property
    def stable_id(self) -> str:
        """
        Description : Getter for the feature's stable Id.
        Return type : String
        Exceptions  : none
        Caller      : Bio::EnsEMBL::IdMapping::Cache
        Status      : At Risk
                    : under development
        """
        return self._stable_id
    
    @property
    def version(self) -> int:
        """
        Description : Getter for the feature's stable Id version.
        Return type : Int
        Exceptions  : none
        Caller      : Bio::EnsEMBL::IdMapping::Cache
        Status      : At Risk
                    : under development
        """
        return self._version

    @property
    def created_date(self) -> int:
        """
        Description : Getter for the feature's creation date.
        Return type : Int
        Exceptions  : none
        Caller      : Bio::EnsEMBL::IdMapping::Cache
        Status      : At Risk
                    : under development
        """
        return self._created_date
    
    @property
    def modified_date(self) -> int:
        """
        Description : Getter for the feature's modification date.
        Return type : Int
        Exceptions  : none
        Caller      : Bio::EnsEMBL::IdMapping::Cache
        Status      : At Risk
                    : under development
        """
        return self._modified_date

