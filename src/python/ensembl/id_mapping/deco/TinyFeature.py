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

from dataclasses import dataclass, field

__all__ = [ 'TinyFeature' ]

@dataclass
class TinyFeature():
    internal_id: int
    stable_id: str
    version: int
    created_date: int = field(repr=False)
    modified_date: int = field(repr=False)
    
    def __str__(self) -> str:
        return f"{self.__class__.__name__}({self.internal_id}:{self.stable_id}.{self.version})"
