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

from time import time
from sqlalchemy.orm import Session
from ensembl.core.models import CoordSystem, SeqRegion

def fetch_species_id_by_seq_region_id(session: Session, seq_region_id: int) -> int:
    species_id = (
        session.query(CoordSystem.species_id)
        .join(SeqRegion.coord_system)
        .where(SeqRegion.seq_region_id == seq_region_id)
        .first()
    )
    if species_id:
        return species_id[0]
    return None

def fetch_species_id_by_seq_region_id(session: Session, seq_region_id: int) -> int:
    species_id = (
        session.query(CoordSystem.species_id)
        .join(SeqRegion.coord_system)
        .where(SeqRegion.seq_region_id == seq_region_id)
        .first()
    )
    if species_id:
        return species_id[0]
    return None

def timeme(func):
    # This function shows the execution time of 
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f'Function {func.__name__!r} executed in {(t2-t1):.4f}s')
        return result
    return wrap_func