from ensembl.database.dbconnection import DBConnection
from sqlalchemy import select, and_
from sqlalchemy.orm import Session, Bundle, aliased
from ensembl.core.models import Gene, SeqRegion, StableIdEvent
import datatable as dt

current_genes = []
history = []

def fetch_current_genes_version(session: Session):
    stmt = (
        select(Gene.stable_id, Gene.version)
        .join(SeqRegion, SeqRegion.seq_region_id == Gene.seq_region_id)
        .where(
            and_(
                Gene.is_current == 1,
                SeqRegion.coord_system_id == 4
            )
        )
    )
    global current_genes
    current_genes = session.execute(stmt).all()

def fetch_old_genes(session: Session):
    stmt = (
            select(
                StableIdEvent
            )
        )
    global history
    history = session.scalars(stmt).all()
    for h in history:
        print(h)

# dbc = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
dbc = DBConnection('mysql://ensro@mysql-ens-mirror-1.ebi.ac.uk:4240/homo_sapiens_core_109_38')
with dbc.session_scope() as session:
    fetch_current_genes_version(session)
    


