from ensembl.database.dbconnection import DBConnection
from ensembl.api.dbsql import ExonAdaptor, SliceAdaptor, CoordSystemAdaptor
# from ensembl.core.models import Exon as ExonORM
# from ensembl.api.core.Exon import Exon

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session


def main2():
    # dbc = DBConnection('mysql://ensro@mysql-ens-core-prod-1.ebi.ac.uk:4524/sgiorgetti_homo_sapiens_core_109_38')
    # dbc = DBConnection('mysql://mysql-ens-sta-1-b.ebi.ac.uk:4685/homo_sapiens_core_109_38')
    dbc = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
    # dbc = DBConnection('mysql://ensro@mysql-ens-sta-3.ebi.ac.uk:4160/protists_amoebozoa1_collection_core_57_110_1')
    # with dbc.session_scope() as session:
    aaa = SliceAdaptor.fetch_by_seq_region(dbc, '19','chromosome', None, 1, 1000, 1)
    for a in aaa:
        print(a)

def main():
    # dbc = DBConnection('mysql://ensro@mysql-ens-core-prod-1.ebi.ac.uk:4524/sgiorgetti_homo_sapiens_core_109_38')
    # dbc = DBConnection('mysql://LOCALHOST:3306/sgiorgetti_havana_human')
    dbc = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
    # with dbc.session_scope() as session:
    #     print(ea.fetch_by_stable_id(session, 'ENSE00002089356'))
    # print(ExonAdaptor.fetch_by_stable_id(dbc, 'ENSE00001544499'))
    e = ExonAdaptor.fetch_by_stable_id_version(dbc, 'ENSE00001544499', 2)
    print(e)
    # existing_databases = mysql_engine.execute("SHOW DATABASES;")
    # existing_databases2 = [d[0] for d in existing_databases]
    # ENSE00002089356

if __name__ == "__main__":
    main()
