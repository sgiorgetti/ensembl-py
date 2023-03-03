from ensembl.database.dbconnection import DBConnection
from ensembl.api.dbsql import ExonAdaptor, SliceAdaptor, CoordSystemAdaptor, TranscriptAdaptor
from ensembl.core.models import Exon as ExonORM
from ensembl.api.core import Exon, SplicedExon, Transcript

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def testDumpGFFTranscript(transcript: Transcript):
    out_file = "your_file.gff"
    qualifiers = {
        "source": "TBD(havana)",
        "rank": "TBD(1)",
        "Parent": "TBD(transcript:ENST00000389680)",
        "constitutive": "TBD(1)",
        "rank": "TBD(1)",
    }
    # ID=transcript:ENST00000641515;Parent=gene:ENSG00000186092;Name=OR4F5-201;biotype=protein_coding;tag=basic,Ensembl_canonical,MANE_Select;transcript_id=ENST00000641515;version=2
    # qualifiers["Name"] = e.unversioned_stable_id
    # qualifiers["exon_id"] = e.unversioned_stable_id
    # qualifiers["version"] = e.version
    # qualifiers["ensembl_phase"] = e.phase
    # qualifiers["ensembl_end_phase"] = e.end_phase
    top_feature = SeqFeature(
        FeatureLocation(transcript.get_slice().location.start, transcript.get_slice().location.end),
        type=transcript.type.lower(),
        strand=transcript.get_slice().strand.value,
        qualifiers=qualifiers
    )
    seq = Seq("")
    rec = SeqRecord(seq, "MT", "name", "descr", ["xref1", "xref2"])
    rec.features = [top_feature]

    with open(out_file, "w") as out_handle:
        GFF.write([rec], out_handle)

def main3():
    out_file = "your_file.gff"

    dbc = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
    e = ExonAdaptor.fetch_by_stable_id_version(dbc, 'ENSE00001544499', 2)
    print(e)

    qualifiers = {
        "source": "TBD(havana)",
        "rank": "TBD(1)",
        "Parent": "TBD(transcript:ENST00000389680)",
        "constitutive": "TBD(1)",
        "rank": "TBD(1)",
    }
    qualifiers["Name"] = e.unversioned_stable_id
    qualifiers["exon_id"] = e.unversioned_stable_id
    qualifiers["version"] = e.version
    qualifiers["ensembl_phase"] = e.phase
    qualifiers["ensembl_end_phase"] = e.end_phase



    top_feature = SeqFeature(
        FeatureLocation(e.get_slice().location.start, e.get_slice().location.end),
        type=e.type.lower(),
        strand=e.get_slice().strand.value,
        qualifiers=qualifiers
    )
    seq = Seq("")
    # # rec = SeqRecord(seq, "ID1")
    # rec = SeqRecord(seq, "ID1", "name", "descr", ["xref1", "xref2"])
    # qualifiers = {
    #     "source": "prediction",
    #     "score": 10.0,
    #     "other": ["Some", "annotations"],
    #     "ID": "gene1",
    # }
    # sub_qualifiers = {"source": "prediction"}
    # top_feature = SeqFeature(
    #     FeatureLocation(0, 20), type="gene", strand=1, qualifiers=qualifiers
    # )
    # top_feature.sub_features = [
    #     SeqFeature(FeatureLocation(0, 5), type="exon", strand=1, qualifiers=sub_qualifiers),
    #     SeqFeature(
    #         FeatureLocation(15, 20), type="exon", strand=1, qualifiers=sub_qualifiers
    #     ),
    # ]
    rec = SeqRecord(seq, "MT", "name", "descr", ["xref1", "xref2"])
    rec.features = [top_feature]

    with open(out_file, "w") as out_handle:
        GFF.write([rec], out_handle)

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
    # t = TranscriptAdaptor.fetch_by_stable_id_version(dbc, 'ENST00000389680', 2)
    with dbc.session_scope() as session:
        t = TranscriptAdaptor.fetch_by_stable_id(session, 'ENST00000641515')
        print(t)
        exons = ExonAdaptor.fetch_all_by_Transcript(session, t)
        # for e in exons:
        #     print(f"{e.index} - {e}")
        t.set_exons(exons)
    # existing_databases = mysql_engine.execute("SHOW DATABASES;")
    # existing_databases2 = [d[0] for d in existing_databases]
    # ENSE00002089356

if __name__ == "__main__":
    main()
