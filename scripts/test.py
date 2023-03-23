from ensembl.database.dbconnection import DBConnection
from ensembl.api.dbsql import ExonAdaptor, SliceAdaptor, CoordSystemAdaptor, TranscriptAdaptor, AssemblyAdaptor, GeneAdaptor, BiotypeAdaptor
from ensembl.core.models import Exon as ExonORM
from ensembl.api.core import Exon, SplicedExon, Transcript, Gene, Assembly

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from ensembl.api.fileio.GFF3Adaptor import GFF3Writer, write as gff_write
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def testDumpGFFGene(gene: Gene, assembly: Assembly):
    out_file = "your_file.gff"
    seq_region_name = f"{gene.get_slice().seq_region_name}"
    ensembl_dirs = {
        "##sequence-region  ": f"{seq_region_name} {gene.get_slice().start} {gene.get_slice().end}",
        "#!genome-build": f"{assembly.accessioning_body} {assembly.id}",
        "#!genome-version": f"{assembly.name}",
        "#!genome-date": f"{assembly.assembly_date}",
        "#!genome-build-accession": f"{assembly.accession_id}",
        "#!genebuild-last-updated": f"{assembly.gb_last_geneset_update}",
    }
    seq = Seq("UFFA")
    rec = SeqRecord(seq, f"{seq_region_name}")
# 1       havana  mRNA    65419   71585   .       +       .       ID=transcript:ENST00000641515;Parent=gene:ENSG00000186092;Name=OR4F5-201;biotype=protein_coding;tag=basic,Ensembl_canonical,MANE_Select;transcript_id=ENST00000641515;version=2
# 1       havana  exon    65419   65433   .       +       .       Parent=transcript:ENST00000641515;Name=ENSE00003812156;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003812156;rank=1;version=1    rec = SeqRecord(seq, f"{seq_region_name}", "name", "descr", ["xref1", "xref2"])
    qualifiers = gene.gff3_qualifiers()
    ### WARNING!!! FeatureLocation assumes 0-based start!!!
    top_feature = SeqFeature(
        FeatureLocation(gene.start-1, gene.end),
        type=f"{gene.type}",
        strand=gene.strand.value,
        qualifiers=qualifiers
    )

    sub_feats = []
    for transcript in gene.get_transcripts():
        sub_qualifiers = transcript.gff3_qualifiers()
        sub_feats.append(
            SeqFeature(FeatureLocation(transcript.start-1, transcript.end),
                       type=f"{transcript.biotype.so_term}",
                       strand=transcript.strand.value,
                       qualifiers=sub_qualifiers
            )
        )
        for exon in transcript.get_exons():
            sub_qualifiers = exon.gff3_qualifiers()
            sub_feats.append(
                SeqFeature(FeatureLocation(exon.start-1, exon.end),
                        type=f"{exon.type}",
                        strand=exon.strand.value,
                        qualifiers=sub_qualifiers
                )
            )
    
    top_feature.sub_features = sub_feats
    rec.features = [top_feature]

    with open(out_file, "w") as out_handle:
        gff_write([rec], out_handle, custom_dirs=ensembl_dirs)


def testDumpGFFTranscript(transcript: Transcript, assembly: Assembly):
    out_file = "your_file.gff"
    seq_region_name = f"{transcript.get_slice().seq_region_name}"
    ensembl_dirs = {
        "##sequence-region  ": f"{seq_region_name} {transcript.get_slice().start} {transcript.get_slice().end}",
        "#!genome-build": f"{assembly.accessioning_body} {assembly.id}",
        "#!genome-version": f"{assembly.name}",
        "#!genome-date": f"{assembly.assembly_date}",
        "#!genome-build-accession": f"{assembly.accession_id}",
        "#!genebuild-last-updated": f"{assembly.gb_last_geneset_update}",
    }
    seq = Seq("UFFA")
    rec = SeqRecord(seq, f"{seq_region_name}")
# 1       havana  mRNA    65419   71585   .       +       .       ID=transcript:ENST00000641515;Parent=gene:ENSG00000186092;Name=OR4F5-201;biotype=protein_coding;tag=basic,Ensembl_canonical,MANE_Select;transcript_id=ENST00000641515;version=2
# 1       havana  exon    65419   65433   .       +       .       Parent=transcript:ENST00000641515;Name=ENSE00003812156;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003812156;rank=1;version=1    rec = SeqRecord(seq, f"{seq_region_name}", "name", "descr", ["xref1", "xref2"])
    qualifiers = transcript.gff3_qualifiers()
    ### WARNING!!! FeatureLocation assumes 0-based start!!!
    top_feature = SeqFeature(
        FeatureLocation(transcript.start-1, transcript.end),
        type=f"{transcript.biotype.so_term}",
        strand=transcript.strand.value,
        qualifiers=qualifiers
    )

    sub_feats = []
    for exon in transcript.get_exons():
        sub_qualifiers = exon.gff3_qualifiers()
        sub_feats.append(
            SeqFeature(FeatureLocation(exon.start-1, exon.end),
                       type=f"{exon.type}",
                       strand=exon.strand.value,
                       qualifiers=sub_qualifiers
            )
        )
    top_feature.sub_features = sub_feats
    rec.features = [top_feature]

    with open(out_file, "w") as out_handle:
        gff_write([rec], out_handle, custom_dirs=ensembl_dirs)


def main2():
    # dbc = DBConnection('mysql://ensro@mysql-ens-core-prod-1.ebi.ac.uk:4524/sgiorgetti_homo_sapiens_core_109_38')
    # dbc = DBConnection('mysql://mysql-ens-sta-1-b.ebi.ac.uk:4685/homo_sapiens_core_109_38')
    dbc = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
    # dbc = DBConnection('mysql://ensro@mysql-ens-sta-3.ebi.ac.uk:4160/protists_amoebozoa1_collection_core_57_110_1')
    with dbc.session_scope() as session:
        sls = SliceAdaptor.fetch_by_seq_region(session, '13', 'chromosome', 'GRCh38', 32315086, 32400268, 1)
        for sl in sls:
            print(sl)
            genes = GeneAdaptor.fetch_all_by_Slice(session, sl, load_transcripts=True)
            for g in genes:
                print(f"\t{g}")
                for t in g.get_transcripts():
                    print(f"\t\t{t}")

def main3():
    dbc = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
    with dbc.session_scope() as session:
        g = GeneAdaptor.fetch_by_stable_id(session, 'ENSG00000186092', True)
        print(f'gene: {g.stable_id} - {g.internal_id}')
        print(f"gene({g.internal_id}) has {len(g.get_transcripts())} transcripts")
        for t in g.get_transcripts():
            print(f"\ttranscript {t.stable_id} - isCanonical {t.is_canonical()}")

def main4():
    dbc = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
    # dbc = DBConnection('mysql://ensro@mysql-ens-mirror-1.ebi.ac.uk:4240/homo_sapiens_core_109_38')
    
    with dbc.session_scope() as session:
        # ss = SliceAdaptor._fetch_all_seq_regions_by_coord_system_id(session, 4)
        # ss = SliceAdaptor.fetch_by_seq_region(session, 1, 'chromosome')
        # ss = SliceAdaptor.fetch_by_seq_region_id(session, 132907)
        # ss = SliceAdaptor.fetch_by_name(session, 'chromosome:GRCh38:MT:1:16569')
        # ss = SliceAdaptor.fetch_all(session, 'chromosome', 'GRCh38')
        # ss = BiotypeAdaptor.fetch_by_name_object_type(session, 'protein_coding', 'gene')
        ss = BiotypeAdaptor.fetch_all_by_object_type(session, 'gene')
        print(ss[0])
        ss = BiotypeAdaptor.fetch_by_name_object_type(session, 'protein_coding', 'transcript')
        print(ss)
        ss = BiotypeAdaptor.fetch_all_by_object_type(session, 'gene')
        print(ss[1])

def main():
    # dbc = DBConnection('mysql://ensro@mysql-ens-core-prod-1.ebi.ac.uk:4524/sgiorgetti_homo_sapiens_core_109_38')
    # dbc = DBConnection('mysql://LOCALHOST:3306/sgiorgetti_havana_human')
    dbc = DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
    with dbc.session_scope() as session:
    #     print(ea.fetch_by_stable_id(session, 'ENSE00002089356'))
    # print(ExonAdaptor.fetch_by_stable_id(dbc, 'ENSE00001544499'))
    # t = TranscriptAdaptor.fetch_by_stable_id_version(dbc, 'ENST00000389680', 2)

        g = GeneAdaptor.fetch_by_stable_id(session, 'ENSG00000186092')
        # t = TranscriptAdaptor.fetch_by_stable_id(session, 'ENST00000641515')
        tr = TranscriptAdaptor.fetch_all_by_gene(session, g)
        for t in tr:
            print(t)
            exons = ExonAdaptor.fetch_all_by_Transcript(session, t)
            t.set_exons(exons)
            assembly = AssemblyAdaptor.fetch_by_coord_system_id(session, t.get_slice().coord_system.internal_id)
        g.set_transcripts(tr)
        testDumpGFFGene(g, assembly)



if __name__ == "__main__":
    main()
