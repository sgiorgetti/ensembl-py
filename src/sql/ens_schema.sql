--
-- EnsEMBL core schema dump for duckdb
--

--
-- Name: ensembl_core_schema; Type: SCHEMA; Schema: -;
--

CREATE SCHEMA ensembl_core_schema;


--
-- Name: alt_allele_attrib_attrib; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE alt_allele_attrib_attrib AS ENUM (
    'IS_REPRESENTATIVE',
    'IS_MOST_COMMON_ALLELE',
    'IN_CORRECTED_ASSEMBLY',
    'HAS_CODING_POTENTIAL',
    'IN_ARTIFICIALLY_DUPLICATED_ASSEMBLY',
    'IN_SYNTENIC_REGION',
    'HAS_SAME_UNDERLYING_DNA_SEQUENCE',
    'IN_BROKEN_ASSEMBLY_REGION',
    'IS_VALID_ALTERNATE',
    'SAME_AS_REPRESENTATIVE',
    'SAME_AS_ANOTHER_ALLELE',
    'MANUALLY_ASSIGNED',
    'AUTOMATICALLY_ASSIGNED',
    'IS_PAR'
);


--
-- Name: biotype_biotype_group; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE biotype_biotype_group AS ENUM (
    'coding',
    'pseudogene',
    'snoncoding',
    'lnoncoding',
    'mnoncoding',
    'LRG',
    'undefined',
    'no_group'
);


--
-- Name: biotype_object_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE biotype_object_type AS ENUM (
    'gene',
    'transcript'
);


--
-- Name: data_file_file_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE data_file_file_type AS ENUM (
    'BAM',
    'BAMCOV',
    'BIGBED',
    'BIGWIG',
    'VCF'
);


--
-- Name: density_type_value_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE density_type_value_type AS ENUM (
    'sum',
    'ratio'
);


--
-- Name: ditag_feature_ditag_side; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE ditag_feature_ditag_side AS ENUM (
    'F',
    'L',
    'R'
);


--
-- Name: dna_align_feature_align_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE dna_align_feature_align_type AS ENUM (
    'ensembl',
    'cigar',
    'vulgar',
    'mdtag'
);


--
-- Name: external_db_status; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE external_db_status AS ENUM (
    'KNOWNXREF',
    'KNOWN',
    'XREF',
    'PRED',
    'ORTH',
    'PSEUDO'
);


--
-- Name: external_db_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE external_db_type AS ENUM (
    'ARRAY',
    'ALT_TRANS',
    'ALT_GENE',
    'MISC',
    'LIT',
    'PRIMARY_DB_SYNONYM',
    'ENSEMBL'
);


--
-- Name: intron_supporting_evidence_score_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE intron_supporting_evidence_score_type AS ENUM (
    'NONE',
    'DEPTH'
);


--
-- Name: marker_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE marker_type AS ENUM (
    'est',
    'microsatellite'
);


--
-- Name: object_xref_ensembl_object_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE object_xref_ensembl_object_type AS ENUM (
    'RawContig',
    'Transcript',
    'Gene',
    'Translation',
    'Operon',
    'OperonTranscript',
    'Marker',
    'RNAProduct'
);


--
-- Name: protein_align_feature_align_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE protein_align_feature_align_type AS ENUM (
    'ensembl',
    'cigar',
    'vulgar',
    'mdtag'
);


--
-- Name: protein_feature_align_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE protein_feature_align_type AS ENUM (
    'ensembl',
    'cigar',
    'cigarplus',
    'vulgar',
    'mdtag'
);


--
-- Name: stable_id_event_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE stable_id_event_type AS ENUM (
    'gene',
    'transcript',
    'translation',
    'rnaproduct'
);


--
-- Name: supporting_feature_feature_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE supporting_feature_feature_type AS ENUM (
    'dna_align_feature',
    'protein_align_feature'
);


--
-- Name: transcript_supporting_feature_feature_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE transcript_supporting_feature_feature_type AS ENUM (
    'dna_align_feature',
    'protein_align_feature'
);


--
-- Name: unmapped_object_ensembl_object_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE unmapped_object_ensembl_object_type AS ENUM (
    'RawContig',
    'Transcript',
    'Gene',
    'Translation'
);


--
-- Name: unmapped_object_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE unmapped_object_type AS ENUM (
    'xref',
    'cDNA',
    'Marker'
);


--
-- Name: xref_info_type; Type: TYPE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TYPE xref_info_type AS ENUM (
    'NONE',
    'PROJECTION',
    'MISC',
    'DEPENDENT',
    'DIRECT',
    'SEQUENCE_MATCH',
    'INFERRED_PAIR',
    'PROBE',
    'UNMAPPED',
    'COORDINATE_OVERLAP',
    'CHECKSUM'
);


--
-- Name: alt_allele; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.alt_allele (
    alt_allele_id UBIGINT PRIMARY KEY,
    alt_allele_group_id UBIGINT NOT NULL,
    gene_id UBIGINT NOT NULL
);


--
-- Name: alt_allele_alt_allele_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.alt_allele_alt_allele_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: alt_allele_attrib; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.alt_allele_attrib (
    alt_allele_id UBIGINT,
    attrib alt_allele_attrib_attrib,
    FOREIGN KEY (alt_allele_id) REFERENCES ensembl_core_schema.alt_allele (alt_allele_id)
);


--
-- Name: alt_allele_group_alt_allele_group_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.alt_allele_group_alt_allele_group_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: analysis; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.analysis (
    analysis_id UINTEGER PRIMARY KEY,
    created timestamp with time zone,
    logic_name VARCHAR(128) NOT NULL,
    db VARCHAR(120),
    db_version VARCHAR(40),
    db_file VARCHAR(120),
    program VARCHAR(80),
    program_version VARCHAR(40),
    program_file VARCHAR(80),
    parameters VARCHAR,
    module VARCHAR(80),
    module_version VARCHAR(40),
    gff_source VARCHAR(40),
    gff_feature VARCHAR(40)
);


--
-- Name: analysis_analysis_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.analysis_analysis_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: analysis_description; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.analysis_description (
    analysis_id UINTEGER PRIMARY KEY,
    description VARCHAR,
    display_label VARCHAR(255) NOT NULL,
    displayable boolean DEFAULT true NOT NULL,
    web_data VARCHAR,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id)
);


--
-- Name: assembly; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.assembly (
    asm_seq_region_id UBIGINT NOT NULL,
    cmp_seq_region_id UBIGINT NOT NULL,
    asm_start UBIGINT NOT NULL,
    asm_end UBIGINT NOT NULL,
    cmp_start UBIGINT NOT NULL,
    cmp_end UBIGINT NOT NULL,
    ori TINYINT NOT NULL
);



--
-- Name: associated_group; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.associated_group (
    associated_group_id UBIGINT PRIMARY KEY,
    description VARCHAR(128)
);


--
-- Name: associated_group_associated_group_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.associated_group_associated_group_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: associated_xref; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.associated_xref (
    associated_xref_id UBIGINT PRIMARY KEY,
    object_xref_id UBIGINT DEFAULT '0' NOT NULL,
    xref_id UBIGINT DEFAULT '0' NOT NULL,
    source_xref_id UBIGINT,
    condition_type VARCHAR(128),
    associated_group_id UBIGINT,
    rank UBIGINT DEFAULT '0'
);


--
-- Name: associated_xref_associated_xref_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.associated_xref_associated_xref_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: attrib_type; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.attrib_type (
    attrib_type_id UINTEGER PRIMARY KEY,
    code VARCHAR(20) DEFAULT '' NOT NULL,
    name VARCHAR(255) DEFAULT '' NOT NULL,
    description VARCHAR
);


--
-- Name: attrib_type_attrib_type_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.attrib_type_attrib_type_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: biotype; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.biotype (
    biotype_id UBIGINT PRIMARY KEY,
    name VARCHAR(64) NOT NULL,
    object_type biotype_object_type DEFAULT 'gene' NOT NULL,
    db_type VARCHAR DEFAULT 'core' NOT NULL,
    attrib_type_id UINTEGER,
    description VARCHAR,
    biotype_group biotype_biotype_group,
    so_acc VARCHAR(64),
    so_term VARCHAR(1023)
);


--
-- Name: biotype_biotype_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.biotype_biotype_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: coord_system; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.coord_system (
    coord_system_id UINTEGER PRIMARY KEY,
    species_id UINTEGER DEFAULT '1' NOT NULL,
    name VARCHAR(40) NOT NULL,
    version VARCHAR(255),
    rank UINTEGER NOT NULL,
    attrib VARCHAR,
);


--
-- Name: coord_system_coord_system_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.coord_system_coord_system_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: data_file; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.data_file (
    data_file_id UBIGINT PRIMARY KEY,
    coord_system_id UINTEGER NOT NULL,
    analysis_id UINTEGER NOT NULL,
    name VARCHAR(100) NOT NULL,
    version_lock boolean DEFAULT false NOT NULL,
    absolute boolean DEFAULT false NOT NULL,
    url VARCHAR,
    file_type data_file_file_type,
    FOREIGN KEY (coord_system_id) REFERENCES ensembl_core_schema.coord_system (coord_system_id),
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id)
);


--
-- Name: data_file_data_file_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.data_file_data_file_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: seq_region; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.seq_region (
    seq_region_id UBIGINT PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    coord_system_id UINTEGER NOT NULL,
    length UBIGINT NOT NULL
);


--
-- Name: density_type; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.density_type (
    density_type_id UBIGINT PRIMARY KEY,
    analysis_id UINTEGER NOT NULL,
    block_size UINTEGER NOT NULL,
    region_features UINTEGER NOT NULL,
    value_type density_type_value_type NOT NULL,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id)
);


--
-- Name: density_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.density_feature (
    density_feature_id UBIGINT PRIMARY KEY,
    density_type_id UBIGINT NOT NULL,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    density_value DOUBLE NOT NULL,
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id),
    FOREIGN KEY (density_type_id) REFERENCES ensembl_core_schema.density_type (density_type_id)
);


--
-- Name: density_feature_density_feature_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.density_feature_density_feature_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: density_type_density_type_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.density_type_density_type_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: dependent_xref; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.dependent_xref (
    object_xref_id UBIGINT NOT NULL,
    master_xref_id UBIGINT NOT NULL,
    dependent_xref_id UBIGINT NOT NULL
);


--
-- Name: ditag; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.ditag (
    ditag_id UBIGINT PRIMARY KEY,
    name VARCHAR(30) NOT NULL,
    type VARCHAR(30) NOT NULL,
    tag_count UINTEGER DEFAULT 1 NOT NULL,
    sequence VARCHAR NOT NULL
);


--
-- Name: ditag_ditag_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.ditag_ditag_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: ditag_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.ditag_feature (
    ditag_feature_id UBIGINT PRIMARY KEY,
    ditag_id UBIGINT DEFAULT '0' NOT NULL,
    ditag_pair_id UBIGINT DEFAULT '0' NOT NULL,
    seq_region_id UBIGINT DEFAULT '0' NOT NULL,
    seq_region_start UBIGINT DEFAULT '0' NOT NULL,
    seq_region_end UBIGINT DEFAULT '0' NOT NULL,
    seq_region_strand boolean DEFAULT false NOT NULL,
    analysis_id UINTEGER DEFAULT 0 NOT NULL,
    hit_start UBIGINT DEFAULT '0' NOT NULL,
    hit_end UBIGINT DEFAULT '0' NOT NULL,
    hit_strand boolean DEFAULT false NOT NULL,
    cigar_line VARCHAR NOT NULL,
    ditag_side ditag_feature_ditag_side NOT NULL,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id),
    FOREIGN KEY (ditag_id) REFERENCES ensembl_core_schema.ditag (ditag_id)
);


--
-- Name: ditag_feature_ditag_feature_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.ditag_feature_ditag_feature_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: dna; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.dna (
    seq_region_id UBIGINT PRIMARY KEY,
    sequence VARCHAR NOT NULL
);


--
-- Name: dna_align_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.dna_align_feature (
    dna_align_feature_id UBIGINT PRIMARY KEY,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand boolean NOT NULL,
    hit_start UINTEGER NOT NULL,
    hit_end UINTEGER NOT NULL,
    hit_strand boolean NOT NULL,
    hit_name VARCHAR(40) NOT NULL,
    analysis_id UINTEGER NOT NULL,
    score DOUBLE,
    evalue DOUBLE,
    perc_ident DOUBLE,
    cigar_line VARCHAR,
    external_db_id UINTEGER,
    hcoverage DOUBLE,
    align_type dna_align_feature_align_type DEFAULT 'ensembl',
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: dna_align_feature_attrib; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.dna_align_feature_attrib (
    dna_align_feature_id UBIGINT NOT NULL,
    attrib_type_id UINTEGER NOT NULL,
    value VARCHAR NOT NULL,
    FOREIGN KEY (dna_align_feature_id) REFERENCES ensembl_core_schema.dna_align_feature (dna_align_feature_id),
    FOREIGN KEY (attrib_type_id) REFERENCES ensembl_core_schema.attrib_type (attrib_type_id)
);


--
-- Name: dna_align_feature_dna_align_feature_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.dna_align_feature_dna_align_feature_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: exon; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.exon (
    exon_id UBIGINT PRIMARY KEY,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand TINYINT NOT NULL,
    phase TINYINT NOT NULL,
    end_phase TINYINT NOT NULL,
    is_current boolean DEFAULT true NOT NULL,
    is_constitutive boolean DEFAULT false NOT NULL,
    stable_id VARCHAR(128),
    version UINTEGER,
    created_date timestamp with time zone,
    modified_date timestamp with time zone,
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: exon_exon_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.exon_exon_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: gene; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.gene (
    gene_id UBIGINT PRIMARY KEY,
    biotype VARCHAR(40) NOT NULL,
    analysis_id UINTEGER NOT NULL,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand TINYINT NOT NULL,
    display_xref_id UBIGINT,
    source VARCHAR(40) NOT NULL,
    description VARCHAR,
    is_current boolean DEFAULT true NOT NULL,
    canonical_transcript_id UBIGINT NOT NULL,
    stable_id VARCHAR(128),
    version UINTEGER,
    created_date timestamp with time zone,
    modified_date timestamp with time zone,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: transcript; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.transcript (
    transcript_id UBIGINT PRIMARY KEY,
    gene_id UBIGINT,
    analysis_id UINTEGER NOT NULL,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand TINYINT NOT NULL,
    display_xref_id UBIGINT,
    source VARCHAR(40) DEFAULT 'ensembl' NOT NULL,
    biotype VARCHAR(40) NOT NULL,
    description VARCHAR,
    is_current boolean DEFAULT true NOT NULL,
    canonical_translation_id UBIGINT,
    stable_id VARCHAR(128),
    version UINTEGER,
    created_date timestamp with time zone,
    modified_date timestamp with time zone,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id),
    FOREIGN KEY (gene_id) REFERENCES ensembl_core_schema.gene (gene_id)
);


--
-- Name: exon_transcript; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.exon_transcript (
    exon_id UBIGINT NOT NULL,
    transcript_id UBIGINT NOT NULL,
    rank UINTEGER NOT NULL,
    PRIMARY KEY(exon_id, transcript_id),
    FOREIGN KEY (exon_id) REFERENCES ensembl_core_schema.exon (exon_id),
    FOREIGN KEY (transcript_id) REFERENCES ensembl_core_schema.transcript (transcript_id)
);


--
-- Name: external_db; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.external_db (
    external_db_id UINTEGER PRIMARY KEY,
    db_name VARCHAR(100) NOT NULL,
    db_release VARCHAR(255),
    status external_db_status NOT NULL,
    priority UINTEGER NOT NULL,
    db_display_name VARCHAR(255),
    type external_db_type NOT NULL,
    secondary_db_name VARCHAR(255),
    secondary_db_table VARCHAR(255),
    description VARCHAR
);


--
-- Name: external_db_external_db_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.external_db_external_db_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: external_synonym; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.external_synonym (
    xref_id UBIGINT NOT NULL,
    synonym VARCHAR(100) NOT NULL
);


--
-- Name: gene_archive; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.gene_archive (
    gene_stable_id VARCHAR(128) NOT NULL,
    gene_version USMALLINT DEFAULT '1' NOT NULL,
    transcript_stable_id VARCHAR(128) NOT NULL,
    transcript_version USMALLINT DEFAULT '1' NOT NULL,
    translation_stable_id VARCHAR(128),
    translation_version USMALLINT DEFAULT '1' NOT NULL,
    peptide_archive_id UBIGINT,
    mapping_session_id UBIGINT NOT NULL
);


--
-- Name: gene_attrib; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.gene_attrib (
    gene_id UBIGINT DEFAULT '0' NOT NULL,
    attrib_type_id UINTEGER DEFAULT 0 NOT NULL,
    value VARCHAR NOT NULL,
    FOREIGN KEY (gene_id) REFERENCES ensembl_core_schema.gene (gene_id),
    FOREIGN KEY (attrib_type_id) REFERENCES ensembl_core_schema.attrib_type (attrib_type_id)
);


--
-- Name: gene_gene_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.gene_gene_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: genome_statistics; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.genome_statistics (
    genome_statistics_id UINTEGER PRIMARY KEY,
    statistic VARCHAR(128) NOT NULL,
    value bigint DEFAULT '0' NOT NULL,
    species_id UINTEGER DEFAULT '1',
    attrib_type_id UINTEGER,
    "timestamp" timestamp with time zone
);


--
-- Name: genome_statistics_genome_statistics_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.genome_statistics_genome_statistics_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: identity_xref; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.identity_xref (
    object_xref_id UBIGINT NOT NULL,
    xref_identity UINTEGER,
    ensembl_identity UINTEGER,
    xref_start UINTEGER,
    xref_end UINTEGER,
    ensembl_start UINTEGER,
    ensembl_end UINTEGER,
    cigar_line VARCHAR,
    score DOUBLE,
    evalue DOUBLE
);


--
-- Name: interpro; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.interpro (
    interpro_ac VARCHAR(40) NOT NULL,
    id VARCHAR(40) NOT NULL
);


--
-- Name: intron_supporting_evidence; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.intron_supporting_evidence (
    intron_supporting_evidence_id UBIGINT PRIMARY KEY,
    analysis_id UINTEGER NOT NULL,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand TINYINT NOT NULL,
    hit_name VARCHAR(100) NOT NULL,
    score numeric(10,3),
    score_type intron_supporting_evidence_score_type DEFAULT 'NONE',
    is_splice_canonical boolean DEFAULT false NOT NULL,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: intron_supporting_evidence_intron_supporting_evidence_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.intron_supporting_evidence_intron_supporting_evidence_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: karyotype; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.karyotype (
    karyotype_id UINTEGER PRIMARY KEY,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    band VARCHAR(40),
    stain VARCHAR(40),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: karyotype_karyotype_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.karyotype_karyotype_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: map; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.map (
    map_id UBIGINT PRIMARY KEY,
    map_name VARCHAR(30) NOT NULL
);


--
-- Name: map_map_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.map_map_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: mapping_session; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.mapping_session (
    mapping_session_id UBIGINT PRIMARY KEY,
    old_db_name VARCHAR(80) DEFAULT '' NOT NULL,
    new_db_name VARCHAR(80) DEFAULT '' NOT NULL,
    old_release VARCHAR(5) DEFAULT '' NOT NULL,
    new_release VARCHAR(5) DEFAULT '' NOT NULL,
    old_assembly VARCHAR(80) DEFAULT '' NOT NULL,
    new_assembly VARCHAR(80) DEFAULT '' NOT NULL,
    created timestamp with time zone NOT NULL
);


--
-- Name: mapping_session_mapping_session_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.mapping_session_mapping_session_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: mapping_set; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.mapping_set (
    mapping_set_id UBIGINT PRIMARY KEY,
    internal_schema_build VARCHAR(20) NOT NULL,
    external_schema_build VARCHAR(20) NOT NULL
);


--
-- Name: marker; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.marker (
    marker_id UBIGINT PRIMARY KEY,
    display_marker_synonym_id UBIGINT,
    left_primer VARCHAR(100) NOT NULL,
    right_primer VARCHAR(100) NOT NULL,
    min_primer_dist UBIGINT NOT NULL,
    max_primer_dist UBIGINT NOT NULL,
    priority UINTEGER,
    type marker_type
);


--
-- Name: marker_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.marker_feature (
    marker_feature_id UBIGINT PRIMARY KEY,
    marker_id UBIGINT NOT NULL,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    analysis_id UINTEGER NOT NULL,
    map_weight UBIGINT
);


--
-- Name: marker_feature_marker_feature_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.marker_feature_marker_feature_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: marker_map_location; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.marker_map_location (
    marker_id UBIGINT NOT NULL,
    map_id UBIGINT NOT NULL,
    chromosome_name VARCHAR(15) NOT NULL,
    marker_synonym_id UBIGINT NOT NULL,
    "position" VARCHAR(15) NOT NULL,
    lod_score DOUBLE
);


--
-- Name: marker_marker_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.marker_marker_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: marker_synonym; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.marker_synonym (
    marker_synonym_id UBIGINT PRIMARY KEY,
    marker_id UBIGINT NOT NULL,
    source VARCHAR(20),
    name VARCHAR(50)
);


--
-- Name: marker_synonym_marker_synonym_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.marker_synonym_marker_synonym_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: meta; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.meta (
    meta_id UINTEGER PRIMARY KEY,
    species_id UINTEGER DEFAULT '1',
    meta_key VARCHAR(40) NOT NULL,
    meta_value VARCHAR(255) NOT NULL
);


--
-- Name: meta_coord; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.meta_coord (
    table_name VARCHAR(40) NOT NULL,
    coord_system_id UINTEGER NOT NULL,
    max_length UINTEGER
);


--
-- Name: meta_meta_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.meta_meta_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: misc_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.misc_feature (
    misc_feature_id UBIGINT PRIMARY KEY,
    seq_region_id UBIGINT DEFAULT '0' NOT NULL,
    seq_region_start UBIGINT DEFAULT '0' NOT NULL,
    seq_region_end UBIGINT DEFAULT '0' NOT NULL,
    seq_region_strand TINYINT DEFAULT '0' NOT NULL,
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: misc_attrib; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.misc_attrib (
    misc_feature_id UBIGINT DEFAULT '0' NOT NULL,
    attrib_type_id UINTEGER DEFAULT 0 NOT NULL,
    value VARCHAR NOT NULL,
    FOREIGN KEY (misc_feature_id) REFERENCES ensembl_core_schema.misc_feature (misc_feature_id),
    FOREIGN KEY (attrib_type_id) REFERENCES ensembl_core_schema.attrib_type (attrib_type_id)
);



--
-- Name: misc_feature_misc_feature_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.misc_feature_misc_feature_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: misc_feature_misc_set; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.misc_feature_misc_set (
    misc_feature_id UBIGINT DEFAULT '0' NOT NULL,
    misc_set_id UINTEGER DEFAULT 0 NOT NULL
);


--
-- Name: misc_set; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.misc_set (
    misc_set_id UINTEGER PRIMARY KEY,
    code VARCHAR(25) DEFAULT '' NOT NULL,
    name VARCHAR(255) DEFAULT '' NOT NULL,
    description VARCHAR NOT NULL,
    max_length UBIGINT NOT NULL
);


--
-- Name: misc_set_misc_set_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.misc_set_misc_set_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: object_xref; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.object_xref (
    object_xref_id UBIGINT PRIMARY KEY,
    ensembl_id UBIGINT NOT NULL,
    ensembl_object_type object_xref_ensembl_object_type NOT NULL,
    xref_id UBIGINT NOT NULL,
    linkage_annotation VARCHAR(255),
    analysis_id UINTEGER
);


--
-- Name: object_xref_object_xref_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.object_xref_object_xref_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: ontology_xref; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.ontology_xref (
    object_xref_id UBIGINT DEFAULT '0' NOT NULL,
    source_xref_id UBIGINT,
    linkage_type VARCHAR(3)
);


--
-- Name: operon; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.operon (
    operon_id UBIGINT PRIMARY KEY,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand TINYINT NOT NULL,
    display_label VARCHAR(255),
    analysis_id UINTEGER NOT NULL,
    stable_id VARCHAR(128),
    version UINTEGER,
    created_date timestamp with time zone,
    modified_date timestamp with time zone,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: operon_operon_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.operon_operon_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: operon_transcript; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.operon_transcript (
    operon_transcript_id UBIGINT PRIMARY KEY,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand TINYINT NOT NULL,
    operon_id UBIGINT NOT NULL,
    display_label VARCHAR(255),
    analysis_id UINTEGER NOT NULL,
    stable_id VARCHAR(128),
    version UINTEGER,
    created_date timestamp with time zone,
    modified_date timestamp with time zone,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (operon_id) REFERENCES ensembl_core_schema.operon (operon_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: operon_transcript_gene; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.operon_transcript_gene (
    operon_transcript_id UBIGINT,
    gene_id UBIGINT,
    FOREIGN KEY (operon_transcript_id) REFERENCES ensembl_core_schema.operon_transcript (operon_transcript_id),
    FOREIGN KEY (gene_id) REFERENCES ensembl_core_schema.gene (gene_id)
);


--
-- Name: operon_transcript_operon_transcript_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.operon_transcript_operon_transcript_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: peptide_archive; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.peptide_archive (
    peptide_archive_id UBIGINT PRIMARY KEY,
    md5_checksum VARCHAR(32),
    peptide_seq VARCHAR NOT NULL
);


--
-- Name: peptide_archive_peptide_archive_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.peptide_archive_peptide_archive_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;



--
-- Name: prediction_transcript; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.prediction_transcript (
    prediction_transcript_id UBIGINT PRIMARY KEY,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand TINYINT NOT NULL,
    analysis_id UINTEGER NOT NULL,
    display_label VARCHAR(255),
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: prediction_exon; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.prediction_exon (
    prediction_exon_id UBIGINT PRIMARY KEY,
    prediction_transcript_id UBIGINT NOT NULL,
    exon_rank UINTEGER NOT NULL,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand TINYINT NOT NULL,
    start_phase TINYINT NOT NULL,
    score DOUBLE,
    p_value DOUBLE,
    FOREIGN KEY (prediction_transcript_id) REFERENCES ensembl_core_schema.prediction_transcript (prediction_transcript_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: prediction_exon_prediction_exon_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.prediction_exon_prediction_exon_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: prediction_transcript_prediction_transcript_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.prediction_transcript_prediction_transcript_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: protein_align_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.protein_align_feature (
    protein_align_feature_id UBIGINT PRIMARY KEY,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand boolean DEFAULT true NOT NULL,
    hit_start UINTEGER NOT NULL,
    hit_end UINTEGER NOT NULL,
    hit_name VARCHAR(40) NOT NULL,
    analysis_id UINTEGER NOT NULL,
    score DOUBLE,
    evalue DOUBLE,
    perc_ident DOUBLE,
    cigar_line VARCHAR,
    external_db_id UINTEGER,
    hcoverage DOUBLE,
    align_type protein_align_feature_align_type DEFAULT 'ensembl',
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: protein_align_feature_protein_align_feature_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.protein_align_feature_protein_align_feature_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: translation; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.translation (
    translation_id UBIGINT PRIMARY KEY,
    transcript_id UBIGINT NOT NULL,
    seq_start UINTEGER NOT NULL,
    start_exon_id UBIGINT NOT NULL,
    seq_end UINTEGER NOT NULL,
    end_exon_id UBIGINT NOT NULL,
    stable_id VARCHAR(128),
    version UINTEGER,
    created_date timestamp with time zone,
    modified_date timestamp with time zone,
    FOREIGN KEY (transcript_id) REFERENCES ensembl_core_schema.transcript (transcript_id),
    FOREIGN KEY (start_exon_id) REFERENCES ensembl_core_schema.exon (exon_id),
    FOREIGN KEY (end_exon_id) REFERENCES ensembl_core_schema.exon (exon_id)
);


--
-- Name: translation_attrib; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.translation_attrib (
    translation_id UBIGINT DEFAULT 0 NOT NULL,
    attrib_type_id UINTEGER DEFAULT 0 NOT NULL,
    value VARCHAR NOT NULL,
    FOREIGN KEY (translation_id) REFERENCES ensembl_core_schema.translation (translation_id),
    FOREIGN KEY (attrib_type_id) REFERENCES ensembl_core_schema.attrib_type (attrib_type_id)
);


--
-- Name: protein_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.protein_feature (
    protein_feature_id UBIGINT PRIMARY KEY,
    translation_id UBIGINT NOT NULL,
    seq_start UINTEGER NOT NULL,
    seq_end UINTEGER NOT NULL,
    hit_start UINTEGER NOT NULL,
    hit_end UINTEGER NOT NULL,
    hit_name VARCHAR(40) NOT NULL,
    analysis_id UINTEGER NOT NULL,
    score DOUBLE,
    evalue DOUBLE,
    perc_ident DOUBLE,
    external_data VARCHAR,
    hit_description VARCHAR,
    cigar_line VARCHAR,
    align_type protein_feature_align_type,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (translation_id) REFERENCES ensembl_core_schema.translation (translation_id)
);


--
-- Name: protein_feature_protein_feature_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.protein_feature_protein_feature_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: repeat_consensus; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.repeat_consensus (
    repeat_consensus_id UBIGINT PRIMARY KEY,
    repeat_name VARCHAR(255) NOT NULL,
    repeat_class VARCHAR(100) NOT NULL,
    repeat_type VARCHAR(40) NOT NULL,
    repeat_consensus VARCHAR
);


--
-- Name: repeat_consensus_repeat_consensus_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.repeat_consensus_repeat_consensus_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: repeat_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.repeat_feature (
    repeat_feature_id UBIGINT PRIMARY KEY,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand boolean DEFAULT true NOT NULL,
    repeat_start UINTEGER NOT NULL,
    repeat_end UINTEGER NOT NULL,
    repeat_consensus_id UBIGINT NOT NULL,
    analysis_id UINTEGER NOT NULL,
    score DOUBLE,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: repeat_feature_repeat_feature_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.repeat_feature_repeat_feature_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: seq_region_attrib; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.seq_region_attrib (
    seq_region_id UBIGINT DEFAULT '0' NOT NULL,
    attrib_type_id UINTEGER DEFAULT 0 NOT NULL,
    value VARCHAR NOT NULL,
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id),
    FOREIGN KEY (attrib_type_id) REFERENCES ensembl_core_schema.attrib_type (attrib_type_id)
);


--
-- Name: seq_region_mapping; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.seq_region_mapping (
    external_seq_region_id UBIGINT NOT NULL,
    internal_seq_region_id UBIGINT NOT NULL,
    mapping_set_id UBIGINT NOT NULL
);


--
-- Name: seq_region_seq_region_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.seq_region_seq_region_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: seq_region_synonym; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.seq_region_synonym (
    seq_region_synonym_id UBIGINT PRIMARY KEY,
    seq_region_id UBIGINT NOT NULL,
    synonym VARCHAR(250) NOT NULL,
    external_db_id UINTEGER,
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: seq_region_synonym_seq_region_synonym_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.seq_region_synonym_seq_region_synonym_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: simple_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.simple_feature (
    simple_feature_id UBIGINT PRIMARY KEY,
    seq_region_id UBIGINT NOT NULL,
    seq_region_start UBIGINT NOT NULL,
    seq_region_end UBIGINT NOT NULL,
    seq_region_strand boolean NOT NULL,
    display_label VARCHAR(255) NOT NULL,
    analysis_id UINTEGER NOT NULL,
    score DOUBLE,
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id),
    FOREIGN KEY (seq_region_id) REFERENCES ensembl_core_schema.seq_region (seq_region_id)
);


--
-- Name: simple_feature_simple_feature_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.simple_feature_simple_feature_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: stable_id_event; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.stable_id_event (
    old_stable_id VARCHAR(128),
    old_version USMALLINT,
    new_stable_id VARCHAR(128),
    new_version USMALLINT,
    mapping_session_id UBIGINT DEFAULT '0' NOT NULL,
    type stable_id_event_type NOT NULL,
    score DOUBLE DEFAULT '0' NOT NULL
);


--
-- Name: supporting_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.supporting_feature (
    exon_id UBIGINT DEFAULT '0' NOT NULL,
    feature_type supporting_feature_feature_type,
    feature_id UBIGINT DEFAULT '0' NOT NULL,
    FOREIGN KEY (exon_id) REFERENCES ensembl_core_schema.exon (exon_id)
);


--
-- Name: transcript_attrib; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.transcript_attrib (
    transcript_id UBIGINT DEFAULT '0' NOT NULL,
    attrib_type_id UINTEGER DEFAULT 0 NOT NULL,
    value VARCHAR NOT NULL,
    FOREIGN KEY (transcript_id) REFERENCES ensembl_core_schema.transcript (transcript_id),
    FOREIGN KEY (attrib_type_id) REFERENCES ensembl_core_schema.attrib_type (attrib_type_id)
);


--
-- Name: transcript_intron_supporting_evidence; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.transcript_intron_supporting_evidence (
    transcript_id UBIGINT NOT NULL,
    intron_supporting_evidence_id UBIGINT NOT NULL,
    previous_exon_id UBIGINT NOT NULL,
    next_exon_id UBIGINT NOT NULL,
    FOREIGN KEY (transcript_id) REFERENCES ensembl_core_schema.transcript (transcript_id)
);


--
-- Name: transcript_supporting_feature; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.transcript_supporting_feature (
    transcript_id UBIGINT DEFAULT '0' NOT NULL,
    feature_type transcript_supporting_feature_feature_type,
    feature_id UBIGINT DEFAULT '0' NOT NULL,
    FOREIGN KEY (transcript_id) REFERENCES ensembl_core_schema.transcript (transcript_id)
);


--
-- Name: transcript_transcript_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.transcript_transcript_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: translation_translation_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.translation_translation_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: unmapped_object; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.unmapped_object (
    unmapped_object_id UBIGINT PRIMARY KEY,
    type unmapped_object_type NOT NULL,
    analysis_id UINTEGER NOT NULL,
    external_db_id UINTEGER,
    identifier VARCHAR(255) NOT NULL,
    unmapped_reason_id UBIGINT NOT NULL,
    query_score DOUBLE,
    target_score DOUBLE,
    ensembl_id UBIGINT DEFAULT '0',
    ensembl_object_type unmapped_object_ensembl_object_type DEFAULT 'RawContig',
    parent VARCHAR(255),
    FOREIGN KEY (analysis_id) REFERENCES ensembl_core_schema.analysis (analysis_id)
);


--
-- Name: unmapped_object_unmapped_object_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.unmapped_object_unmapped_object_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: unmapped_reason; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.unmapped_reason (
    unmapped_reason_id UBIGINT NOT NULL,
    summary_description VARCHAR(255),
    full_description VARCHAR(255)
);


--
-- Name: unmapped_reason_unmapped_reason_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.unmapped_reason_unmapped_reason_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: xref; Type: TABLE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE TABLE ensembl_core_schema.xref (
    xref_id UBIGINT PRIMARY KEY,
    external_db_id UINTEGER NOT NULL,
    dbprimary_acc VARCHAR(512) NOT NULL,
    display_label VARCHAR(512) NOT NULL,
    version VARCHAR(10),
    description VARCHAR,
    info_type xref_info_type DEFAULT 'NONE',
    info_text VARCHAR(255) DEFAULT '' NOT NULL
);


--
-- Name: xref_xref_id_seq; Type: SEQUENCE; Schema: ensembl_core_schema; Owner: stefano
--

CREATE SEQUENCE ensembl_core_schema.xref_xref_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE;


--
-- Name: alt_allele alt_allele_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.alt_allele ALTER COLUMN alt_allele_id SET DEFAULT nextval('ensembl_core_schema.alt_allele_alt_allele_id_seq');


--
-- Name: analysis analysis_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.analysis ALTER COLUMN analysis_id SET DEFAULT nextval('ensembl_core_schema.analysis_analysis_id_seq');


--
-- Name: associated_group associated_group_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.associated_group ALTER COLUMN associated_group_id SET DEFAULT nextval('ensembl_core_schema.associated_group_associated_group_id_seq');


--
-- Name: associated_xref associated_xref_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.associated_xref ALTER COLUMN associated_xref_id SET DEFAULT nextval('ensembl_core_schema.associated_xref_associated_xref_id_seq');


--
-- Name: attrib_type attrib_type_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.attrib_type ALTER COLUMN attrib_type_id SET DEFAULT nextval('ensembl_core_schema.attrib_type_attrib_type_id_seq');


--
-- Name: biotype biotype_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.biotype ALTER COLUMN biotype_id SET DEFAULT nextval('ensembl_core_schema.biotype_biotype_id_seq');


--
-- Name: coord_system coord_system_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.coord_system ALTER COLUMN coord_system_id SET DEFAULT nextval('ensembl_core_schema.coord_system_coord_system_id_seq');


--
-- Name: data_file data_file_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.data_file ALTER COLUMN data_file_id SET DEFAULT nextval('ensembl_core_schema.data_file_data_file_id_seq');


--
-- Name: density_feature density_feature_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.density_feature ALTER COLUMN density_feature_id SET DEFAULT nextval('ensembl_core_schema.density_feature_density_feature_id_seq');


--
-- Name: density_type density_type_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.density_type ALTER COLUMN density_type_id SET DEFAULT nextval('ensembl_core_schema.density_type_density_type_id_seq');


--
-- Name: ditag ditag_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.ditag ALTER COLUMN ditag_id SET DEFAULT nextval('ensembl_core_schema.ditag_ditag_id_seq');


--
-- Name: ditag_feature ditag_feature_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.ditag_feature ALTER COLUMN ditag_feature_id SET DEFAULT nextval('ensembl_core_schema.ditag_feature_ditag_feature_id_seq');


--
-- Name: dna_align_feature dna_align_feature_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.dna_align_feature ALTER COLUMN dna_align_feature_id SET DEFAULT nextval('ensembl_core_schema.dna_align_feature_dna_align_feature_id_seq');


--
-- Name: exon exon_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.exon ALTER COLUMN exon_id SET DEFAULT nextval('ensembl_core_schema.exon_exon_id_seq');


--
-- Name: external_db external_db_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.external_db ALTER COLUMN external_db_id SET DEFAULT nextval('ensembl_core_schema.external_db_external_db_id_seq');


--
-- Name: gene gene_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.gene ALTER COLUMN gene_id SET DEFAULT nextval('ensembl_core_schema.gene_gene_id_seq');


--
-- Name: genome_statistics genome_statistics_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.genome_statistics ALTER COLUMN genome_statistics_id SET DEFAULT nextval('ensembl_core_schema.genome_statistics_genome_statistics_id_seq');


--
-- Name: intron_supporting_evidence intron_supporting_evidence_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.intron_supporting_evidence ALTER COLUMN intron_supporting_evidence_id SET DEFAULT nextval('ensembl_core_schema.intron_supporting_evidence_intron_supporting_evidence_id_seq');


--
-- Name: karyotype karyotype_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.karyotype ALTER COLUMN karyotype_id SET DEFAULT nextval('ensembl_core_schema.karyotype_karyotype_id_seq');


--
-- Name: map map_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.map ALTER COLUMN map_id SET DEFAULT nextval('ensembl_core_schema.map_map_id_seq');


--
-- Name: mapping_session mapping_session_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.mapping_session ALTER COLUMN mapping_session_id SET DEFAULT nextval('ensembl_core_schema.mapping_session_mapping_session_id_seq');


--
-- Name: marker marker_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.marker ALTER COLUMN marker_id SET DEFAULT nextval('ensembl_core_schema.marker_marker_id_seq');


--
-- Name: marker_feature marker_feature_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.marker_feature ALTER COLUMN marker_feature_id SET DEFAULT nextval('ensembl_core_schema.marker_feature_marker_feature_id_seq');


--
-- Name: marker_synonym marker_synonym_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.marker_synonym ALTER COLUMN marker_synonym_id SET DEFAULT nextval('ensembl_core_schema.marker_synonym_marker_synonym_id_seq');


--
-- Name: meta meta_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.meta ALTER COLUMN meta_id SET DEFAULT nextval('ensembl_core_schema.meta_meta_id_seq');


--
-- Name: misc_feature misc_feature_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.misc_feature ALTER COLUMN misc_feature_id SET DEFAULT nextval('ensembl_core_schema.misc_feature_misc_feature_id_seq');


--
-- Name: misc_set misc_set_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.misc_set ALTER COLUMN misc_set_id SET DEFAULT nextval('ensembl_core_schema.misc_set_misc_set_id_seq');


--
-- Name: object_xref object_xref_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.object_xref ALTER COLUMN object_xref_id SET DEFAULT nextval('ensembl_core_schema.object_xref_object_xref_id_seq');


--
-- Name: operon operon_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.operon ALTER COLUMN operon_id SET DEFAULT nextval('ensembl_core_schema.operon_operon_id_seq');


--
-- Name: operon_transcript operon_transcript_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.operon_transcript ALTER COLUMN operon_transcript_id SET DEFAULT nextval('ensembl_core_schema.operon_transcript_operon_transcript_id_seq');


--
-- Name: peptide_archive peptide_archive_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.peptide_archive ALTER COLUMN peptide_archive_id SET DEFAULT nextval('ensembl_core_schema.peptide_archive_peptide_archive_id_seq');


--
-- Name: prediction_exon prediction_exon_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.prediction_exon ALTER COLUMN prediction_exon_id SET DEFAULT nextval('ensembl_core_schema.prediction_exon_prediction_exon_id_seq');


--
-- Name: prediction_transcript prediction_transcript_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.prediction_transcript ALTER COLUMN prediction_transcript_id SET DEFAULT nextval('ensembl_core_schema.prediction_transcript_prediction_transcript_id_seq');


--
-- Name: protein_align_feature protein_align_feature_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.protein_align_feature ALTER COLUMN protein_align_feature_id SET DEFAULT nextval('ensembl_core_schema.protein_align_feature_protein_align_feature_id_seq');


--
-- Name: protein_feature protein_feature_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.protein_feature ALTER COLUMN protein_feature_id SET DEFAULT nextval('ensembl_core_schema.protein_feature_protein_feature_id_seq');


--
-- Name: repeat_consensus repeat_consensus_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.repeat_consensus ALTER COLUMN repeat_consensus_id SET DEFAULT nextval('ensembl_core_schema.repeat_consensus_repeat_consensus_id_seq');


--
-- Name: repeat_feature repeat_feature_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.repeat_feature ALTER COLUMN repeat_feature_id SET DEFAULT nextval('ensembl_core_schema.repeat_feature_repeat_feature_id_seq');


--
-- Name: seq_region seq_region_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.seq_region ALTER COLUMN seq_region_id SET DEFAULT nextval('ensembl_core_schema.seq_region_seq_region_id_seq');


--
-- Name: seq_region_synonym seq_region_synonym_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.seq_region_synonym ALTER COLUMN seq_region_synonym_id SET DEFAULT nextval('ensembl_core_schema.seq_region_synonym_seq_region_synonym_id_seq');


--
-- Name: simple_feature simple_feature_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.simple_feature ALTER COLUMN simple_feature_id SET DEFAULT nextval('ensembl_core_schema.simple_feature_simple_feature_id_seq');


--
-- Name: transcript transcript_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.transcript ALTER COLUMN transcript_id SET DEFAULT nextval('ensembl_core_schema.transcript_transcript_id_seq');


--
-- Name: translation translation_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.translation ALTER COLUMN translation_id SET DEFAULT nextval('ensembl_core_schema.translation_translation_id_seq');


--
-- Name: unmapped_object unmapped_object_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.unmapped_object ALTER COLUMN unmapped_object_id SET DEFAULT nextval('ensembl_core_schema.unmapped_object_unmapped_object_id_seq');


--
-- Name: unmapped_reason unmapped_reason_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.unmapped_reason ALTER COLUMN unmapped_reason_id SET DEFAULT nextval('ensembl_core_schema.unmapped_reason_unmapped_reason_id_seq');


--
-- Name: xref xref_id; Type: DEFAULT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.xref ALTER COLUMN xref_id SET DEFAULT nextval('ensembl_core_schema.xref_xref_id_seq');


/*
--
-- Name: alt_allele idx_16668_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.alt_allele
    ADD CONSTRAINT idx_16668_primary PRIMARY KEY (alt_allele_id);


--
-- Name: alt_allele_group idx_16676_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.alt_allele_group
    ADD CONSTRAINT idx_16676_primary PRIMARY KEY (alt_allele_group_id);


--
-- Name: analysis idx_16681_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.analysis
    ADD CONSTRAINT idx_16681_primary PRIMARY KEY (analysis_id);


--
-- Name: associated_group idx_16702_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.associated_group
    ADD CONSTRAINT idx_16702_primary PRIMARY KEY (associated_group_id);


--
-- Name: associated_xref idx_16707_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.associated_xref
    ADD CONSTRAINT idx_16707_primary PRIMARY KEY (associated_xref_id);


--
-- Name: attrib_type idx_16715_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.attrib_type
    ADD CONSTRAINT idx_16715_primary PRIMARY KEY (attrib_type_id);


--
-- Name: biotype idx_16724_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.biotype
    ADD CONSTRAINT idx_16724_primary PRIMARY KEY (biotype_id);


--
-- Name: coord_system idx_16733_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.coord_system
    ADD CONSTRAINT idx_16733_primary PRIMARY KEY (coord_system_id);


--
-- Name: data_file idx_16741_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.data_file
    ADD CONSTRAINT idx_16741_primary PRIMARY KEY (data_file_id);


--
-- Name: density_feature idx_16750_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.density_feature
    ADD CONSTRAINT idx_16750_primary PRIMARY KEY (density_feature_id);


--
-- Name: density_type idx_16755_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.density_type
    ADD CONSTRAINT idx_16755_primary PRIMARY KEY (density_type_id);


--
-- Name: dependent_xref idx_16759_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.dependent_xref
    ADD CONSTRAINT idx_16759_primary PRIMARY KEY (object_xref_id);


--
-- Name: ditag idx_16763_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.ditag
    ADD CONSTRAINT idx_16763_primary PRIMARY KEY (ditag_id);


--
-- Name: ditag_feature idx_16771_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.ditag_feature
    ADD CONSTRAINT idx_16771_primary PRIMARY KEY (ditag_feature_id);


--
-- Name: dna idx_16787_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.dna
    ADD CONSTRAINT idx_16787_primary PRIMARY KEY (seq_region_id);


--
-- Name: dna_align_feature idx_16793_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.dna_align_feature
    ADD CONSTRAINT idx_16793_primary PRIMARY KEY (dna_align_feature_id);


--
-- Name: exon idx_16806_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.exon
    ADD CONSTRAINT idx_16806_primary PRIMARY KEY (exon_id);


--
-- Name: exon_transcript idx_16812_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.exon_transcript
    ADD CONSTRAINT idx_16812_primary PRIMARY KEY (exon_id, transcript_id, rank);


--
-- Name: external_db idx_16816_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.external_db
    ADD CONSTRAINT idx_16816_primary PRIMARY KEY (external_db_id);


--
-- Name: external_synonym idx_16822_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.external_synonym
    ADD CONSTRAINT idx_16822_primary PRIMARY KEY (xref_id, synonym);


--
-- Name: gene idx_16826_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.gene
    ADD CONSTRAINT idx_16826_primary PRIMARY KEY (gene_id);


--
-- Name: genome_statistics idx_16847_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.genome_statistics
    ADD CONSTRAINT idx_16847_primary PRIMARY KEY (genome_statistics_id);


--
-- Name: identity_xref idx_16853_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.identity_xref
    ADD CONSTRAINT idx_16853_primary PRIMARY KEY (object_xref_id);


--
-- Name: intron_supporting_evidence idx_16862_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.intron_supporting_evidence
    ADD CONSTRAINT idx_16862_primary PRIMARY KEY (intron_supporting_evidence_id);


--
-- Name: karyotype idx_16869_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.karyotype
    ADD CONSTRAINT idx_16869_primary PRIMARY KEY (karyotype_id);


--
-- Name: map idx_16874_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.map
    ADD CONSTRAINT idx_16874_primary PRIMARY KEY (map_id);


--
-- Name: mapping_session idx_16879_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.mapping_session
    ADD CONSTRAINT idx_16879_primary PRIMARY KEY (mapping_session_id);


--
-- Name: mapping_set idx_16889_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.mapping_set
    ADD CONSTRAINT idx_16889_primary PRIMARY KEY (mapping_set_id);


--
-- Name: marker idx_16893_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.marker
    ADD CONSTRAINT idx_16893_primary PRIMARY KEY (marker_id);


--
-- Name: marker_feature idx_16898_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.marker_feature
    ADD CONSTRAINT idx_16898_primary PRIMARY KEY (marker_feature_id);


--
-- Name: marker_map_location idx_16902_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.marker_map_location
    ADD CONSTRAINT idx_16902_primary PRIMARY KEY (marker_id, map_id);


--
-- Name: marker_synonym idx_16906_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.marker_synonym
    ADD CONSTRAINT idx_16906_primary PRIMARY KEY (marker_synonym_id);


--
-- Name: meta idx_16911_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.meta
    ADD CONSTRAINT idx_16911_primary PRIMARY KEY (meta_id);


--
-- Name: misc_feature idx_16927_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.misc_feature
    ADD CONSTRAINT idx_16927_primary PRIMARY KEY (misc_feature_id);


--
-- Name: misc_feature_misc_set idx_16935_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.misc_feature_misc_set
    ADD CONSTRAINT idx_16935_primary PRIMARY KEY (misc_feature_id, misc_set_id);


--
-- Name: misc_set idx_16941_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.misc_set
    ADD CONSTRAINT idx_16941_primary PRIMARY KEY (misc_set_id);


--
-- Name: object_xref idx_16950_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.object_xref
    ADD CONSTRAINT idx_16950_primary PRIMARY KEY (object_xref_id);


--
-- Name: operon idx_16959_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.operon
    ADD CONSTRAINT idx_16959_primary PRIMARY KEY (operon_id);


--
-- Name: operon_transcript idx_16964_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.operon_transcript
    ADD CONSTRAINT idx_16964_primary PRIMARY KEY (operon_transcript_id);


--
-- Name: peptide_archive idx_16972_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.peptide_archive
    ADD CONSTRAINT idx_16972_primary PRIMARY KEY (peptide_archive_id);


--
-- Name: prediction_exon idx_16979_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.prediction_exon
    ADD CONSTRAINT idx_16979_primary PRIMARY KEY (prediction_exon_id);


--
-- Name: prediction_transcript idx_16984_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.prediction_transcript
    ADD CONSTRAINT idx_16984_primary PRIMARY KEY (prediction_transcript_id);


--
-- Name: protein_align_feature idx_16989_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.protein_align_feature
    ADD CONSTRAINT idx_16989_primary PRIMARY KEY (protein_align_feature_id);


--
-- Name: protein_feature idx_16998_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.protein_feature
    ADD CONSTRAINT idx_16998_primary PRIMARY KEY (protein_feature_id);


--
-- Name: repeat_consensus idx_17005_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.repeat_consensus
    ADD CONSTRAINT idx_17005_primary PRIMARY KEY (repeat_consensus_id);


--
-- Name: repeat_feature idx_17012_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.repeat_feature
    ADD CONSTRAINT idx_17012_primary PRIMARY KEY (repeat_feature_id);


--
-- Name: seq_region idx_17039_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.seq_region
    ADD CONSTRAINT idx_17039_primary PRIMARY KEY (seq_region_id);


--
-- Name: seq_region_synonym idx_17054_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.seq_region_synonym
    ADD CONSTRAINT idx_17054_primary PRIMARY KEY (seq_region_synonym_id);


--
-- Name: simple_feature idx_17059_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.simple_feature
    ADD CONSTRAINT idx_17059_primary PRIMARY KEY (simple_feature_id);


--
-- Name: transcript idx_17074_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.transcript
    ADD CONSTRAINT idx_17074_primary PRIMARY KEY (transcript_id);


--
-- Name: transcript_intron_supporting_evidence idx_17089_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.transcript_intron_supporting_evidence
    ADD CONSTRAINT idx_17089_primary PRIMARY KEY (intron_supporting_evidence_id, transcript_id);


--
-- Name: translation idx_17098_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.translation
    ADD CONSTRAINT idx_17098_primary PRIMARY KEY (translation_id);


--
-- Name: unmapped_object idx_17110_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.unmapped_object
    ADD CONSTRAINT idx_17110_primary PRIMARY KEY (unmapped_object_id);


--
-- Name: unmapped_reason idx_17119_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.unmapped_reason
    ADD CONSTRAINT idx_17119_primary PRIMARY KEY (unmapped_reason_id);


--
-- Name: xref idx_17126_primary; Type: CONSTRAINT; Schema: ensembl_core_schema; Owner: stefano
--

ALTER TABLE ensembl_core_schema.xref
    ADD CONSTRAINT idx_17126_primary PRIMARY KEY (xref_id);

*/

--
-- Name: idx_16668_gene_id; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16668_gene_id ON ensembl_core_schema.alt_allele (gene_id, alt_allele_group_id);


--
-- Name: idx_16668_gene_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16668_gene_idx ON ensembl_core_schema.alt_allele (gene_id);


--
-- Name: idx_16672_aa_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16672_aa_idx ON ensembl_core_schema.alt_allele_attrib (alt_allele_id, attrib);


--
-- Name: idx_16681_logic_name_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16681_logic_name_idx ON ensembl_core_schema.analysis (logic_name);


--
-- Name: idx_16687_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16687_analysis_idx ON ensembl_core_schema.analysis_description (analysis_id);


--
-- Name: idx_16693_all_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16693_all_idx ON ensembl_core_schema.assembly (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start, cmp_end, ori);


--
-- Name: idx_16693_asm_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16693_asm_seq_region_idx ON ensembl_core_schema.assembly (asm_seq_region_id, asm_start);


--
-- Name: idx_16693_cmp_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16693_cmp_seq_region_idx ON ensembl_core_schema.assembly (cmp_seq_region_id);


--
-- Name: idx_16707_associated_group_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16707_associated_group_idx ON ensembl_core_schema.associated_xref (associated_group_id);


--
-- Name: idx_16707_associated_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16707_associated_idx ON ensembl_core_schema.associated_xref (xref_id);


--
-- Name: idx_16707_associated_object_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16707_associated_object_idx ON ensembl_core_schema.associated_xref (object_xref_id);


--
-- Name: idx_16707_associated_source_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16707_associated_source_idx ON ensembl_core_schema.associated_xref (source_xref_id);


--
-- Name: idx_16707_object_associated_source_type_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16707_object_associated_source_type_idx ON ensembl_core_schema.associated_xref (object_xref_id, xref_id, source_xref_id, condition_type, associated_group_id);


--
-- Name: idx_16715_code_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16715_code_idx ON ensembl_core_schema.attrib_type (code);


--
-- Name: idx_16724_name_type_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16724_name_type_idx ON ensembl_core_schema.biotype (name, object_type);


--
-- Name: idx_16733_name_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16733_name_idx ON ensembl_core_schema.coord_system (name, version, species_id);


--
-- Name: idx_16733_rank_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16733_rank_idx ON ensembl_core_schema.coord_system (rank, species_id);


--
-- Name: idx_16733_species_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16733_species_idx ON ensembl_core_schema.coord_system (species_id);


--
-- Name: idx_16741_df_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16741_df_analysis_idx ON ensembl_core_schema.data_file (analysis_id);


--
-- Name: idx_16741_df_name_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16741_df_name_idx ON ensembl_core_schema.data_file (name);


--
-- Name: idx_16741_df_unq_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16741_df_unq_idx ON ensembl_core_schema.data_file (coord_system_id, analysis_id, name, file_type);


--
-- Name: idx_16750_seq_region_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16750_seq_region_id_idx ON ensembl_core_schema.density_feature (seq_region_id);


--
-- Name: idx_16750_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16750_seq_region_idx ON ensembl_core_schema.density_feature (density_type_id, seq_region_id, seq_region_start);


--
-- Name: idx_16755_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16755_analysis_idx ON ensembl_core_schema.density_type (analysis_id, block_size, region_features);


--
-- Name: idx_16759_dependent; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16759_dependent ON ensembl_core_schema.dependent_xref (dependent_xref_id);


--
-- Name: idx_16759_master_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16759_master_idx ON ensembl_core_schema.dependent_xref (master_xref_id);


--
-- Name: idx_16771_ditag_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16771_ditag_idx ON ensembl_core_schema.ditag_feature (ditag_id);


--
-- Name: idx_16771_ditag_pair_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16771_ditag_pair_idx ON ensembl_core_schema.ditag_feature (ditag_pair_id);


--
-- Name: idx_16771_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16771_seq_region_idx ON ensembl_core_schema.ditag_feature (seq_region_id, seq_region_start, seq_region_end);


--
-- Name: idx_16793_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16793_analysis_idx ON ensembl_core_schema.dna_align_feature (analysis_id);


--
-- Name: idx_16793_external_db_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16793_external_db_idx ON ensembl_core_schema.dna_align_feature (external_db_id);


--
-- Name: idx_16793_hit_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16793_hit_idx ON ensembl_core_schema.dna_align_feature (hit_name);


--
-- Name: idx_16793_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16793_seq_region_idx ON ensembl_core_schema.dna_align_feature (seq_region_id, analysis_id, seq_region_start, score);


--
-- Name: idx_16793_seq_region_idx_2; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16793_seq_region_idx_2 ON ensembl_core_schema.dna_align_feature (seq_region_id, seq_region_start);


--
-- Name: idx_16800_dna_align_feature_attribx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16800_dna_align_feature_attribx ON ensembl_core_schema.dna_align_feature_attrib (dna_align_feature_id, attrib_type_id, value);


--
-- Name: idx_16800_dna_align_feature_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16800_dna_align_feature_idx ON ensembl_core_schema.dna_align_feature_attrib (dna_align_feature_id);


--
-- Name: idx_16800_type_val_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16800_type_val_idx ON ensembl_core_schema.dna_align_feature_attrib (attrib_type_id, value);


--
-- Name: idx_16800_val_only_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16800_val_only_idx ON ensembl_core_schema.dna_align_feature_attrib (value);


--
-- Name: idx_16806_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16806_seq_region_idx ON ensembl_core_schema.exon (seq_region_id, seq_region_start);


--
-- Name: idx_16806_stable_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16806_stable_id_idx ON ensembl_core_schema.exon (stable_id, version);


--
-- Name: idx_16812_exon; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16812_exon ON ensembl_core_schema.exon_transcript (exon_id);


--
-- Name: idx_16812_transcript; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16812_transcript ON ensembl_core_schema.exon_transcript (transcript_id);


--
-- Name: idx_16816_db_name_db_release_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16816_db_name_db_release_idx ON ensembl_core_schema.external_db (db_name, db_release);


--
-- Name: idx_16822_name_index; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16822_name_index ON ensembl_core_schema.external_synonym (synonym);


--
-- Name: idx_16826_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16826_analysis_idx ON ensembl_core_schema.gene (analysis_id);


--
-- Name: idx_16826_canonical_transcript_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16826_canonical_transcript_id_idx ON ensembl_core_schema.gene (canonical_transcript_id);


--
-- Name: idx_16826_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16826_seq_region_idx ON ensembl_core_schema.gene (seq_region_id, seq_region_start);


--
-- Name: idx_16826_stable_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16826_stable_id_idx ON ensembl_core_schema.gene (stable_id, version);


--
-- Name: idx_16826_xref_id_index; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16826_xref_id_index ON ensembl_core_schema.gene (display_xref_id);


--
-- Name: idx_16833_gene_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16833_gene_idx ON ensembl_core_schema.gene_archive (gene_stable_id, gene_version);


--
-- Name: idx_16833_peptide_archive_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16833_peptide_archive_id_idx ON ensembl_core_schema.gene_archive (peptide_archive_id);


--
-- Name: idx_16833_transcript_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16833_transcript_idx ON ensembl_core_schema.gene_archive (transcript_stable_id, transcript_version);


--
-- Name: idx_16833_translation_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16833_translation_idx ON ensembl_core_schema.gene_archive (translation_stable_id, translation_version);


--
-- Name: idx_16839_gene_attribx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16839_gene_attribx ON ensembl_core_schema.gene_attrib (gene_id, attrib_type_id, value);


--
-- Name: idx_16839_gene_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16839_gene_idx ON ensembl_core_schema.gene_attrib (gene_id);


--
-- Name: idx_16839_type_val_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16839_type_val_idx ON ensembl_core_schema.gene_attrib (attrib_type_id, value);


--
-- Name: idx_16839_val_only_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16839_val_only_idx ON ensembl_core_schema.gene_attrib (value);


--
-- Name: idx_16847_stats_uniq; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16847_stats_uniq ON ensembl_core_schema.genome_statistics (statistic, attrib_type_id, species_id);


--
-- Name: idx_16858_accession_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16858_accession_idx ON ensembl_core_schema.interpro (interpro_ac, id);


--
-- Name: idx_16858_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16858_id_idx ON ensembl_core_schema.interpro (id);


--
-- Name: idx_16862_analysis_id; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16862_analysis_id ON ensembl_core_schema.intron_supporting_evidence (analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name);


--
-- Name: idx_16862_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16862_seq_region_idx ON ensembl_core_schema.intron_supporting_evidence (seq_region_id, seq_region_start);


--
-- Name: idx_16869_region_band_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16869_region_band_idx ON ensembl_core_schema.karyotype (seq_region_id, band);


--
-- Name: idx_16889_mapping_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16889_mapping_idx ON ensembl_core_schema.mapping_set (internal_schema_build, external_schema_build);


--
-- Name: idx_16893_display_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16893_display_idx ON ensembl_core_schema.marker (display_marker_synonym_id);


--
-- Name: idx_16893_marker_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16893_marker_idx ON ensembl_core_schema.marker (marker_id, priority);


--
-- Name: idx_16898_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16898_analysis_idx ON ensembl_core_schema.marker_feature (analysis_id);


--
-- Name: idx_16898_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16898_seq_region_idx ON ensembl_core_schema.marker_feature (seq_region_id, seq_region_start);


--
-- Name: idx_16902_map_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16902_map_idx ON ensembl_core_schema.marker_map_location (map_id, chromosome_name, "position");


--
-- Name: idx_16906_marker_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16906_marker_idx ON ensembl_core_schema.marker_synonym (marker_id);


--
-- Name: idx_16906_marker_synonym_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16906_marker_synonym_idx ON ensembl_core_schema.marker_synonym (marker_synonym_id, name);


--
-- Name: idx_16911_species_key_value_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16911_species_key_value_idx ON ensembl_core_schema.meta (species_id, meta_key, meta_value);


--
-- Name: idx_16911_species_value_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16911_species_value_idx ON ensembl_core_schema.meta (species_id, meta_value);


--
-- Name: idx_16916_cs_table_name_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16916_cs_table_name_idx ON ensembl_core_schema.meta_coord (coord_system_id, table_name);


--
-- Name: idx_16919_misc_attribx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16919_misc_attribx ON ensembl_core_schema.misc_attrib (misc_feature_id, attrib_type_id, value);


--
-- Name: idx_16919_misc_feature_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16919_misc_feature_idx ON ensembl_core_schema.misc_attrib (misc_feature_id);


--
-- Name: idx_16919_type_val_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16919_type_val_idx ON ensembl_core_schema.misc_attrib (attrib_type_id, value);


--
-- Name: idx_16919_val_only_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16919_val_only_idx ON ensembl_core_schema.misc_attrib (value);


--
-- Name: idx_16927_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16927_seq_region_idx ON ensembl_core_schema.misc_feature (seq_region_id, seq_region_start);


--
-- Name: idx_16935_reverse_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16935_reverse_idx ON ensembl_core_schema.misc_feature_misc_set (misc_set_id, misc_feature_id);


--
-- Name: idx_16941_code_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16941_code_idx ON ensembl_core_schema.misc_set (code);


--
-- Name: idx_16950_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16950_analysis_idx ON ensembl_core_schema.object_xref (analysis_id);


--
-- Name: idx_16950_ensembl_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16950_ensembl_idx ON ensembl_core_schema.object_xref (ensembl_object_type, ensembl_id);


--
-- Name: idx_16950_xref_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16950_xref_idx ON ensembl_core_schema.object_xref (xref_id, ensembl_object_type, ensembl_id, analysis_id);


--
-- Name: idx_16954_object_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16954_object_idx ON ensembl_core_schema.ontology_xref (object_xref_id);


--
-- Name: idx_16954_object_source_type_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16954_object_source_type_idx ON ensembl_core_schema.ontology_xref (object_xref_id, source_xref_id, linkage_type);


--
-- Name: idx_16954_source_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16954_source_idx ON ensembl_core_schema.ontology_xref (source_xref_id);


--
-- Name: idx_16959_name_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16959_name_idx ON ensembl_core_schema.operon (display_label);


--
-- Name: idx_16959_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16959_seq_region_idx ON ensembl_core_schema.operon (seq_region_id, seq_region_start);


--
-- Name: idx_16959_stable_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16959_stable_id_idx ON ensembl_core_schema.operon (stable_id, version);


--
-- Name: idx_16964_operon_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16964_operon_idx ON ensembl_core_schema.operon_transcript (operon_id);


--
-- Name: idx_16964_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16964_seq_region_idx ON ensembl_core_schema.operon_transcript (seq_region_id, seq_region_start);


--
-- Name: idx_16964_stable_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16964_stable_id_idx ON ensembl_core_schema.operon_transcript (stable_id, version);


--
-- Name: idx_16968_operon_transcript_gene_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16968_operon_transcript_gene_idx ON ensembl_core_schema.operon_transcript_gene (operon_transcript_id, gene_id);


--
-- Name: idx_16972_checksum; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16972_checksum ON ensembl_core_schema.peptide_archive (md5_checksum);


--
-- Name: idx_16979_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16979_seq_region_idx ON ensembl_core_schema.prediction_exon (seq_region_id, seq_region_start);


--
-- Name: idx_16979_transcript_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16979_transcript_idx ON ensembl_core_schema.prediction_exon (prediction_transcript_id);


--
-- Name: idx_16984_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16984_analysis_idx ON ensembl_core_schema.prediction_transcript (analysis_id);


--
-- Name: idx_16984_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16984_seq_region_idx ON ensembl_core_schema.prediction_transcript (seq_region_id, seq_region_start);


--
-- Name: idx_16989_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16989_analysis_idx ON ensembl_core_schema.protein_align_feature (analysis_id);


--
-- Name: idx_16989_external_db_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16989_external_db_idx ON ensembl_core_schema.protein_align_feature (external_db_id);


--
-- Name: idx_16989_hit_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16989_hit_idx ON ensembl_core_schema.protein_align_feature (hit_name);


--
-- Name: idx_16989_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16989_seq_region_idx ON ensembl_core_schema.protein_align_feature (seq_region_id, analysis_id, seq_region_start, score);


--
-- Name: idx_16989_seq_region_idx_2; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16989_seq_region_idx_2 ON ensembl_core_schema.protein_align_feature (seq_region_id, seq_region_start);


--
-- Name: idx_16998_aln_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_16998_aln_idx ON ensembl_core_schema.protein_feature (translation_id, hit_name, seq_start, seq_end, hit_start, hit_end, analysis_id);


--
-- Name: idx_16998_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16998_analysis_idx ON ensembl_core_schema.protein_feature (analysis_id);


--
-- Name: idx_16998_hitname_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16998_hitname_idx ON ensembl_core_schema.protein_feature (hit_name);


--
-- Name: idx_16998_translation_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_16998_translation_idx ON ensembl_core_schema.protein_feature (translation_id);


--
-- Name: idx_17005_class; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17005_class ON ensembl_core_schema.repeat_consensus (repeat_class);


--
-- Name: idx_17005_consensus; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17005_consensus ON ensembl_core_schema.repeat_consensus (repeat_consensus);


--
-- Name: idx_17005_name; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17005_name ON ensembl_core_schema.repeat_consensus (repeat_name);


--
-- Name: idx_17005_type; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17005_type ON ensembl_core_schema.repeat_consensus (repeat_type);


--
-- Name: idx_17012_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17012_analysis_idx ON ensembl_core_schema.repeat_feature (analysis_id);


--
-- Name: idx_17012_repeat_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17012_repeat_idx ON ensembl_core_schema.repeat_feature (repeat_consensus_id);


--
-- Name: idx_17012_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17012_seq_region_idx ON ensembl_core_schema.repeat_feature (seq_region_id, seq_region_start);


--
-- Name: idx_17039_cs_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17039_cs_idx ON ensembl_core_schema.seq_region (coord_system_id);


--
-- Name: idx_17039_name_cs_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17039_name_cs_idx ON ensembl_core_schema.seq_region (name, coord_system_id);


--
-- Name: idx_17043_region_attribx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17043_region_attribx ON ensembl_core_schema.seq_region_attrib (seq_region_id, attrib_type_id, value);


--
-- Name: idx_17043_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17043_seq_region_idx ON ensembl_core_schema.seq_region_attrib (seq_region_id);


--
-- Name: idx_17043_type_val_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17043_type_val_idx ON ensembl_core_schema.seq_region_attrib (attrib_type_id, value);


--
-- Name: idx_17043_val_only_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17043_val_only_idx ON ensembl_core_schema.seq_region_attrib (value);


--
-- Name: idx_17050_mapping_set_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17050_mapping_set_idx ON ensembl_core_schema.seq_region_mapping (mapping_set_id);


--
-- Name: idx_17050_seq_region_mapping_uindex; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17050_seq_region_mapping_uindex ON ensembl_core_schema.seq_region_mapping (external_seq_region_id, internal_seq_region_id, mapping_set_id);


--
-- Name: idx_17054_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17054_seq_region_idx ON ensembl_core_schema.seq_region_synonym (seq_region_id);


--
-- Name: idx_17054_syn_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17054_syn_idx ON ensembl_core_schema.seq_region_synonym (synonym, seq_region_id);


--
-- Name: idx_17059_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17059_analysis_idx ON ensembl_core_schema.simple_feature (analysis_id);


--
-- Name: idx_17059_hit_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17059_hit_idx ON ensembl_core_schema.simple_feature (display_label);


--
-- Name: idx_17059_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17059_seq_region_idx ON ensembl_core_schema.simple_feature (seq_region_id, seq_region_start);


--
-- Name: idx_17063_new_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17063_new_idx ON ensembl_core_schema.stable_id_event (new_stable_id);


--
-- Name: idx_17063_old_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17063_old_idx ON ensembl_core_schema.stable_id_event (old_stable_id);


--
-- Name: idx_17063_uni_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17063_uni_idx ON ensembl_core_schema.stable_id_event (mapping_session_id, old_stable_id, new_stable_id, type);


--
-- Name: idx_17068_all_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17068_all_idx ON ensembl_core_schema.supporting_feature (exon_id, feature_type, feature_id);


--
-- Name: idx_17068_feature_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17068_feature_idx ON ensembl_core_schema.supporting_feature (feature_type, feature_id);


--
-- Name: idx_17074_analysis_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17074_analysis_idx ON ensembl_core_schema.transcript (analysis_id);


--
-- Name: idx_17074_canonical_translation_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17074_canonical_translation_idx ON ensembl_core_schema.transcript (canonical_translation_id);


--
-- Name: idx_17074_gene_index; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17074_gene_index ON ensembl_core_schema.transcript (gene_id);


--
-- Name: idx_17074_seq_region_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17074_seq_region_idx ON ensembl_core_schema.transcript (seq_region_id, seq_region_start);


--
-- Name: idx_17074_stable_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17074_stable_id_idx ON ensembl_core_schema.transcript (stable_id, version);


--
-- Name: idx_17074_xref_id_index; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17074_xref_id_index ON ensembl_core_schema.transcript (display_xref_id);


--
-- Name: idx_17082_transcript_attribx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17082_transcript_attribx ON ensembl_core_schema.transcript_attrib (transcript_id, attrib_type_id, value);


--
-- Name: idx_17082_transcript_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17082_transcript_idx ON ensembl_core_schema.transcript_attrib (transcript_id);


--
-- Name: idx_17082_type_val_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17082_type_val_idx ON ensembl_core_schema.transcript_attrib (attrib_type_id, value);


--
-- Name: idx_17082_val_only_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17082_val_only_idx ON ensembl_core_schema.transcript_attrib (value);


--
-- Name: idx_17089_transcript_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17089_transcript_idx ON ensembl_core_schema.transcript_intron_supporting_evidence (transcript_id);


--
-- Name: idx_17092_all_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17092_all_idx ON ensembl_core_schema.transcript_supporting_feature (transcript_id, feature_type, feature_id);


--
-- Name: idx_17092_feature_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17092_feature_idx ON ensembl_core_schema.transcript_supporting_feature (feature_type, feature_id);


--
-- Name: idx_17098_stable_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17098_stable_id_idx ON ensembl_core_schema.translation (stable_id, version);


--
-- Name: idx_17098_transcript_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17098_transcript_idx ON ensembl_core_schema.translation (transcript_id);


--
-- Name: idx_17102_translation_attribx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17102_translation_attribx ON ensembl_core_schema.translation_attrib (translation_id, attrib_type_id, value);


--
-- Name: idx_17102_translation_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17102_translation_idx ON ensembl_core_schema.translation_attrib (translation_id);


--
-- Name: idx_17102_type_val_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17102_type_val_idx ON ensembl_core_schema.translation_attrib (attrib_type_id, value);


--
-- Name: idx_17102_val_only_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17102_val_only_idx ON ensembl_core_schema.translation_attrib (value);


--
-- Name: idx_17110_anal_exdb_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17110_anal_exdb_idx ON ensembl_core_schema.unmapped_object (analysis_id, external_db_id);


--
-- Name: idx_17110_ext_db_identifier_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17110_ext_db_identifier_idx ON ensembl_core_schema.unmapped_object (external_db_id, identifier);


--
-- Name: idx_17110_id_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17110_id_idx ON ensembl_core_schema.unmapped_object (identifier);


--
-- Name: idx_17110_unique_unmapped_obj_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17110_unique_unmapped_obj_idx ON ensembl_core_schema.unmapped_object (ensembl_id, ensembl_object_type, identifier, unmapped_reason_id, parent, external_db_id);


--
-- Name: idx_17126_display_index; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17126_display_index ON ensembl_core_schema.xref (display_label);


--
-- Name: idx_17126_id_index; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE UNIQUE INDEX idx_17126_id_index ON ensembl_core_schema.xref (dbprimary_acc, external_db_id, info_type, info_text, version);


--
-- Name: idx_17126_info_type_idx; Type: INDEX; Schema: ensembl_core_schema; Owner: stefano
--

CREATE INDEX idx_17126_info_type_idx ON ensembl_core_schema.xref (info_type);


