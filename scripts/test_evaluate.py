from .evaluate import *
import pytest
import pysam
from pathlib import Path

TEST_CASES = Path("scripts/test_cases")
TEST_VCF = TEST_CASES / "test.vcf"
TEST_PANEL = TEST_CASES / "test_panel.fa"
TEST_REF_SEQ = TEST_CASES / "test_reference.fa"
TEST_TMP_PANEL = "/tmp/deleteme.fa"
TEST_MAKE_PROBE_VCF = TEST_CASES / "test_make_probe.vcf"
TEST_QUERY_VCF = TEST_CASES / "test_query.vcf"
TEST_QUERY_REF = TEST_CASES / "test_query.fa"


def retrieve_entry_from_test_vcf(idx: int) -> pysam.VariantRecord:
    with pysam.VariantFile(TEST_VCF) as vcf:
        for i, record in enumerate(vcf):
            if i == idx:
                return record
    raise IndexError("You asked for an index that is beyond the number in the test VCF")


def retrieve_entry_from_test_query_vcf(idx: int) -> pysam.VariantRecord:
    with pysam.VariantFile(TEST_QUERY_VCF) as vcf:
        for i, record in enumerate(vcf):
            if i == idx:
                return record
    raise IndexError("You asked for an index that is beyond the number in the test VCF")


def create_sam_header(name: str, length: int) -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_text(
        f"@SQ	SN:{name}	LN:{length}\n@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -t 1 panel.fa -"
    )


@pytest.fixture()
def sequence():
    return Sequence("0123456789")


class TestSequence:
    def test_getLeftFlank_atLeftBorder_returnSequenceFromRightBorder(self, sequence):
        assert sequence.get_left_flank(0, 4) == "6789"

    def test_getLeftFlank_closeToLeftBorder_returnSequenceFromLeftAndRightBorder(
        self, sequence
    ):
        assert sequence.get_left_flank(2, 4) == "8901"

    def test_getLeftFlank_atRightBorder_returnSequenceToLeft(self, sequence):
        assert sequence.get_left_flank(9, 4) == "5678"

    def test_getRight_Flank_inMiddle_returnSequenceToRight(self, sequence):
        assert sequence.get_right_flank(3, 4) == "4567"

    def test_getRightFlank_atRightBorder_returnSequenceFromLeftBorder(self, sequence):
        assert sequence.get_right_flank(8, 4) == "9012"

    def test_getRightFlank_atLeftBorder_returnSequenceToRight(self, sequence):
        assert sequence.get_right_flank(0, 4) == "1234"

    def test_getProbe_inMiddle_returnSequenceWithoutCrossingBorders(self, sequence):
        assert sequence.get_probe((5, 6), 4) == "123456789"

    def test_getProbe_nearLeftBorder_returnSequenceWithCrossingLeftBorder(
        self, sequence
    ):
        assert sequence.get_probe((2, 3), 4) == "890123456"

    def test_getProbe_nearRightBorder_returnSequenceWithCrossingRightBorder(
        self, sequence
    ):
        assert sequence.get_probe((7, 8), 4) == "345678901"

    def test_getProbe_atLeftBorder_returnSequenceWithCrossingLeftBorder(self, sequence):
        assert sequence.get_probe((0, 2), 4) == "6789012345"

    def test_getProbe_nearWithFlankLenOne_returnSequenceWithOneCharEitherSide(
        self, sequence
    ):
        assert sequence.get_probe((8, 9), 1) == "789"

    def test_getProbe_nearWithFlankLenZero_returnSequenceWithNoCharEitherSide(
        self, sequence
    ):
        assert sequence.get_probe((4, 6), 0) == "45"


class TestBwa:
    def test_align_refNotIndexed_raiseIndexError(self):
        with pytest.raises(IndexError) as excinfo:
            bwa = BWA()
            query = ">test\nTACGACACAGTGACGACATAAC"
            bwa.align(query)

        assert "indexed by BWA" in str(excinfo.value)

    def test_align_validQuery_returnValidSam(self):
        bwa = BWA()
        bwa.index(TEST_PANEL)
        query = ">test\nGACGTTAAATGCAAAAATCGCACGTCTTGAGCAGGATATAAAA"
        _, sam = bwa.align(query)
        expected = "test\t0\tC15154T\t55\t60\t43M\t*\t0\t0\tGACGTTAAATGCAAAAATCGCACGTCTTGAGCAGGATATAAAA\t*\tNM:i:0\tMD:Z:43\tAS:i:43\tXS:i:0"
        actual = sam[0].to_string()

        assert len(sam) == 1
        assert actual == expected

    def test_align_validQuery_returnValidHeader(self):
        bwa = BWA()
        bwa.index(TEST_PANEL)
        query = ">test\nGACGTTAAATGCAAAAATCGCACGTCTTGAGCAGGATATAAAA"
        header, _ = bwa.align(query)
        expected_start = (
            "@SQ\tSN:C15154T\tLN:201\n@SQ\tSN:T16509G\tLN:201\n@PG\tID:bwa\tPN:bwa\tVN:"
        )
        expected_end = f"CL:bwa mem -t 1 {TEST_PANEL} -\n"
        actual = str(header)

        # have to do it this way as the bwa version number is in the header and may vary
        assert actual.startswith(expected_start)
        assert actual.endswith(expected_end)

    def test_getOptions_defaultOneThread(self):
        bwa = BWA()

        actual = bwa.get_options()
        expected = ["-t", "1"]

        assert actual == expected

    def test_getOptions_setThreeThreads_returnThreeThreads(self):
        bwa = BWA(3)

        actual = bwa.get_options()
        expected = ["-t", "3"]

        assert actual == expected


class TestReference:
    def test_makePanel_invalidVcf_raisesFileNotFoundError(self):
        vcf = Path("doesntexist.vcf")
        sequence = Path(TEST_REF_SEQ)
        panel = Path(TEST_TMP_PANEL)
        ref = Reference(vcf, sequence)

        with pytest.raises(FileNotFoundError) as excinfo:
            ref.make_panel(panel)

    def test_makePanel_invalidReferenceSequence_raisesFileNotFoundError(self):
        vcf = Path(TEST_MAKE_PROBE_VCF)
        sequence = Path("doesntexist.fa")
        panel = Path(TEST_TMP_PANEL)
        ref = Reference(vcf, sequence)

        with pytest.raises(OSError) as excinfo:
            ref.make_panel(panel)

    def test_makePanel_oneVariantInMiddle_returnOneProbe(self):
        flank_width = 2
        vcf = Path(TEST_MAKE_PROBE_VCF)
        sequence = Path(TEST_REF_SEQ)
        panel = Path(TEST_TMP_PANEL)
        ref = Reference(vcf, sequence)

        ref.make_panel(panel, flank_width)

        expected_seq = "GATAC"
        expected = f">G6T flank_width={flank_width}\n{expected_seq}\n"
        actual = panel.read_text()

        assert actual == expected


class TestQuery:
    def test_calculateProbeBoundariesForEntry_variantShorterThanMinLen_returnProbeOfMinLen(
        self
    ):
        query = Query(TEST_QUERY_VCF, TEST_PANEL)
        query._min_probe_length = 10
        variant = retrieve_entry_from_test_query_vcf(1)

        expected = [16, 25]
        actual = query.calculate_probe_boundaries_for_entry(variant)

        assert actual == expected

    def test_calculateProbeBoundariesForEntry_variantAtStartOfGene_returnZeroLenLeftProbe(
        self
    ):
        query = Query(TEST_QUERY_VCF, TEST_PANEL)
        query._min_probe_length = 10
        variant = retrieve_entry_from_test_query_vcf(0)

        expected = [0, 6]
        actual = query.calculate_probe_boundaries_for_entry(variant)

        assert actual == expected

    def test_calculateProbeBoundariesForEntry_variantLongerThanMinLen_returnZeroLenProbes(
        self
    ):
        query = Query(TEST_QUERY_VCF, TEST_PANEL)
        query._min_probe_length = 2
        variant = retrieve_entry_from_test_query_vcf(1)

        expected = [19, 22]
        actual = query.calculate_probe_boundaries_for_entry(variant)

        assert actual == expected

    def test_calculateProbeBoundariesForEntry_nonRefvariant_returnProbesRelativeToRef(
        self
    ):
        query = Query(TEST_QUERY_VCF, TEST_PANEL)
        query._min_probe_length = 10
        variant = retrieve_entry_from_test_query_vcf(1)

        expected = [16, 25]
        actual = query.calculate_probe_boundaries_for_entry(variant)

        assert actual == expected

    def test_createProbeForVariant_nameNotInSeen_returnEntry0(self):
        query = Query(TEST_QUERY_VCF, TEST_PANEL)
        query._min_probe_length = 2
        variant = retrieve_entry_from_test_query_vcf(0)

        expected = pysam.FastxRecord()
        expected.set_name(
            f"GC00000001_155_pos1_entry0_CONF{get_genotype_confidence(variant)}"
        )
        actual = query.create_probe_for_variant(variant)

        assert actual.name == expected.name

    def test_createProbeForVariant_nameInSeen_returnEntry1(self):
        query = Query(TEST_QUERY_VCF, TEST_PANEL)
        query._min_probe_length = 2
        query._probe_names.add("GC00000001_155_pos1")
        variant = retrieve_entry_from_test_query_vcf(0)

        expected = pysam.FastxRecord()
        expected.set_name(
            f"GC00000001_155_pos1_entry1_CONF{get_genotype_confidence(variant)}"
        )
        actual = query.create_probe_for_variant(variant)

        assert actual.name == expected.name

    @pytest.mark.xfail(reason="Changing implementation")
    def test_createProbesForGeneVariants_emptyVariants_returnEmptyProbes(self):
        query = Query(TEST_QUERY_VCF, TEST_QUERY_REF)

        expected = ""
        actual = query.create_probes_for_gene_variants(pysam.FastxRecord(), [])

        assert actual == expected

    @pytest.mark.xfail(reason="Changing implementation")
    def test_makeProbes_twoValidVariantsOneInvalid_returnTwoProbes(self):
        query = Query(TEST_QUERY_VCF, TEST_QUERY_REF)
        query._min_probe_length = 6

        expected = ">GC00000001_155_pos1_entry0_CONF262.757\nCTGG\n>GC00000001_155_pos20_entry0_CONF262.757\nCTTGGC\n"
        actual = query.make_probes()

        assert actual == expected


def test_isInvalidVcfEntry_withNoneGenotype_returnTrue():
    entry = retrieve_entry_from_test_vcf(0)
    assert is_invalid_vcf_entry(entry)


def test_isInvalidVcfEntry_withGenotype1_returnFalse():
    entry = retrieve_entry_from_test_vcf(1)
    assert not is_invalid_vcf_entry(entry)


def test_getGenotypeConfidence():
    entry = retrieve_entry_from_test_vcf(0)
    assert get_genotype_confidence(entry) == 262.757


def test_getGenotype_genotypeNone_returnNone():
    entry = retrieve_entry_from_test_vcf(0)
    assert get_genotype(entry) is None


def test_getGenotype_genotype1_return1():
    entry = retrieve_entry_from_test_vcf(1)
    assert get_genotype(entry) == 1


def test_getVariantSequence_genotypeNone_returnRef():
    entry = retrieve_entry_from_test_vcf(0)
    expected = "CTGCCCGTTGGC"
    actual = get_variant_sequence(entry)

    assert actual == expected


def test_getVariantSequence_genotypeOne_returnFirstAlt():
    entry = retrieve_entry_from_test_vcf(1)
    expected = "TTGGGGGAAGGCTCTGCACTGCCCGTTGGC"
    actual = get_variant_sequence(entry)

    assert actual == expected


def test_getVariantSequence_genotypeZero_returnRef():
    entry = retrieve_entry_from_test_vcf(2)
    expected = "CTGCCCGTTGGC"
    actual = get_variant_sequence(entry)

    assert actual == expected


def test_getVariantLength_genotypeNone_returnRef():
    entry = retrieve_entry_from_test_vcf(0)
    expected = 12
    actual = get_variant_length(entry)

    assert actual == expected


def test_getVariantLength_genotypeOne_returnFirstAlt():
    entry = retrieve_entry_from_test_vcf(1)
    expected = 30
    actual = get_variant_length(entry)

    assert actual == expected


def test_getVariantLength_genotypeZero_returnRef():
    entry = retrieve_entry_from_test_vcf(2)
    expected = 12
    actual = get_variant_length(entry)

    assert actual == expected


def test_isMappingInvalid_unmappedEntry_returnTrue():
    header = create_sam_header("C15154T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00004785_pos200_entry0	4	C15154T	124	48	69M	*	0	0	CAAATCGGAAGCTAACAGAGCCAATACGCGCCTTGACGCCCAGGACTATTTTGATTGCCTGCGCTGCTT	*	NM:i:0	MD:Z:69	AS:i:69	XS:i:53",
        header,
    )

    assert is_mapping_invalid(record)


def test_isMappingInvalid_mappedEntry_returnFalse():
    header = create_sam_header("C15154T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00004785_pos200_entry0	0	C15154T	124	48	69M	*	0	0	CAAATCGGAAGCTAACAGAGCCAATACGCGCCTTGACGCCCAGGACTATTTTGATTGCCTGCGCTGCTT	*	NM:i:0	MD:Z:69	AS:i:69	XS:i:53",
        header,
    )
    is_mapping_valid = not is_mapping_invalid(record)

    assert is_mapping_valid


def test_isMappingInvalid_supplementaryEntry_returnTrue():
    header = create_sam_header("C15154T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00004785_pos200_entry0	2048	C15154T	124	48	69M	*	0	0	CAAATCGGAAGCTAACAGAGCCAATACGCGCCTTGACGCCCAGGACTATTTTGATTGCCTGCGCTGCTT	*	NM:i:0	MD:Z:69	AS:i:69	XS:i:53",
        header,
    )

    assert is_mapping_invalid(record)


def test_isMappingInvalid_secondaryEntry_returnTrue():
    header = create_sam_header("C15154T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00004785_pos200_entry0	256	C15154T	124	48	69M	*	0	0	CAAATCGGAAGCTAACAGAGCCAATACGCGCCTTGACGCCCAGGACTATTTTGATTGCCTGCGCTGCTT	*	NM:i:0	MD:Z:69	AS:i:69	XS:i:53",
        header,
    )

    assert is_mapping_invalid(record)


def test_isSnpCalledCorrectly_correctBaseInQuery_returnTrue():
    header = create_sam_header("C15154T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00004785_pos168_entry0	0	C15154T	94	0	70M	*	0	0	GAAAGCATTACGCCCACAAATCTATGCTGCCAAATCGGAAGCTAACAGAGCCAATACGCGCCTTGACGCC	*	NM:i:1	MD:Z:0A69	AS:i:69	XS:i:69	XA:Z:C15129A,+119,70M,1;",
        header,
    )

    assert is_snp_called_correctly(record)


def test_isSnpCalledCorrectly_incorrectBaseInQueryByChangingQuery_returnFalse():
    header = create_sam_header("C15154T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00004785_pos168_entry0	0	C15154T	94	0	70M	*	0	0	GAAAGCAATACGCCCACAAATCTATGCTGCCAAATCGGAAGCTAACAGAGCCAATACGCGCCTTGACGCC	*	NM:i:2	MD:Z:0A6T62	AS:i:69	XS:i:69	XA:Z:C15129A,+119,70M,1;",
        header,
    )

    assert not is_snp_called_correctly(record)


def test_isSnpCalledCorrectly_recordWithDeletionsButIncorrectBase_returnFalse():
    header = create_sam_header("G27981T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00002440_3_pos672_entry0_CONF60.7373	0	G27981T	31	60	33M3D36M	*	0	0	GTAATTGGTTCTTAGAGCAGCATATTTATGATCTTTTTGCTATTGTGGATGAACGACACTATAAGCAGT	*	NM:i:4	MD:Z:33^TTT34T1	AS:i:58	XS:i:0",
        header,
    )

    assert not is_snp_called_correctly(record)


def test_isSnpCalledCorrectly_recordWithDeletionsButCorrectBase_returnTrue():
    header = create_sam_header("G27981T", 201)
    record = pysam.AlignedSegment.fromstring(
        "GC00002440_3_pos672_entry0_CONF60.7373	0	G27981T	31	60	33M3D36M	*	0	0	GTAATTGGTTCTTAGAGCAGCATATTTATGATCTTTTTGCTATTGTGGATGAACGACACTATAAGCATT	*	NM:i:3	MD:Z:33^TTT36	AS:i:58	XS:i:0",
        header,
    )

    assert is_snp_called_correctly(record)


def test_mapProbesToPanel_oneRecordPerfectMapping():
    probe = pysam.FastxRecord()
    name = "GC00004785_pos168_entry0_CONF123.45"
    probe.set_name(name)
    probe.set_sequence("ACGTCTTGAGCAGGATATAAAAGCATTACGCCCACAAATCTATGCTGCCA")
    panel = Path(TEST_PANEL)

    actual = map_probes_to_panel(str(probe), panel)
    expected = {
        "snps_called_correctly": [True],
        "mismatches": [0],
        "ids": [name],
        "ref_ids": ["C15154T"],
        "total_pandora_calls": 1,
        "pandora_calls_crossing_ref_site": 1,
        "reference_sites_called": 1,
    }

    assert actual == expected


def test_mapProbesToPanel_oneRecordSnpCalledOneMismatch():
    probe = pysam.FastxRecord()
    name = "GC00004785_pos168_entry0_CONF123.45"
    probe.set_name(name)
    probe.set_sequence("ACGTCGTGAGCAGGATATAAAAGCATTACGCCCACAAATCTATGCTGCCA")
    panel = Path(TEST_PANEL)

    actual = map_probes_to_panel(str(probe), panel)
    expected = {
        "snps_called_correctly": [True],
        "mismatches": [1],
        "ids": [name],
        "ref_ids": ["C15154T"],
        "total_pandora_calls": 1,
        "pandora_calls_crossing_ref_site": 1,
        "reference_sites_called": 1,
    }

    assert actual == expected


def test_mapProbesToPanel_oneRecordSnpCalledTwoMismatches():
    probe = pysam.FastxRecord()
    name = "GC00004785_pos168_entry0_CONF123.45"
    probe.set_name(name)
    probe.set_sequence("ACGTCGTGAGCAGGATATAAAAGCATTACGCCCACAAATCTATGCTCCCA")
    panel = Path(TEST_PANEL)

    actual = map_probes_to_panel(str(probe), panel)
    expected = {
        "snps_called_correctly": [True],
        "mismatches": [2],
        "ids": [name],
        "ref_ids": ["C15154T"],
        "total_pandora_calls": 1,
        "pandora_calls_crossing_ref_site": 1,
        "reference_sites_called": 1,
    }

    assert actual == expected


def test_mapProbesToPanel_oneRecordSnpNotCalledNoOtherMismatches():
    probe = pysam.FastxRecord()
    name = "GC00004785_pos168_entry0_CONF123.45"
    probe.set_name(name)
    probe.set_sequence("TTAACGCCCTCAATTTTGAGGACGTAACCTACACCGACGCCAAAGGCGAA")
    panel = Path(TEST_PANEL)

    actual = map_probes_to_panel(str(probe), panel)
    expected = {
        "snps_called_correctly": [False],
        "mismatches": [1],
        "ids": [name],
        "ref_ids": ["T16509G"],
        "total_pandora_calls": 1,
        "pandora_calls_crossing_ref_site": 1,
        "reference_sites_called": 1,
    }

    assert actual == expected


def test_mapProbesToPanel_oneRecordDoesntMap():
    probe = pysam.FastxRecord()
    name = "GC00004785_pos168_entry0_CONF123.45"
    probe.set_name(name)
    probe.set_sequence("T" * 60)
    panel = Path(TEST_PANEL)

    actual = map_probes_to_panel(str(probe), panel)
    expected = {
        "snps_called_correctly": [],
        "mismatches": [],
        "ids": [],
        "ref_ids": [],
        "total_pandora_calls": 1,
        "pandora_calls_crossing_ref_site": 0,
        "reference_sites_called": 0,
    }

    assert actual == expected


def test_mapProbesToPanel_oneRecordMapsToPanelButToLeftOfVariant():
    probe = pysam.FastxRecord()
    name = "GC00004785_pos168_entry0_CONF123.45"
    probe.set_name(name)
    probe.set_sequence(
        "ACACGATAACGTATTACGTGCCATTGCAAATATTGAATGCTCAGAGAAATTTAACGCCCTCAATTTTGAGGACGT"
    )
    panel = Path(TEST_PANEL)

    actual = map_probes_to_panel(str(probe), panel)
    expected = {
        "snps_called_correctly": [],
        "mismatches": [],
        "ids": [],
        "ref_ids": [],
        "total_pandora_calls": 1,
        "pandora_calls_crossing_ref_site": 0,
        "reference_sites_called": 0,
    }

    assert actual == expected


def test_mapProbesToPanel_oneRecordMapsToPanelButToRightOfVariant():
    probe = pysam.FastxRecord()
    name = "GC00004785_pos168_entry0_CONF123.45"
    probe.set_name(name)
    probe.set_sequence(
        "ACCTACACCGACGCCAAAGGCGAAAAACGCCCAATGTACCAAATCACCAAAAACGGCTTCGTCTTCCTGGTGATGGGATTCACT"
    )
    panel = Path(TEST_PANEL)

    actual = map_probes_to_panel(str(probe), panel)
    expected = {
        "snps_called_correctly": [],
        "mismatches": [],
        "ids": [],
        "ref_ids": [],
        "total_pandora_calls": 1,
        "pandora_calls_crossing_ref_site": 0,
        "reference_sites_called": 0,
    }

    assert actual == expected


def test_mapProbesToPanel_oneRecordLastBaseIsVariantSite():
    probe = pysam.FastxRecord()
    name = "GC00004785_pos168_entry0_CONF123.45"
    probe.set_name(name)
    probe.set_sequence("TATTGAATGCTCAGAGAAATTTAACGCCCTCAATTTTGAGGACGTG")
    panel = Path(TEST_PANEL)

    actual = map_probes_to_panel(str(probe), panel)
    expected = {
        "snps_called_correctly": [True],
        "mismatches": [0],
        "ids": [name],
        "ref_ids": ["T16509G"],
        "total_pandora_calls": 1,
        "pandora_calls_crossing_ref_site": 1,
        "reference_sites_called": 1,
    }

    assert actual == expected


def test_mapProbesToPanel_oneRecordFirstBaseIsVariantSite():
    probe = pysam.FastxRecord()
    name = "GC00004785_pos168_entry0_CONF123.45"
    probe.set_name(name)
    probe.set_sequence("GACCTACACCGACGCCAAAGGCGAAAAACGCCCAATGTACCAAATCACCAAAAAC")
    panel = Path(TEST_PANEL)

    actual = map_probes_to_panel(str(probe), panel)
    expected = {
        "snps_called_correctly": [True],
        "mismatches": [0],
        "ids": [name],
        "ref_ids": ["T16509G"],
        "total_pandora_calls": 1,
        "pandora_calls_crossing_ref_site": 1,
        "reference_sites_called": 1,
    }

    assert actual == expected


def test_writeResults():
    data = {"a": [1, 2, 3], "foo": "bar", "bool": [True, False]}
    output = Path("/tmp/output.json")

    write_results(data, output)

    expected = '{\n    "a": [\n        1,\n        2,\n        3\n    ],\n    "foo": "bar",\n    "bool": [\n        true,\n        false\n    ]\n}'
    actual = output.read_text()

    assert actual == expected


def test_candidateContainsExpectedSnp_correctSnpCalled_returnTrue():
    header = create_sam_header("GC00001889_1.0", 978)
    record = pysam.AlignedSegment.fromstring(
        "T36234G	0	GC00001889_1.0	710	60	201M	*	0	0	TTGCCGGTCTGTTTATTGCTATCGGTCACAGCCCGAATACTGCGATTTTCGGAGGGCAGCTGGAACTGGAAAGCGGCTACATCAAAGTACAGTCGGGCATGCATGGTAATGCCACCCAGACCAGCATCCCTGGCGTCTTTGCCGCAGGCGACGTGATGGATCACATTTATCGCCAGGCTATTACCTCGGCCGGTACAGGCT	*	NM:i:1	MD:Z:25C175	AS:i:196	XS:i:0",
        header,
    )

    assert candidate_contains_expected_snp(record)


def test_candidateContainsExpectedSnp_incorrectSnpCalled_returnFalse():
    header = create_sam_header("GC00000001_27.0", 1824)
    record = pysam.AlignedSegment.fromstring(
        "T25983C	0	GC00000001_27.0	440	60	201M	*	0	0	CCATGGTGATGGCCACAGTGGATTACGTGGATGACAAACTGAAAGAGCATGAACAGTCACGACGTCACCCGGATGCCTCGCTGACCACAAAAGGTTTTACCCAGTTAACCAGCGACACCACCAGTACATCAGAGACCCGGGCAGCCACGCCAAAGGCCATAAAGATCGCCATGGACAACGCTGATGGTCGGCTTGCGAAAA	*	NM:i:2	MD:Z:100T7G92	AS:i:191	XS:i:0",
        header,
    )

    assert not candidate_contains_expected_snp(record)


def test_candidateContainsExpectedSnp_variantInSoftClippedRegion_returnFalse():
    header = create_sam_header("GC00000001_27.0", 1824)
    record = pysam.AlignedSegment.fromstring(
        "A25435C	0	GC00000001_27.0	1	60	109S92M	*	0	0	CGGTAAAAGTAAAGTCTATAGCTACTTCTGGTTTTTTATCTGGTTTTTTGTTTTTAACTGTGCTGTTTATTCTTCTCTGAAGACAATAAATATAAAGGATCTAATCAGCATGAGTACCACAACACGAAAATTTAAAACCGTTATCACCGATACGGGTGCAAAAAAATTAGCTCAGGCAGCCGCGCCAGATGGTAAGCCTGT	*	NM:i:0	MD:Z:92	AS:i:92	XS:i:0",
        header,
    )

    assert not candidate_contains_expected_snp(record)


def test_candidateContainsExpectedSnp_variantPositionNotInAlignment_returnFalse():
    header = create_sam_header("GC00000001_155.0", 1824)
    record = pysam.AlignedSegment.fromstring(
        "G6912T	2048	GC00000001_155.0	1698	0	112H89M	*	0	0	TCAGCTCAGTAGCGACACCAACAGCACTTCTGAAACATTAGCTGCAACTCCAAAAGCGGTTAAAGCTGCATACGATCTAGCAGCCGGAA	*	NM:i:1	MD:Z:48G40	AS:i:84	XS:i:84	SA:Z:GC00000001_155.0,1754,+,146M55S,3,6;	XA:Z:GC00000001_155.0,+1098,112S89M,1;GC00000001_155.0,+1398,112S89M,1;GC00000001_155.0,+798,112S89M,1;",
        header,
    )

    assert not candidate_contains_expected_snp(record)


def test_NoOverlappingIntervals_NoChange():
    intervals = [[2, 4], [6, 9], [11, 12]]

    result = merge_overlap_intervals(intervals)
    expected = [[2, 4], [6, 9], [11, 12]]
    assert result == expected


def test_TwoIntervalsEqualEndStart_NoChange():
    intervals = [[6, 9], [9, 12]]

    result = merge_overlap_intervals(intervals)
    expected = [[6, 9], [9, 12]]
    assert result == expected


def test_TwoIntervalsOverlap_Merge():
    intervals = [[6, 9], [8, 12]]

    result = merge_overlap_intervals(intervals)
    expected = [[6, 12]]
    assert result == expected


def test_ThreeIntervalsOverlap_Merge():
    intervals = [[6, 9], [8, 12], [11, 14]]

    result = merge_overlap_intervals(intervals)
    expected = [[6, 14]]
    assert result == expected


def test_ThreeIntervalsOverlapTwoEqualsEndStart_MergeOverlapDontMergeEquals():
    intervals = [[6, 9], [8, 12], [11, 14], [14, 16]]

    result = merge_overlap_intervals(intervals)
    expected = [[6, 14], [14, 16]]

    assert result == expected


def test_FindIndexInIntervals_emptyIntervalsReturnsNegativeOne():
    intervals = []
    query = 2

    actual = find_index_in_intervals(intervals, query)
    expected = -1

    assert actual == expected


def test_FindIndexInIntervals_queryNotInIntervalsReturnsNegativeOne():
    intervals = [[3, 7], [9, 20]]
    query = 2

    actual = find_index_in_intervals(intervals, query)
    expected = -1

    assert actual == expected


def test_FindIndexInIntervals_queryInFirstIntervalsReturnsZero():
    intervals = [[3, 7], [9, 20]]
    query = 5

    actual = find_index_in_intervals(intervals, query)
    expected = 0

    assert actual == expected


def test_FindIndexInIntervals_queryInSecondIntervalsReturnsOne():
    intervals = [[3, 7], [9, 20]]
    query = 14

    actual = find_index_in_intervals(intervals, query)
    expected = 1

    assert actual == expected


def test_FindIndexInIntervals_queryEqualsStartOfFirstIntervalReturnsZero():
    intervals = [[3, 7], [9, 20]]
    query = 3

    actual = find_index_in_intervals(intervals, query)
    expected = 0

    assert actual == expected


def test_FindIndexInIntervals_queryEqualsEndOfFirstIntervalReturnsZero():
    intervals = [[3, 7], [9, 20]]
    query = 7

    actual = find_index_in_intervals(intervals, query)
    expected = 0

    assert actual == expected


def test_FindIndexInIntervals_queryGreaterthanLastIntervalReturnsNegativeOne():
    intervals = [[3, 7], [9, 20]]
    query = 70

    actual = find_index_in_intervals(intervals, query)
    expected = -1

    assert actual == expected


def test_AnyEntriesOverlap_emptyEntriesReturnsFalse():
    entries = []

    assert not any_entries_overlap(entries)


def test_AnyEntriesOverlap_oneEntryReturnsFalse():
    entries = [retrieve_entry_from_test_query_vcf(1)]

    assert not any_entries_overlap(entries)


def test_AnyEntriesOverlap_twoNonOverlappingEntriesReturnsFalse():
    entries = [
        retrieve_entry_from_test_query_vcf(1),
        retrieve_entry_from_test_query_vcf(2),
    ]

    assert not any_entries_overlap(entries)

def test_AnyEntriesOverlap_twoOverlappingEntriesReturnsTrue():
    entries = [
        retrieve_entry_from_test_vcf(1),
        retrieve_entry_from_test_query_vcf(1),
    ]

    assert any_entries_overlap(entries)