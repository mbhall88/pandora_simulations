from .results import *
from pathlib import Path
import pytest

TEST_CASES = Path("scripts/test_cases")
EMPTY_JSON = TEST_CASES / "empty.json"
TEST_JSON = TEST_CASES / "test.json"


class TestVariantCall:
    def test_equality_twoEqualVariantCalls_returnTrue(self):
        var1 = VariantCall("name_pos1_entry0_CONF8.8", True, 0, "A44C")
        var2 = VariantCall("name_pos1_entry0_CONF8.8", True, 7, "A44C")

        assert var1 == var2

    def test_equality_twoNonEqualVariantCalls_returnFalse(self):
        var1 = VariantCall("name_pos1_entry0_CONF8.8", True, 0, "A44C")
        var2 = VariantCall("name_pos1_entry0_CONF8.7", True, 0, "A44C")

        assert not var1 == var2

    def test_nonEquality_twoEqualVariantCalls_returnFalse(self):
        var1 = VariantCall("name_pos1_entry0_CONF8.8", True, 0, "A44C")
        var2 = VariantCall("name_pos1_entry0_CONF8.8", True, 0, "A44C")

        assert not var1 != var2

    def test_nonEquality_twoNonEqualVariantCalls_returnTrue(self):
        var1 = VariantCall("name_pos1_entry0_CONF8.8", True, 0, "A44C")
        var2 = VariantCall("name_pos1_entry0_CONF8.7", True, 0, "A44C")

        assert var1 != var2


class TestResult:
    def test_equality_twoEqualResults_returnTrue(self):
        result1 = Result(coverage=55, num_snps=8)
        result1.data = dict(test="foobar")
        result2 = Result(coverage=55, num_snps=8)
        result2.data = dict(test="foobar")

        assert result1 == result2

    def test_equality_twoNonEqualResults_returnFalse(self):
        result1 = Result(coverage=55, num_snps=8, read_quality="perfect")
        result1.data = dict(test="foobar")
        result2 = Result(coverage=55, num_snps=8)
        result2.data = dict(test="foobar")

        assert not result1 == result2

    def test_nonEquality_twoEqualResults_returnFalse(self):
        result1 = Result(coverage=55, num_snps=8)
        result1.data = dict(test="foobar")
        result2 = Result(coverage=55, num_snps=8)
        result2.data = dict(test="foobar")

        assert not result1 != result2

    def test_nonEquality_twoNonEqualResults_returnTrue(self):
        result1 = Result(coverage=55, num_snps=8, read_quality="perfect")
        result1.data = dict(test="foobar")
        result2 = Result(coverage=55, num_snps=8)
        result2.data = dict(test="foobar")

        assert result1 != result2

    def test_loadDataFromJson_emptyJsonReturnsEmptyData(self):
        result = Result()
        result.load_data_from_json(EMPTY_JSON)

        expected = dict()
        actual = result.data

        assert actual == expected

    def test_loadDataFromJson_jsonWithContentsReturnsDictWithSameContents(self):
        result = Result()
        result.load_data_from_json(TEST_JSON)

        expected = {
            "reference_sites_called": 7,
            "ids": ["id1_pos1_entry0_CONF3", "id2_pos2_entry1_CONF5.6"],
            "snps_called_correctly": [True, False],
            "mismatches": [0, 3],
            "ref_ids": ["A44C", "A44C"],
        }
        actual = result.data

        assert actual == expected

    def test_uniqueSnpsCalled_noData_returnsZero(self):
        result = Result()

        actual = result.unique_snps_called()
        expected = 0

        assert actual == expected

    def test_uniqueSnpsCalled_withSevenUniqueSnps_returnsSeven(self):
        result = Result()
        result.load_data_from_json(TEST_JSON)

        actual = result.unique_snps_called()
        expected = 7

        assert actual == expected

    def test_getVariantCalls_noData_returnEmptyList(self):
        result = Result()

        actual = result._get_variant_calls()
        expected = []

        assert actual == expected

    def test_getVariantCalls_twoVariantCalls_returnListWithTwoCalls(self):
        result = Result()
        result.load_data_from_json(TEST_JSON)

        actual = result.variant_calls
        expected = [
            VariantCall("id1_pos1_entry0_CONF3", True, 0, "A44C"),
            VariantCall("id2_pos2_entry1_CONF5.6", False, 3, "A44C"),
        ]

        assert actual == expected

    def test_extractParametersFromPath_invalidPath_raiseIndexError(self):
        path = Path("invalid/path")

        with pytest.raises(IndexError) as excinfo:
            actual = Result.extract_parameters_from_path(path)

        assert "index out of range" in str(excinfo.value)

    def test_extractParametersFromPath_validPath_returnValidResult(self):
        path = Path("analysis/1/6/perfect/50/11/evaluate/result.json")

        actual = Result.extract_parameters_from_path(path)
        expected = dict(
            denovo_kmer_size=11,
            coverage=50,
            read_quality="perfect",
            num_snps=6,
            max_nesting=1,
        )

        assert actual == expected

    def test_fromDict_emptyDict_returnEmptyResult(self):
        actual = Result.from_dict(dict())
        expected = Result()

        assert actual == expected

    def test_fromDict_dictWithAllParams_returnResultWithSameParams(self):
        actual = Result.from_dict(
            dict(
                denovo_kmer_size=11,
                coverage=50,
                read_quality="perfect",
                num_snps=6,
                max_nesting=1,
            )
        )
        expected = Result(11, 50, "perfect", 6, 1)

        assert actual == expected

    def test_denovoRecall_noVariantSitesCalledCorrectly_returnZero(self):
        result = Result.from_dict(dict(num_snps=50))

        actual = result.denovo_recall()
        expected = 0.0

        assert actual == expected

    def test_denovoRecall_numSnpsIsZero_returnZero(self):
        result = Result.from_dict(dict(num_snps=0))
        result.data = dict(variant_sites_denovo_correctly_discovered=4)

        actual = result.denovo_recall()
        expected = 0.0

        assert actual == expected

    def test_denovoRecall_halfSnpsCalledCorrectly_returnHalf(self):
        result = Result.from_dict(dict(num_snps=50))
        result.data = dict(variant_sites_denovo_correctly_discovered=25)

        actual = result.denovo_recall()
        expected = 0.5

        assert actual == expected

    def test_denovoPrecision_allZeroes_returnZero(self):
        result = Result()

        actual = result.denovo_precision()
        expected = 0.0

        assert actual == expected

    def test_denovoPrecision_halfSlicesWIthCorectSnp_returnHalf(self):
        result = Result()
        result.data = dict(
            variant_sites_denovo_correctly_discovered=25, total_slices=50
        )

        actual = result.denovo_precision()
        expected = 0.5

        assert actual == expected

    def test_overallRecall_noFalseNegatives_returnOne(self):
        result = Result.from_dict(dict(num_snps=2))
        result.data = dict(
            reference_sites_called=2,
            ids=[
                "id1_pos1_entry0_CONF3",
                "id2_pos1_entry0_CONF3",
                "id3_pos5_entry0_CONF3",
            ],
            snps_called_correctly=[True, False, True],
            mismatches=[0, 2, 0],
            ref_ids=["A1C", "T1C", "T5G"],
        )
        result.variant_calls = result._get_variant_calls()

        actual = result.overall_recall()
        expected = 1.0

        assert actual == expected

    def test_overallRecall_noTruePositives_returnZero(self):
        result = Result.from_dict(dict(num_snps=5))
        result.data = dict(
            reference_sites_called=0,
            ids=["id1_pos1_entry0_CONF3", "id2_pos1_entry0_CONF3"],
            snps_called_correctly=[False, False],
            mismatches=[1, 2],
        )
        result.variant_calls = result._get_variant_calls()

        actual = result.overall_recall()
        expected = 0.0

        assert actual == expected

    def test_overallRecall_allZeroes_returnZero(self):
        result = Result.from_dict(dict(num_snps=0))

        actual = result.overall_recall()
        expected = 0.0

        assert actual == expected

    def test_overallRecall_halfTrueVariants_returnHalf(self):
        result = Result.from_dict(dict(num_snps=4))
        result.data = dict(
            reference_sites_called=2,
            ids=["id1_pos1_entry0_CONF3", "id2_pos2_entry0_CONF3"],
            snps_called_correctly=[True, True],
            mismatches=[1, 2],
            ref_ids=["A1T", "C2T"],
        )
        result.variant_calls = result._get_variant_calls()

        actual = result.overall_recall()
        expected = 0.5

        assert actual == expected

    def test_overallPrecision_halfTrueVariants_returnHalf(self):
        result = Result.from_dict(dict(num_snps=4))
        result.data = dict(
            reference_sites_called=2,
            ids=["id1_pos1_entry0_CONF3", "id2_pos1_entry0_CONF3"],
            snps_called_correctly=[True, False],
            mismatches=[1, 2],
            ref_ids=["A1T", "T2G"],
        )
        result.variant_calls = result._get_variant_calls()

        actual = result.overall_precision()
        expected = 0.5

        assert actual == expected

    def test_overallPrecision_noFalsePositives_returnOne(self):
        result = Result.from_dict(dict(num_snps=5))
        result.data = dict(
            reference_sites_called=2,
            ids=["id1_pos1_entry0_CONF3", "id2_pos1_entry0_CONF3"],
            snps_called_correctly=[True, True],
            mismatches=[0, 2],
            ref_ids=["A1T", "C5T"],
        )
        result.variant_calls = result._get_variant_calls()

        actual = result.overall_precision()
        expected = 1.0

        assert actual == expected

    def test_overallPrecision_noTruePositives_returnZero(self):
        result = Result.from_dict(dict(num_snps=5))
        result.data = dict(
            reference_sites_called=0,
            ids=["id1_pos1_entry0_CONF3", "id2_pos1_entry0_CONF3"],
            snps_called_correctly=[False, False],
            mismatches=[1, 2],
        )
        result.variant_calls = result._get_variant_calls()

        actual = result.overall_precision()
        expected = 0.0

        assert actual == expected

    def test_overallPrecision_allZeroes_returnZero(self):
        result = Result.from_dict(dict(num_snps=0))

        actual = result.overall_precision()
        expected = 0.0

        assert actual == expected
