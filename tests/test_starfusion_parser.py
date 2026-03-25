"""Tests for StarFusionParser."""

from pathlib import Path

from fusql.parsers.starfusion import StarFusionParser


class TestStarFusionParser:
    """Tests for StarFusionParser.parse()."""

    def test_parse_returns_list(self, sample_starfusion_path):
        """Parser returns a list of StarFusionEntry records."""
        parser = StarFusionParser()
        result = parser.parse(sample_starfusion_path)
        assert isinstance(result, list)
        assert len(result) > 0

    def test_parse_extracts_fusion_name(self, sample_starfusion_path):
        """Parser extracts the FusionName field correctly."""
        parser = StarFusionParser()
        fusions = parser.parse(sample_starfusion_path)
        # All rows in fixture are FIP1L1--PDGFRA
        for fusion in fusions:
            assert fusion.fusion_name == "FIP1L1--PDGFRA"

    def test_parse_extracts_gene_symbols(self, sample_starfusion_path):
        """Parser strips ENSG suffix from gene fields."""
        parser = StarFusionParser()
        fusions = parser.parse(sample_starfusion_path)
        # LeftGene column is "FIP1L1^ENSG00000145216.16"
        # RightGene column is "PDGFRA^ENSG00000134853.12"
        first = fusions[0]
        assert first.gene1 == "FIP1L1"
        assert first.gene2 == "PDGFRA"

    def test_parse_extracts_breakpoints(self, sample_starfusion_path):
        """Parser extracts numeric breakpoint positions from chr:pos:strand."""
        parser = StarFusionParser()
        fusions = parser.parse(sample_starfusion_path)
        # First row: LeftBreakpoint="chr4:53414722:+", RightBreakpoint="chr4:54274885:+"
        first = fusions[0]
        assert first.left_breakpos == 53414722
        assert first.right_breakpos == 54274885

    def test_parse_read_counts(self, sample_starfusion_path):
        """Parser extracts junction and spanning read counts."""
        parser = StarFusionParser()
        fusions = parser.parse(sample_starfusion_path)
        # First row: JunctionReadCount=22747, SpanningFragCount=0
        first = fusions[0]
        assert first.junction_reads == 22747
        assert first.spanning_reads == 0

    def test_sample_id_extraction(self, sample_starfusion_path):
        """Parser extracts sample_id from file path."""
        parser = StarFusionParser()
        fusions = parser.parse(sample_starfusion_path)
        # All fusions from same file share the same sample_id
        sample_ids = {f.sample_id for f in fusions}
        assert len(sample_ids) == 1
        assert isinstance(list(sample_ids)[0], str)
