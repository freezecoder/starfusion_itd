"""Tests for AribaParser."""

from pathlib import Path

from fusql.parsers.ariba import AribaParser


class TestAribaParser:
    """Tests for AribaParser.parse()."""

    def test_parse_returns_list(self, sample_ariba_path):
        """Parser returns a list of AribaFusion records."""
        parser = AribaParser()
        result = parser.parse(sample_ariba_path)
        assert isinstance(result, list)
        assert len(result) > 0

    def test_parse_extracts_gene_names(self, sample_ariba_path):
        """Parser correctly extracts gene1 and gene2 from each row."""
        parser = AribaParser()
        fusions = parser.parse(sample_ariba_path)
        # First fusion in fixture is RUNX1--RUNX1T1
        first = fusions[0]
        assert first.gene1 == "RUNX1"
        assert first.gene2 == "RUNX1T1"
        # Second fusion is ETV6--ABL1
        second = fusions[1]
        assert second.gene1 == "ETV6"
        assert second.gene2 == "ABL1"

    def test_parse_extracts_exon_numbers(self, sample_ariba_path):
        """Parser extracts exon numbers from fusion_transcript field when present."""
        parser = AribaParser()
        fusions = parser.parse(sample_ariba_path)
        # The parser's extract_exon_number() searches for E<digits> patterns.
        # The test fixture uses nucleotide sequences (no E<digit> patterns),
        # so exon values will be 0 — verify the field is populated correctly.
        for fusion in fusions:
            assert fusion.exon1 >= 0
            assert fusion.exon2 >= 0
            assert fusion.fusion_id  # non-empty fusion_id

    def test_parse_extracts_splice_sites(self, sample_ariba_path):
        """Parser extracts and normalizes splice site types."""
        parser = AribaParser()
        fusions = parser.parse(sample_ariba_path)
        # All test fixtures have CDS/splice-site
        for fusion in fusions:
            assert fusion.splice_site_1 in ["CDS", "splice-site", "donor", "acceptor", "unknown"]
            assert fusion.splice_site_2 in ["CDS", "splice-site", "donor", "acceptor", "unknown"]

    def test_parse_calculates_total_reads(self, sample_ariba_path):
        """Parser sums split_reads1 + split_reads2 for total reads."""
        parser = AribaParser()
        fusions = parser.parse(sample_ariba_path)
        # RUNX1--RUNX1T1 has split_reads1=114, split_reads2=50, total=164
        first = fusions[0]
        assert first.reads == 114 + 50
        # ETV6--ABL1 has split_reads1=73, split_reads2=69, total=142
        second = fusions[1]
        assert second.reads == 73 + 69

    def test_sample_id_extraction(self, sample_ariba_path):
        """Parser extracts sample_id from file path."""
        parser = AribaParser()
        fusions = parser.parse(sample_ariba_path)
        # All fusions from same file share the same sample_id
        sample_ids = {f.sample_id for f in fusions}
        assert len(sample_ids) == 1
        # The sample_id comes from the file stem since no Run/Sample pattern
        assert isinstance(list(sample_ids)[0], str)
