#!/usr/bin/env python3
# coding=utf-8

import sys
import os
import tempfile
import pytest
from unittest.mock import patch, MagicMock, call

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin')))

# Patch EnaApiHandler before importing so the module-level instantiation doesn't fail
with patch('ena_portal_api.ena_handler.EnaApiHandler'):
    import generate_inputs


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_assembly(accession):
    return {'analysis_accession': accession, 'generated_ftp': ''}


def make_xml(*, run_ref=None, assembly_name=None, filename=None):
    """Build the minimal dict structure returned by load_xml."""
    analysis = {}

    if run_ref is not None:
        analysis['RUN_REF'] = {'@accession': run_ref}

    seq_assembly = {}
    if assembly_name is not None:
        seq_assembly['NAME'] = assembly_name
    if seq_assembly:
        analysis['ANALYSIS_TYPE'] = {'SEQUENCE_ASSEMBLY': seq_assembly}

    if filename is not None:
        analysis['FILES'] = {'FILE': {'@filename': filename}}

    return {'ANALYSIS_SET': {'ANALYSIS': analysis}}


# ---------------------------------------------------------------------------
# get_run_from_assembly_xml
# ---------------------------------------------------------------------------

class TestGetRunFromAssemblyXml:

    def test_run_ref_single_dict(self):
        """Path 1: RUN_REF is a single dict with @accession."""
        xml = make_xml(run_ref='ERR3063408')
        with patch.object(generate_inputs, 'load_xml', return_value=xml):
            result = generate_inputs.get_run_from_assembly_xml(make_assembly('ERZ1234567'))
        assert result == 'ERR3063408'

    def test_run_ref_list(self):
        """Path 1: RUN_REF is a list; first element is used."""
        xml = make_xml()
        xml['ANALYSIS_SET']['ANALYSIS']['RUN_REF'] = [
            {'@accession': 'ERR3063408'},
            {'@accession': 'ERR9999999'},
        ]
        with patch.object(generate_inputs, 'load_xml', return_value=xml):
            result = generate_inputs.get_run_from_assembly_xml(make_assembly('ERZ1234567'))
        assert result == 'ERR3063408'

    def test_fallback_to_assembly_name(self):
        """Path 2: no RUN_REF; run accession embedded in SEQUENCE_ASSEMBLY/NAME."""
        xml = make_xml(assembly_name='ERR3063408_METASPADES')
        with patch.object(generate_inputs, 'load_xml', return_value=xml):
            result = generate_inputs.get_run_from_assembly_xml(make_assembly('ERZ1234567'))
        assert result == 'ERR3063408'

    def test_fallback_to_filename(self):
        """Path 3: no RUN_REF, no NAME; run accession in FILE @filename."""
        xml = make_xml(filename='webin-data-mgmt/ERR3063408_assembly.fasta.gz')
        with patch.object(generate_inputs, 'load_xml', return_value=xml):
            result = generate_inputs.get_run_from_assembly_xml(make_assembly('ERZ1234567'))
        assert result == 'ERR3063408'

    def test_drr_accession(self):
        """DRR (DDBJ) accessions are accepted."""
        xml = make_xml(run_ref='DRR123456')
        with patch.object(generate_inputs, 'load_xml', return_value=xml):
            result = generate_inputs.get_run_from_assembly_xml(make_assembly('ERZ1234567'))
        assert result == 'DRR123456'

    def test_srr_accession(self):
        """SRR (NCBI SRA) accessions are accepted."""
        xml = make_xml(run_ref='SRR9876543')
        with patch.object(generate_inputs, 'load_xml', return_value=xml):
            result = generate_inputs.get_run_from_assembly_xml(make_assembly('ERZ1234567'))
        assert result == 'SRR9876543'

    def test_no_run_found_exits(self):
        """No run accession in any field → sys.exit(1)."""
        xml = make_xml(assembly_name='no_run_here', filename='no_run_here.fasta.gz')
        with patch.object(generate_inputs, 'load_xml', return_value=xml):
            with pytest.raises(SystemExit):
                generate_inputs.get_run_from_assembly_xml(make_assembly('ERZ1234567'))

    def test_load_xml_failure_returns_none(self):
        """load_xml returning None causes the function to return None early."""
        with patch.object(generate_inputs, 'load_xml', return_value=None):
            result = generate_inputs.get_run_from_assembly_xml(make_assembly('ERZ1234567'))
        assert result is None

    def test_run_ref_takes_priority_over_name(self):
        """RUN_REF is preferred even when NAME also contains a run accession."""
        xml = make_xml(run_ref='ERR1111111', assembly_name='ERR2222222_assembly')
        with patch.object(generate_inputs, 'load_xml', return_value=xml):
            result = generate_inputs.get_run_from_assembly_xml(make_assembly('ERZ1234567'))
        assert result == 'ERR1111111'

    def test_name_takes_priority_over_filename(self):
        """NAME match is preferred over filename when RUN_REF is absent."""
        xml = make_xml(assembly_name='ERR1111111_assembly', filename='ERR2222222.fasta.gz')
        with patch.object(generate_inputs, 'load_xml', return_value=xml):
            result = generate_inputs.get_run_from_assembly_xml(make_assembly('ERZ1234567'))
        assert result == 'ERR1111111'


# ---------------------------------------------------------------------------
# transform_paths
# ---------------------------------------------------------------------------

class TestTransformPaths:

    def test_public_ftp_becomes_https(self):
        path = 'ftp.sra.ebi.ac.uk/vol1/fastq/ERR306/ERR3063408_1.fastq.gz'
        result, privacy = generate_inputs.transform_paths(path)
        assert result.startswith('https://')
        assert privacy == 'public'

    def test_private_ftp_unchanged(self):
        path = 'ftp.dcc-private.ebi.ac.uk/vol1/fastq/ERR306/ERR3063408_1.fastq.gz'
        result, privacy = generate_inputs.transform_paths(path)
        assert result == path
        assert privacy == 'private'

    def test_invalid_path_raises(self):
        with pytest.raises(ValueError):
            generate_inputs.transform_paths('s3://some-bucket/file.fastq.gz')


# ---------------------------------------------------------------------------
# fetch_data — run-coverage checks
# ---------------------------------------------------------------------------

class TestFetchDataRunCoverage:
    """Tests for the check that all assemblies have a run detected."""

    RUNS = [
        {'run_accession': 'ERR1000001', 'library_source': 'METAGENOMIC',
         'library_strategy': 'WGS', 'fastq_ftp': 'ftp.sra.ebi.ac.uk/vol1/ERR1000001.fastq.gz'},
        {'run_accession': 'ERR1000002', 'library_source': 'METAGENOMIC',
         'library_strategy': 'WGS', 'fastq_ftp': 'ftp.sra.ebi.ac.uk/vol1/ERR1000002.fastq.gz'},
    ]
    ASSEMBLIES = [
        {'analysis_accession': 'ERZ0000001', 'generated_ftp': 'ftp.sra.ebi.ac.uk/vol1/ERZ0000001.fasta.gz'},
        {'analysis_accession': 'ERZ0000002', 'generated_ftp': 'ftp.sra.ebi.ac.uk/vol1/ERZ0000002.fasta.gz'},
        {'analysis_accession': 'ERZ0000003', 'generated_ftp': 'ftp.sra.ebi.ac.uk/vol1/ERZ0000003.fasta.gz'},
    ]

    def _run_fetch(self, run_map, tmpdir):
        """
        run_map: dict of analysis_accession -> run accession (or None).
        Patches handler + get_run_from_assembly_xml and calls fetch_data.
        """
        def fake_get_run(assembly):
            return run_map.get(assembly['analysis_accession'])

        with patch.object(generate_inputs.handler, 'get_study_runs', return_value=self.RUNS), \
             patch.object(generate_inputs.handler, 'get_study_assemblies', return_value=self.ASSEMBLIES), \
             patch.object(generate_inputs, 'get_run_from_assembly_xml', side_effect=fake_get_run):
            return generate_inputs.fetch_data(
                raw_reads_study='PRJEB31790',
                assembly_study='ERP190907',
                outdir=tmpdir,
                input_scientific_name=None,
                input_env_biome=None,
                keep_metat=False,
            )

    def test_all_resolved_no_warning(self, capsys):
        run_map = {
            'ERZ0000001': 'ERR1000001',
            'ERZ0000002': 'ERR1000002',
            'ERZ0000003': 'ERR1000001',
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            assembly_run, _, _ = self._run_fetch(run_map, tmpdir)

        out = capsys.readouterr().out
        assert 'WARNING' not in out
        assert 'Resolved runs for 3/3' in out
        assert len(assembly_run) == 3

    def test_some_missing_run_prints_warning(self, capsys):
        run_map = {
            'ERZ0000001': 'ERR1000001',
            'ERZ0000002': None,           # XML could not be parsed
            'ERZ0000003': 'ERR1000002',
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            assembly_run, _, _ = self._run_fetch(run_map, tmpdir)

        out = capsys.readouterr().out
        assert 'WARNING' in out
        assert 'ERZ0000002' in out
        assert 'Resolved runs for 2/3' in out
        assert 'ERZ0000002' not in assembly_run

    def test_all_missing_run_exits(self):
        run_map = {
            'ERZ0000001': None,
            'ERZ0000002': None,
            'ERZ0000003': None,
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(SystemExit):
                self._run_fetch(run_map, tmpdir)

    def test_missing_run_not_in_assembly_run_dict(self, capsys):
        """Assemblies with no detected run must not appear in the returned dict."""
        run_map = {
            'ERZ0000001': 'ERR1000001',
            'ERZ0000002': None,
            'ERZ0000003': None,
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            assembly_run, _, _ = self._run_fetch(run_map, tmpdir)

        assert 'ERZ0000001' in assembly_run
        assert 'ERZ0000002' not in assembly_run
        assert 'ERZ0000003' not in assembly_run


# ---------------------------------------------------------------------------
# Integration test — real ENA API call (skipped by default)
# ---------------------------------------------------------------------------

class TestIntegration:
    """
    Calls the live ENA XML API.  Run with:
        pytest -m integration tests/unit/test_generate_inputs.py
    """

    # A known assembly from study ERP190907
    ASSEMBLY_ACC = 'ERZ3834802'

    @pytest.fixture(autouse=True)
    def real_handler(self):
        """Replace the module-level mock handler with a real EnaApiHandler instance."""
        from ena_portal_api.ena_handler import EnaApiHandler
        with patch.object(generate_inputs, 'handler', EnaApiHandler()):
            yield

    def test_get_run_from_real_assembly_xml(self):
        """get_run_from_assembly_xml returns a valid run accession for a real ERZ."""
        result = generate_inputs.get_run_from_assembly_xml(make_assembly(self.ASSEMBLY_ACC))
        assert result is not None
        assert result.startswith(('ERR', 'DRR', 'SRR')), f"Unexpected accession: {result}"

    def test_fetch_data_erp190907_prjeb31790(self):
        """
        Smoke-test fetch_data with the real studies.
        Checks that at least one assembly-run pair is resolved and written to disk.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            assembly_run, assembly_path, run_paths = generate_inputs.fetch_data(
                raw_reads_study='PRJEB31790',
                assembly_study='ERP190907',
                outdir=tmpdir,
                input_scientific_name=None,
                input_env_biome=None,
                keep_metat=False,
            )
        assert len(assembly_run) > 0, "Expected at least one assembly-run pair"
        for erz, run in assembly_run.items():
            assert run.startswith(('ERR', 'DRR', 'SRR')), f"Bad run for {erz}: {run}"
