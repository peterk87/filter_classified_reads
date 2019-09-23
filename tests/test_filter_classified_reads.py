#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `filter_classified_reads` package."""
import os

from click.testing import CliRunner

from filter_classified_reads.const import VIRUSES_TAXID
from filter_classified_reads import cli
from filter_classified_reads.target_classified_reads import \
    common_unclassified_reads, \
    TargetClassifiedReads
from filter_classified_reads.io import read_kraken_report
from filter_classified_reads.tax_node import TaxNode

r1 = os.path.abspath(
    'tests/data/SRR8207674_1.viral_unclassified.seqtk_seed42_n10000.fastq.gz')
r2 = os.path.abspath(
    'tests/data/SRR8207674_2.viral_unclassified.seqtk_seed42_n10000.fastq.gz')
k2_report = os.path.abspath('tests/data/SRR8207674-kraken2_report.tsv')
k2_results = os.path.abspath('tests/data/SRR8207674-kraken2_results.tsv')
c_report = os.path.abspath('tests/data/SRR8207674-centrifuge_kreport.tsv')
c_results = os.path.abspath('tests/data/SRR8207674-centrifuge_results.tsv')


def count_lines(path: str) -> int:
    with os.popen(f'zcat < {path}') as f:
        return sum([1 for l in f])


def test_try_parse_taxids():
    assert cli.try_parse_taxids('') is None, 'Empty strings must return None'
    assert cli.try_parse_taxids('invalid') is None, \
        'Invalid strings must return None'
    assert cli.try_parse_taxids('1') == [1], \
        'Must parse a singular integer string into a list of int'
    assert cli.try_parse_taxids('1,2,3') == [1, 2, 3], \
        'Must parse multiple integer values from a string delimited by commas'
    assert cli.try_parse_taxids('\t1 \n,2 , 3\t ') == [1, 2, 3], \
        'Whitespace must be stripped away!'


def test_common_unclassified_reads():
    tcr = TargetClassifiedReads(centrifuge_unclassified={x for x in 'abc'},
                                kraken2_unclassified={x for x in 'bcd'})
    assert common_unclassified_reads(tcr) == {'b', 'c'}
    tcr.__dict__['test_unclassified'] = {'c'}
    assert common_unclassified_reads(tcr) == {'c'}, \
        'Must return common unclassified read IDs for all "*_unclassified" ' \
        'attributes'


def test_build_taxonomy_tree():
    kreport_fields = 'perc n_reads n_reads_specific rank taxid sciname'.split()
    df_c_kreport = read_kraken_report(c_report)
    assert list(df_c_kreport.columns) == kreport_fields, \
        f'kreport DataFrame must contain the following columns: ' \
        f'{kreport_fields}'
    assert df_c_kreport.shape[0] > 1, \
        'kreport DataFrame must have more than 1 record'
    tax_tree = TaxNode.build_taxonomy_tree(df_c_kreport)
    assert tax_tree.name == 'root', \
        'The name of the TaxNode returned by TaxNode.build_taxonomy_tree() ' \
        'must be "root"'
    assert tax_tree.taxid == 1, 'The taxid of the root TaxNode must be 1'
    assert tax_tree.spaces == 0, \
        'The root node should have no spaces in its scientific name'
    virus_node = tax_tree.search(VIRUSES_TAXID)
    assert virus_node is not None, \
        'There must be a "Viruses" TaxNode in the taxonomy tree'
    assert virus_node.taxid == VIRUSES_TAXID, \
        f'The "Viruses" TaxNode should have a taxid of {VIRUSES_TAXID}'
    assert virus_node.name == 'Viruses', \
        'The name of the "Viruses" TaxNode must be "Viruses"'
    assert virus_node.spaces > 0, \
        'The "Viruses" TaxNode must have more than 0 spaces in its ' \
        'scientific name'
    assert len(virus_node.taxids_set()) > 2, \
        'The "Viruses" TaxNode must have more than 2 descendants'


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 2
    assert 'Missing option' in result.output
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert 'Show this message and exit.' in help_result.output
    with runner.isolated_filesystem():
        out1 = 'R1.fq.gz'
        out2 = 'R2.fq.gz'
        test_run = runner.invoke(cli.main, ['-i', r1, '-I', r2,
                                            '-o', out1, '-O', out2,
                                            '-c', c_results, '-C', c_report,
                                            '-k', k2_results, '-K', k2_report])
        assert test_run.exit_code == 0
        assert os.path.exists(out1)
        assert count_lines(out1) > 4
        assert os.path.exists(out2)
        assert count_lines(out2) > 4
