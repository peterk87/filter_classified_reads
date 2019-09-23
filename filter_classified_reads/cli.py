#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Console script for filter_classified_reads."""
import logging
import sys
from typing import Optional, List

import click

from filter_classified_reads.util import \
    parse_taxids_string, \
    compare_kraken2_and_centrifuge
from filter_classified_reads.io import write_reads_seqtk
from filter_classified_reads.target_classified_reads import \
    TargetClassifiedReads, \
    common_unclassified_reads, \
    find_target_read_ids
from filter_classified_reads.const import CENTRIFUGE, KRAKEN2, LOG_FORMAT


@click.command()
@click.option('-i', '--reads1', type=click.Path(exists=True),
              required=True,
              help='Single end reads or Forward reads if paired-end specified')
@click.option('-I', '--reads2', type=click.Path(exists=True),
              help='Reverse reads [optional]')
@click.option('-c', '--centrifuge-results', type=click.Path(exists=True),
              help='Centrifuge classification results')
@click.option('-C', '--centrifuge-kreport', type=click.Path(exists=True),
              help='Centrifuge Kraken-style report')
@click.option('-k', '--kraken2-results', type=click.Path(exists=True),
              help='Kraken2 classification results')
@click.option('-K', '--kraken2-kreport', type=click.Path(exists=True),
              help='Kraken2 report')
@click.option('-o', '--output1', required=True,
              help='Filtered forward reads or single-end reads.')
@click.option('-O', '--output2',
              help='Filtered reverse reads. Must be specified if providing '
                   'paired end read input!')
@click.option('--exclude-unclassified', is_flag=True,
              help='Do not include unclassified reads in the final output.')
@click.option('--taxids', default=None,
              help=('Optional NCBI Taxonomy ID(s). Comma-delimited with no '
                    'whitespace if more than one to filter for, '
                    'e.g. "1,2,3,4"'))
def main(reads1: str,
         reads2: Optional[str],
         centrifuge_results: Optional[str],
         centrifuge_kreport: Optional[str],
         kraken2_results: Optional[str],
         kraken2_kreport: Optional[str],
         output1: str,
         output2: Optional[str],
         exclude_unclassified: bool,
         taxids: Optional[str]):
    """Filter viral reads and unclassified based on classification results.

    Requires either Kraken2 or Centrifuge classification results or both of a
    FASTQ or 2 FASTQs if paired-end used.

    Note: Kraken-style reports are required along with by read classification
          results for extracting taxonomic hierarchy.
    """

    if not (centrifuge_kreport or centrifuge_results
            or kraken2_kreport or kraken2_results):
        raise click.exceptions.UsageError(
            'No Centrifuge or Kraken2 results and reports specified! Cannot '
            'filter on classification results.')
    if bool(centrifuge_kreport) ^ bool(centrifuge_results):
        raise click.exceptions.UsageError(
            'Both the Centrifuge results and Kraken-style report must be '
            'specified with `-c` for the results file and `-C` for the '
            'Kraken report file!')
    if bool(kraken2_results) ^ bool(kraken2_kreport):
        raise click.exceptions.UsageError(
            'Both the Kraken2 results and report files must be specified '
            'with `-k` for the results file and `-K` for the Kraken report '
            'file!')
    logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)
    parsed_taxids = try_parse_taxids(taxids)
    if reads2:
        if output2 is None:
            raise click.UsageError(f'If paired reads are specified, you must '
                                   f'specify an output file for the filtered '
                                   f'reverse reads with `-O/--output2`!')

    tcr = TargetClassifiedReads()
    if centrifuge_results and centrifuge_kreport:
        tcr = find_target_read_ids(tcr=tcr,
                                   kreport=centrifuge_kreport,
                                   results=centrifuge_results,
                                   method=CENTRIFUGE,
                                   taxids=parsed_taxids)
    if kraken2_results and kraken2_kreport:
        tcr = find_target_read_ids(tcr=tcr,
                                   kreport=kraken2_kreport,
                                   results=kraken2_results,
                                   method=KRAKEN2,
                                   taxids=parsed_taxids)

    target_read_ids = tcr.centrifuge_targets | tcr.kraken2_targets

    unclassified_read_ids = common_unclassified_reads(tcr)
    logging.info(f'Found N={len(unclassified_read_ids)} common unclassified '
                 f'reads by all classification methods.')

    compare_kraken2_and_centrifuge(centrifuge_results,
                                   kraken2_results,
                                   target_read_ids,
                                   tcr)

    if exclude_unclassified:
        filtered_read_ids = list(target_read_ids)
    else:
        filtered_read_ids = list(target_read_ids | unclassified_read_ids)
    # sorting in case that speeds up reading from FASTQ file
    filtered_read_ids.sort()
    if len(filtered_read_ids) == 0:
        logging.warning('No reads found for taxa of interest' +
                        " including unclassified" if not exclude_unclassified
                        else "" + '!')
    else:
        logging.info(f'Writing n={len(filtered_read_ids)} filtered reads '
                     f'from "{reads1}" to "{output1}"')

        write_reads_seqtk(reads1, filtered_read_ids, output1)
        if reads2 and output2 is not None:
            logging.info(f'Writing n={len(filtered_read_ids)} filtered reads '
                         f'from "{reads2}" to "{output2}"')
            write_reads_seqtk(reads2, filtered_read_ids, output2)
    logging.info('Done!')


def try_parse_taxids(taxids: Optional[str]) -> Optional[List[int]]:
    if taxids is None or taxids == '':
        return None
    try:
        return parse_taxids_string(taxids)
    except Exception:
        logging.error(f'Could not parse "taxids" command-line argument '
                      f'({taxids}). Setting to default value of "None"')
        return None


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
