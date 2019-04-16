from typing import TextIO, Iterable

import pandas as pd
import screed


def read_kraken_report(path):
    fields = 'perc n_reads n_reads_specific rank taxid sciname'.split()
    return pd.read_csv(path, sep='\t', header=None, names=fields)


def read_kraken2_results(path: str) -> pd.DataFrame:
    kraken2_fields = [('is_classified', 'category'),
                      ('readID', str),
                      ('taxID', 'uint32'),
                      ('queryLength', 'uint16'),
                      ('LCA_mapping', str)]
    return pd.read_csv(path, sep='\t', header=None,
                       names=[k for k, v in kraken2_fields],
                       dtype={k: v for k, v in kraken2_fields}) \
             .set_index('readID')


def read_centrifuge_results(path: str) -> pd.DataFrame:
    centrifuge_results_dtypes = {
        'readID': str,
        'seqID': 'category',
        'taxID': 'uint32',
        'score': 'uint32',
        '2ndBestScore': 'uint32',
        'hitLength': 'uint16',
        'queryLength': 'uint16',
        'numMatches': 'uint8', }
    return pd.read_csv(path,
                       sep='\t',
                       dtype=centrifuge_results_dtypes) \
        .set_index('readID')


def read_fai(path: str) -> pd.DataFrame:
    """Read samtools faidx tab-delimited table into pd.DataFrame."""
    cols = 'name length offset linebases linewidth qualoffset'.split()
    return pd.read_csv(path,
                       names=cols,
                       sep='\t')\
             .set_index('name')


def fq_entry(fh: TextIO, read_id: str, read_index: pd.Series) -> str:
    """Get the FASTQ read entry corresponding to samtools faidx information."""
    fh.seek(read_index['offset'])
    seq = fh.read(read_index['length'])
    fh.seek(read_index['qualoffset'])
    qual = fh.read(read_index['length'])
    return f'@{read_id}\n{seq}\n+\n{qual}\n'


def write_reads(screed_db: screed.openscreed.ScreedDB,
                read_ids: Iterable[str],
                out: str) -> None:
    with open(out, 'w') as fh_out:
        for read_id in read_ids:
            rec: screed.screedRecord.Record = screed_db[read_id]
            fh_out.write(f'@{read_id} {rec["annotations"]}\n'
                         f'{rec["sequence"]}\n'
                         f'+\n'
                         f'{rec["quality"]}\n')
