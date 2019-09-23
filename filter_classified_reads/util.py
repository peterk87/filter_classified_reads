from typing import List, Optional, Set, TYPE_CHECKING
import logging
import re
import subprocess as sp

import pandas as pd

if TYPE_CHECKING:
    from filter_classified_reads.target_classified_reads import \
        TargetClassifiedReads  # noqa


def prefix_spaces(s: str) -> int:
    """Count number of prefix spaces in a string."""
    count = 0
    for x in s:
        if x != ' ':
            return count
        count += 1
    return count


def parse_taxids_string(taxids_string: str) -> List[int]:
    return [int(x.strip()) for x in taxids_string.split(',')]


def compare_kraken2_and_centrifuge(centrifuge_results: Optional[str],
                                   kraken2_results: Optional[str],
                                   target_read_ids: Set[str],
                                   tcr: 'TargetClassifiedReads') -> None:
    if centrifuge_results is not None \
        and kraken2_results is not None \
        and tcr.kraken2_unclassified is not None \
        and tcr.kraken2_targets is not None \
        and tcr.centrifuge_unclassified is not None \
        and tcr.centrifuge_targets is not None:
        n_target_uq_c = len(tcr.centrifuge_targets - tcr.kraken2_targets)
        n_target_uq_k2 = len(tcr.kraken2_targets - tcr.centrifuge_targets)
        n_target_total = len(target_read_ids)
        logging.info(f'Total viral reads={n_target_total}')
        logging.info(f'Centrifuge found n={n_target_uq_c} target reads not '
                     f'found with Kraken2')
        logging.info(f'Kraken2 found n={n_target_uq_k2} target reads not found'
                     f' with Centrifuge')

        uc_uq_k2 = tcr.kraken2_unclassified - tcr.centrifuge_unclassified
        uc_uq_c = tcr.centrifuge_unclassified - tcr.kraken2_unclassified
        if tcr.centrifuge_df_results is not None \
            and isinstance(tcr.centrifuge_df_results, pd.DataFrame):
            c_read_ids = set(tcr.centrifuge_df_results.index)
            n_k2_not_in_centrifuge = len(uc_uq_k2 - c_read_ids)
            if n_k2_not_in_centrifuge:
                logging.info(f'N={n_k2_not_in_centrifuge} Unclassified reads '
                             f'by Kraken2 not in Centrifuge results')
        if tcr.kraken2_df_results is not None \
            and isinstance(tcr.kraken2_df_results, pd.DataFrame):
            k2_read_ids = set(tcr.kraken2_df_results.index)
            n_c_not_in_k2 = len(uc_uq_c - k2_read_ids)
            if n_c_not_in_k2:
                logging.info(f'N={n_c_not_in_k2} Unclassified reads by '
                             f'Centrifuge not in Kraken2 results')
        if tcr.centrifuge_unclassified and tcr.kraken2_unclassified:
            uc_intersect = tcr.centrifuge_unclassified \
                           & tcr.kraken2_unclassified
            logging.info(f'N={len(uc_intersect)} reads unclassified by both '
                         f'Centrifuge and Kraken2.')

def check_bin(bin, version_pattern=None) -> Optional[str]:
    """Check if a binary app exists else raise a FileNotFoundError"""
    try:
        p = sp.Popen([bin], stderr=sp.PIPE)
        _, stderr = p.communicate()
        if version_pattern:
            m = re.search(version_pattern, stderr.decode())
            if m:
                return m.group(1)
            else:
                raise FileNotFoundError()
        else:
            return None
    except FileNotFoundError:
        raise FileNotFoundError(
            'You must have "seqtk" installed to run this program!')

def check_seqtk() -> str:
    return check_bin('seqtk', r'Version:\s*(\S+)')

def check_pbgzip() -> None:
    return check_bin('pbgzip')
