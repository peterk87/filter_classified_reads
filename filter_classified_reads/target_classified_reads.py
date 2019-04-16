import logging
from typing import Set, Optional, List

import pandas as pd
import attr

from filter_classified_reads.const import \
    classification_methods, \
    CENTRIFUGE, \
    VIRUSES_TAXID
from filter_classified_reads.io import \
    read_centrifuge_results, \
    read_kraken2_results, \
    read_kraken_report
from filter_classified_reads.tax_node import TaxNode


@attr.s
class TargetClassifiedReads:
    centrifuge_targets: Set[str] = attr.ib(factory=set)
    centrifuge_unclassified: Optional[Set[str]] = attr.ib(default=None)
    centrifuge_df_results: Optional[pd.DataFrame] = attr.ib(default=None)
    kraken2_targets: Set[str] = attr.ib(factory=set)
    kraken2_unclassified: Optional[Set[str]] = attr.ib(default=None)
    kraken2_df_results: Optional[pd.DataFrame] = attr.ib(default=None)


def common_unclassified_reads(tcr: TargetClassifiedReads) -> Set[str]:
    """Get common unclassified read IDs for all classification methods"""
    all_unclassified_by_method = [getattr(tcr, x, None) for x in
                                  tcr.__dict__.keys()
                                  if x.endswith('_unclassified')]
    filt = [x for x in all_unclassified_by_method if x is not None]
    unclassified_read_ids, *rest_uc = filt
    for uc in rest_uc:
        unclassified_read_ids &= uc
    return unclassified_read_ids


def find_target_read_ids(tcr: TargetClassifiedReads,
                         kreport: str,
                         results: str,
                         taxids: List[int] = None,
                         method: str = 'centrifuge') -> TargetClassifiedReads:
    assert method in classification_methods, (f'Cannot handle classification '
                                              f'results of method="{method}"! '
                                              f'Can only handle one of these: '
                                              f'{classification_methods}')
    logging.info(f'Parsing {method} results into DataFrame')
    if method == CENTRIFUGE:
        df_results = read_centrifuge_results(results)
    else:
        df_results = read_kraken2_results(results)

    tcr.__dict__[f'{method}_df_results'] = df_results

    logging.info(f'Parsed n={df_results.shape[0]} {method} '
                 f'result records into DataFrame from "{results}"')
    df_kreport = read_kraken_report(kreport)
    logging.info(f'Parsed n={df_kreport.shape[0]} {method} '
                 f'Kraken-style report records into DataFrame from '
                 f'"{kreport}"')
    df_unclassified = subset_unclassified(df_results)
    unclassified_read_ids = set(df_unclassified.index)
    logging.info(f'Found {len(unclassified_read_ids)} unclassified reads from '
                 f'Centrifuge results')
    tcr.__dict__[f'{method}_unclassified'] = unclassified_read_ids

    if taxids:
        if (df_kreport.taxid.isin(taxids)).sum() == 0:
            logging.warning(f'No taxonomic classification matches to '
                            f'taxids={taxids} in {method} results: {results}')
            tcr.__dict__[f'{method}_targets'] = set()
            return tcr
    else:
        if (df_kreport.taxid.isin([VIRUSES_TAXID])).sum() == 0:
            logging.warning(
                f'No taxonomic classification matches to Viruses '
                f'taxid={VIRUSES_TAXID} in {method} results: {results}')
            tcr.__dict__[f'{method}_targets'] = set()
            return tcr

    tax_tree = TaxNode.build_taxonomy_tree(df_kreport)
    if taxids:
        all_taxids: Set[int] = set()
        for taxid in taxids:
            node = tax_tree.search(taxid)
            if node is not None:
                all_taxids |= node.taxids_set()
        logging.info(f'From input taxids ({taxids}), found {len(all_taxids)} '
                     f'unique taxids including descendants.')
    else:
        node = tax_tree.viral_tax_node()
        if node is not None:
            all_taxids = node.taxids_set()
        else:
            all_taxids = set()
        logging.info(f'Found {len(all_taxids)} unique viral Taxonomy IDs')
    df_target_taxids = subset_classifications_by_taxids(df_results, all_taxids)
    target_read_ids = set(df_target_taxids.index)
    logging.info(f'Found {len(target_read_ids)} target reads from {method} '
                 f'results')
    tcr.__dict__[f'{method}_targets'] = target_read_ids
    return tcr


def subset_classifications_by_taxids(df: pd.DataFrame,
                                     taxids: Set[int]) -> pd.DataFrame:
    return df[df.taxID.isin(taxids)]


def subset_unclassified(df: pd.DataFrame) -> pd.DataFrame:
    return subset_classifications_by_taxids(df, {0})
