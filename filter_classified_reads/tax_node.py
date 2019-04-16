from collections.__init__ import defaultdict
from typing import Optional, List, Mapping, Iterator, Set

import attr
import pandas as pd

from filter_classified_reads.const import VIRUSES_TAXID
from filter_classified_reads.util import prefix_spaces


@attr.s
class TaxNode:
    parent: Optional['TaxNode'] = attr.ib(default=None)
    name: Optional[str] = attr.ib(default=None)
    taxid: Optional[int] = attr.ib(default=None)
    spaces: int = attr.ib(default=0)
    children: List['TaxNode'] = attr.ib(factory=list)
    rank: Optional[str] = attr.ib(default=None)

    @classmethod
    def build_taxonomy_tree(cls, df_kreport: pd.DataFrame) -> 'TaxNode':
        """Construct a taxonomy tree from Kraken-style report."""
        root_node = cls(name='root', taxid=1)
        # keep track of previously added nodes by number of spaces
        nodes: Mapping[int, List['TaxNode']] = defaultdict(list)
        nodes[0] += [root_node]
        for idx, row in df_kreport.iterrows():
            sciname = row.sciname
            if sciname == 'unclassified' or sciname == 'root':
                continue
            spaces = prefix_spaces(sciname)
            # print(sciname, spaces)
            parent_node: 'TaxNode' = nodes[spaces - 2][-1]
            node = cls(name=sciname.strip(),
                       taxid=row['taxid'],
                       parent=parent_node,
                       spaces=spaces,
                       rank=row['rank'])
            nodes[spaces] += [node]
            parent_node.children.append(node)
        return root_node

    def search(self, taxid: int) -> Optional['TaxNode']:
        if self.taxid == taxid:
            return self
        if len(self.children) > 0:
            for child in self.children:
                found_node = child.search(taxid)
                if found_node:
                    return found_node
        return None

    def iter_children(self) -> Iterator['TaxNode']:
        for child in self.children:
            yield child
            if len(child.children) > 0:
                yield from child.iter_children()

    def viral_tax_node(self: 'TaxNode') -> Optional['TaxNode']:
        # superkingdom, viruses: https://www.ncbi.nlm.nih.gov/taxonomy/10239
        return self.search(VIRUSES_TAXID)

    def taxids_set(self) -> Set[int]:
        """Collect specified TaxNode and descendent taxids into a set."""
        taxids = {x.taxid for x in self.iter_children() if x.taxid is not None}
        if self.taxid is not None:
            taxids.add(self.taxid)
        return taxids
