import os
from ._tcrtrie import Trie, AIRREntity

__all__ = ["Trie", "AIRREntity", "VDJdb"]

_here = os.path.dirname(__file__)
VDJdb = Trie(os.path.join(_here, "vdjdb_airr.tsv"))
