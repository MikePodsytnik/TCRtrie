from __future__ import annotations

from pathlib import Path
import sqlite3
from typing import Optional

import pandas as pd

from ._tcrtrie import Trie


class VDJdbClient:
    def __init__(self, *, trie: Trie, sqlite_path: Path):
        self._trie = trie
        self._sqlite_path = sqlite_path
        self._df = self._load_table()

    def _load_table(self) -> pd.DataFrame:
        con = sqlite3.connect(self._sqlite_path)
        try:
            df = pd.read_sql_query("SELECT * FROM vdjdb", con)
        finally:
            con.close()
        return df

    @property
    def trie(self) -> Trie:
        return self._trie

    def search(
        self,
        *,
        query: str,
        maxSubstitution: int = 0,
        maxInsertion: int = 0,
        maxDeletion: int = 0,
        maxEdits: Optional[int] = None,
        vGeneFilter: Optional[str] = None,
        jGeneFilter: Optional[str] = None,
    ) -> pd.DataFrame:
        raw = self._trie.SearchIndices(
            query=query,
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
            vGeneFilter=vGeneFilter,
            jGeneFilter=jGeneFilter,
        )

        if not raw:
            return pd.DataFrame()

        idxs = [int(i) for i, _ in raw]
        dists = [int(d) for _, d in raw]

        df = self._df[self._df["idx"].isin(set(idxs))].copy()

        df["_distance"] = df["idx"].map(dict(zip(idxs, dists)))
        order = {idx: pos for pos, idx in enumerate(idxs)}
        df["_order"] = df["idx"].map(order)
        df = df.sort_values("_order").drop(columns=["_order"])
        df = df.rename(columns={"_distance": "distance"})
        df = df.drop('idx', axis=1)
        df = df.drop('complex.id', axis=1)
        col = df.pop('distance')
        df.insert(0, 'distance', col)

        return df
