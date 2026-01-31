from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sqlite3
from typing import Any

import pandas as pd

from ._tcrtrie import Trie


class VDJdbClient:
    def __init__(self, *, trie: Trie, sqlite_path: Path):
        self._trie = trie
        self._sqlite_path = sqlite_path

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
            maxEdits: int | None = None,
            vGeneFilter: str | None = None,
            jGeneFilter: str | None = None,
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

        con = sqlite3.connect(self._sqlite_path)
        try:
            q = f"SELECT * FROM vdjdb WHERE idx IN ({','.join(['?']*len(idxs))})"
            df = pd.read_sql_query(q, con, params=idxs)
        finally:
            con.close()

        df["_distance"] = df["idx"].map(dict(zip(idxs, dists)))
        order = {idx: pos for pos, idx in enumerate(idxs)}
        df["_order"] = df["idx"].map(order)
        df = df.sort_values("_order").drop(columns=["_order"])

        df = df.rename(columns={"_distance": "distance"})

        return df
