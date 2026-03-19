from __future__ import annotations

from pathlib import Path
import sqlite3
from typing import Optional, Sequence
import importlib.resources as ir

import numpy as np
import pandas as pd

from ._tcrtrie import Trie


class VDJdbClient:
    def __init__(self, *, trie: Trie, sqlite_path: Path):
        self._trie = trie
        self._sqlite_path = sqlite_path

        self._df = self._load_table()
        self._columns = list(self._df.columns)
        self._values = self._df.to_numpy(copy=False)

        self._empty = pd.DataFrame(columns=["distance", *self._columns])
        self._empty_all = pd.DataFrame(columns=["query", "distance", *self._columns])

        matrix_res = ir.files(__package__).joinpath("data/matrices/blosum62.txt")
        with ir.as_file(matrix_res) as p:
            self._trie.LoadSubstitutionMatrix(str(p), "", 1.5)

    def _load_table(self) -> pd.DataFrame:
        con = sqlite3.connect(self._sqlite_path)
        try:
            return pd.read_sql_query("SELECT * FROM vdjdb", con)
        finally:
            con.close()

    def _normalize_query(self, query) -> str:
        if not isinstance(query, str):
            raise TypeError("query must be a string")
        query = query.strip().upper()
        if not query:
            raise ValueError("query must not be empty")
        return query

    def _normalize_queries(self, queries: Sequence[str]) -> list[str]:
        return [self._normalize_query(q) for q in queries]

    @property
    def trie(self) -> Trie:
        return self._trie

    @property
    def df(self) -> pd.DataFrame:
        return self._df.copy(deep=False)

    def to_pandas(self) -> pd.DataFrame:
        return self._df.copy(deep=False)

    def __dataframe__(self, nan_as_null: bool = False, allow_copy: bool = True):
        return self._df.copy(deep=False).__dataframe__(nan_as_null=nan_as_null, allow_copy=allow_copy)

    def load_matrix(
            self,
            path: str | Path,
            delimiter: str = "",
            gapFactor: float = 1.5,
    ) -> None:
        self._trie.LoadSubstitutionMatrix(str(path), delimiter, gapFactor)

    def print_matrix(self) -> None:
        self._trie.PrintMatrix()

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
        query = self._normalize_query(query)

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
            return self._empty.copy()

        arr = np.asarray(raw, dtype=np.int64)
        idxs = arr[:, 0]
        dists = arr[:, 1]

        selected = self._values[idxs]
        out = pd.DataFrame(selected, columns=self._columns, copy=False)
        out.insert(0, "distance", dists)
        return out

    def search_for_all(
        self,
        *,
        queries: Sequence[str],
        maxSubstitution: int = 0,
        maxInsertion: int = 0,
        maxDeletion: int = 0,
        maxEdits: Optional[int] = None,
        vGeneFilters: Optional[Sequence[str]] = None,
        jGeneFilters: Optional[Sequence[str]] = None,
    ) -> pd.DataFrame:
        qs = self._normalize_queries(queries)

        if not qs:
            return self._empty_all.copy()

        if vGeneFilters is not None and len(vGeneFilters) != len(qs):
            raise ValueError("vGeneFilters length must match queries length")
        if jGeneFilters is not None and len(jGeneFilters) != len(qs):
            raise ValueError("jGeneFilters length must match queries length")

        raw = self._trie.SearchIndicesForAll(
            queries=qs,
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
            vGeneFilters=None if vGeneFilters is None else list(vGeneFilters),
            jGeneFilters=None if jGeneFilters is None else list(jGeneFilters),
        )

        if not raw:
            return self._empty_all.copy()

        frames: list[pd.DataFrame] = []

        for i, q in enumerate(qs):
            hits = raw.get(q) or []
            if not hits:
                continue

            arr = np.asarray(hits, dtype=np.int64)
            idxs = arr[:, 0]
            dists = arr[:, 1]

            selected = self._values[idxs]
            out = pd.DataFrame(selected, columns=self._columns, copy=False)

            out.insert(0, "distance", dists)
            out.insert(0, "query", q)

            frames.append(out)

        if not frames:
            return self._empty_all.copy()

        return pd.concat(frames, ignore_index=True, copy=False)

    def search_with_matrix(
        self,
        *,
        query: str,
        maxCost: float,
        vGeneFilter: Optional[str] = None,
        jGeneFilter: Optional[str] = None,
    ) -> pd.DataFrame:
        query = self._normalize_query(query)

        raw = self._trie.SearchIndicesWithMatrix(
            query=query,
            maxCost=maxCost,
            vGeneFilter=vGeneFilter,
            jGeneFilter=jGeneFilter,
        )

        if not raw:
            return pd.DataFrame(columns=["cost", *self._columns])

        arr = np.asarray(raw, dtype=np.float64)
        idxs = arr[:, 0].astype(np.int64, copy=False)
        costs = arr[:, 1]

        selected = self._values[idxs]
        out = pd.DataFrame(selected, columns=self._columns, copy=False)
        out.insert(0, "cost", costs)
        return out

    def search_for_all_with_matrix(
        self,
        *,
        queries: Sequence[str],
        maxCost: float,
        vGeneFilters: Optional[Sequence[str]] = None,
        jGeneFilters: Optional[Sequence[str]] = None,
    ) -> pd.DataFrame:
        qs = self._normalize_queries(queries)

        if not qs:
            return pd.DataFrame(columns=["query", "cost", *self._columns])

        if vGeneFilters is not None and len(vGeneFilters) != len(qs):
            raise ValueError("vGeneFilters length must match queries length")
        if jGeneFilters is not None and len(jGeneFilters) != len(qs):
            raise ValueError("jGeneFilters length must match queries length")

        raw = self._trie.SearchIndicesForAllWithMatrix(
            queries=qs,
            maxCost=maxCost,
            vGeneFilters=None if vGeneFilters is None else list(vGeneFilters),
            jGeneFilters=None if jGeneFilters is None else list(jGeneFilters),
        )

        if not raw:
            return pd.DataFrame(columns=["query", "cost", *self._columns])

        frames: list[pd.DataFrame] = []

        for i, q in enumerate(qs):
            hits = raw.get(q) or []
            if not hits:
                continue

            arr = np.asarray(hits, dtype=np.float64)
            idxs = arr[:, 0].astype(np.int64, copy=False)
            costs = arr[:, 1]

            selected = self._values[idxs]
            out = pd.DataFrame(selected, columns=self._columns, copy=False)

            out.insert(0, "cost", costs)
            out.insert(0, "query", q)

            frames.append(out)

        if not frames:
            return pd.DataFrame(columns=["query", "cost", *self._columns])

        return pd.concat(frames, ignore_index=True, copy=False)