from __future__ import annotations

from pathlib import Path
import importlib.resources as ir
import json
import sqlite3
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd

from ._tcrtrie import Trie


class VDJdbClient:
    def __init__(self, *, trie: Trie, sqlitePath: Path):
        self._trie = trie
        self._sqlitePath = sqlitePath

        self._df = self._loadTable()
        self._columns = list(self._df.columns)
        self._values = self._df.to_numpy(copy=False)

        self._empty = pd.DataFrame(columns=["target_index", "distance", *self._columns])
        self._emptyAll = pd.DataFrame(columns=["query_index", "query", "target_index", "distance", *self._columns])
        self._emptyMatrix = pd.DataFrame(columns=["target_index", "cost", *self._columns])
        self._emptyAllMatrix = pd.DataFrame(columns=["query_index", "query", "target_index", "cost", *self._columns])

        matrixRes = ir.files(__package__).joinpath("data/matrices/blosum62.txt")
        with ir.as_file(matrixRes) as path:
            self._trie.LoadSubstitutionMatrix(str(path), "", 1.5)

    def _loadTable(self) -> pd.DataFrame:
        con = sqlite3.connect(self._sqlitePath)
        try:
            return pd.read_sql_query("SELECT * FROM vdjdb", con)
        finally:
            con.close()

    def _normalizeQuery(self, query: str) -> str:
        if not isinstance(query, str):
            raise TypeError("query must be a string")
        normalized = query.strip().upper()
        if not normalized:
            raise ValueError("query must not be empty")
        return normalized

    def _normalizeQueries(self, queries: Sequence[str]) -> list[str]:
        return [self._normalizeQuery(query) for query in queries]

    def _normalizeOptionalFilters(self, filters: Optional[Sequence[str]]) -> Optional[list[str]]:
        if filters is None:
            return None
        return [self._normalizeQuery(value) for value in filters]

    def _resolveNumThreads(self, numThreads: Optional[int]) -> int:
        if numThreads is None:
            return 4
        value = int(numThreads)
        if value < 1:
            raise ValueError("numThreads must be >= 1")
        return value

    def _buildHitsFrame(self, indices: np.ndarray, values: np.ndarray, scoreName: str) -> pd.DataFrame:
        selected = self._values[indices]
        out = pd.DataFrame(selected, columns=self._columns, copy=False)
        out.insert(0, scoreName, values)
        out.insert(0, "target_index", indices)
        return out

    def _withAlignmentColumns(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy(deep=False)
        if "alignment_query" not in out.columns:
            out.insert(len(out.columns), "alignment_query", pd.Series(dtype="object"))
        if "alignment_cdr3" not in out.columns:
            out.insert(len(out.columns), "alignment_cdr3", pd.Series(dtype="object"))
        if "alignment_meta" not in out.columns:
            out.insert(len(out.columns), "alignment_meta", pd.Series(dtype="object"))
        return out

    def _opTypeName(self, opType: Any) -> str:
        name = getattr(opType, "name", None)
        if isinstance(name, str) and name:
            return name
        return str(opType).split(".")[-1]

    def _alignmentToJson(self, alignment: Any) -> str:
        ops = []
        for op in alignment.ops:
            ops.append(
                {
                    "type": self._opTypeName(op.type),
                    "query_pos": int(op.queryPos),
                    "query_char": op.queryChar,
                    "target_char": op.targetChar,
                }
            )

        distanceValue = float(alignment.distance)
        if distanceValue.is_integer():
            distanceJson: int | float = int(distanceValue)
        else:
            distanceJson = distanceValue

        payload = {
            "distance": distanceJson,
            "substitutions": int(alignment.substitutions),
            "insertions": int(alignment.insertions),
            "deletions": int(alignment.deletions),
            "ops": ops,
        }
        return json.dumps(payload, ensure_ascii=False)

    def _alignHits(
        self,
        *,
        query: str,
        hits: Sequence[tuple[int, int]],
        maxSubstitution: Optional[int] = None,
        maxInsertion: Optional[int] = None,
        maxDeletion: Optional[int] = None,
        maxEdits: Optional[int] = None,
        numThreads: int = 4,
    ):
        if not hits:
            return []

        return self._trie.AlignIndexHits(
            query=query,
            hits=[(int(index), int(distance)) for index, distance in hits],
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
            numThreads=self._resolveNumThreads(numThreads),
        )

    def _alignHitsWithMatrix(
        self,
        *,
        query: str,
        hits: Sequence[tuple[int, float]],
        maxCost: Optional[float] = None,
        numThreads: int = 4,
    ):
        if not hits:
            return []

        return self._trie.AlignIndexHitsWithMatrix(
            query=query,
            hits=[(int(index), float(cost)) for index, cost in hits],
            maxCost=maxCost,
            numThreads=self._resolveNumThreads(numThreads),
        )

    def _attachDetailedAlignment(
        self,
        df: pd.DataFrame,
        *,
        maxWorkers: int,
        useMatrix: bool,
        maxSubstitution: Optional[int] = None,
        maxInsertion: Optional[int] = None,
        maxDeletion: Optional[int] = None,
        maxEdits: Optional[int] = None,
        maxCost: Optional[float] = None,
    ) -> pd.DataFrame:
        out = self._withAlignmentColumns(df)
        if out.empty:
            return out

        alignmentQuery: list[Optional[str]] = [None] * len(out)
        alignmentCdr3: list[Optional[str]] = [None] * len(out)
        alignmentMeta: list[Optional[str]] = [None] * len(out)

        groupKeys = ["query_index", "query"] if "query_index" in out.columns else ["query"]

        for _, group in out.groupby(groupKeys, sort=False):
            queryValue = str(group.iloc[0]["query"])
            rowPositions = [int(position) for position in group.index]

            if useMatrix:
                hits = list(zip(group["target_index"].astype(int).tolist(), group["cost"].astype(float).tolist()))
                alignments = self._alignHitsWithMatrix(
                    query=queryValue,
                    hits=hits,
                    maxCost=maxCost,
                    numThreads=maxWorkers,
                )
            else:
                hits = list(zip(group["target_index"].astype(int).tolist(), group["distance"].astype(int).tolist()))
                alignments = self._alignHits(
                    query=queryValue,
                    hits=hits,
                    maxSubstitution=maxSubstitution,
                    maxInsertion=maxInsertion,
                    maxDeletion=maxDeletion,
                    maxEdits=maxEdits,
                    numThreads=maxWorkers,
                )

            for rowPosition, alignment in zip(rowPositions, alignments):
                if alignment is None:
                    continue
                alignmentQuery[rowPosition] = alignment.queryAligned
                alignmentCdr3[rowPosition] = alignment.targetAligned
                alignmentMeta[rowPosition] = self._alignmentToJson(alignment)

        out["alignment_query"] = alignmentQuery
        out["alignment_cdr3"] = alignmentCdr3
        out["alignment_meta"] = alignmentMeta
        return out

    @property
    def df(self) -> pd.DataFrame:
        return self._df.copy(deep=False)

    def toDataFrame(self) -> pd.DataFrame:
        return self._df.copy(deep=False)

    def __dataframe__(self, nan_as_null: bool = False, allow_copy: bool = True):
        return self._df.copy(deep=False).__dataframe__(nan_as_null=nan_as_null, allow_copy=allow_copy)

    def loadMatrix(
        self,
        path: str | Path,
        delimiter: str = "",
        gapFactor: float = 1.5,
    ) -> None:
        self._trie.LoadSubstitutionMatrix(str(path), delimiter, gapFactor)

    def printMatrix(self) -> None:
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
        detailed: bool = False,
        numThreads: int = 4,
    ) -> pd.DataFrame:
        normalizedQuery = self._normalizeQuery(query)
        resolvedThreads = self._resolveNumThreads(numThreads)

        raw = self._trie.SearchIndices(
            query=normalizedQuery,
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
            vGeneFilter=None if vGeneFilter is None else self._normalizeQuery(vGeneFilter),
            jGeneFilter=None if jGeneFilter is None else self._normalizeQuery(jGeneFilter),
        )

        if not raw:
            base = self._empty.copy()
            base.insert(0, "query", pd.Series(dtype="object"))
            return self._withAlignmentColumns(base) if detailed else base

        arr = np.asarray(raw, dtype=np.int64)
        indices = arr[:, 0]
        distances = arr[:, 1]
        out = self._buildHitsFrame(indices, distances, "distance")
        out.insert(0, "query", normalizedQuery)

        if not detailed:
            return out

        return self._attachDetailedAlignment(
            out,
            maxWorkers=resolvedThreads,
            useMatrix=False,
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
        )

    def searchForAll(
        self,
        *,
        queries: Sequence[str],
        maxSubstitution: int = 0,
        maxInsertion: int = 0,
        maxDeletion: int = 0,
        maxEdits: Optional[int] = None,
        vGeneFilters: Optional[Sequence[str]] = None,
        jGeneFilters: Optional[Sequence[str]] = None,
        detailed: bool = False,
        numThreads: int = 4,
    ) -> pd.DataFrame:
        normalizedQueries = self._normalizeQueries(queries)
        resolvedThreads = self._resolveNumThreads(numThreads)

        if not normalizedQueries:
            return self._withAlignmentColumns(self._emptyAll.copy()) if detailed else self._emptyAll.copy()

        if vGeneFilters is not None and len(vGeneFilters) != len(normalizedQueries):
            raise ValueError("vGeneFilters length must match queries length")
        if jGeneFilters is not None and len(jGeneFilters) != len(normalizedQueries):
            raise ValueError("jGeneFilters length must match queries length")

        raw = self._trie.SearchIndicesForAll(
            queries=normalizedQueries,
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
            vGeneFilters=self._normalizeOptionalFilters(vGeneFilters),
            jGeneFilters=self._normalizeOptionalFilters(jGeneFilters),
            numThreads=resolvedThreads,
        )

        if not raw:
            return self._withAlignmentColumns(self._emptyAll.copy()) if detailed else self._emptyAll.copy()

        frames: list[pd.DataFrame] = []

        for queryIndex, hits in enumerate(raw):
            if not hits:
                continue

            arr = np.asarray(hits, dtype=np.int64)
            indices = arr[:, 0]
            distances = arr[:, 1]

            out = self._buildHitsFrame(indices, distances, "distance")
            out.insert(0, "query", normalizedQueries[queryIndex])
            out.insert(0, "query_index", queryIndex)
            frames.append(out)

        if not frames:
            return self._withAlignmentColumns(self._emptyAll.copy()) if detailed else self._emptyAll.copy()

        out = pd.concat(frames, ignore_index=True, copy=False)

        if not detailed:
            return out

        return self._attachDetailedAlignment(
            out,
            maxWorkers=resolvedThreads,
            useMatrix=False,
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
        )

    def searchWithMatrix(
        self,
        *,
        query: str,
        maxCost: float,
        vGeneFilter: Optional[str] = None,
        jGeneFilter: Optional[str] = None,
        detailed: bool = False,
        numThreads: int = 4,
    ) -> pd.DataFrame:
        normalizedQuery = self._normalizeQuery(query)
        resolvedThreads = self._resolveNumThreads(numThreads)

        raw = self._trie.SearchIndicesWithMatrix(
            query=normalizedQuery,
            maxCost=maxCost,
            vGeneFilter=None if vGeneFilter is None else self._normalizeQuery(vGeneFilter),
            jGeneFilter=None if jGeneFilter is None else self._normalizeQuery(jGeneFilter),
        )

        if not raw:
            base = self._emptyMatrix.copy()
            base.insert(0, "query", pd.Series(dtype="object"))
            return self._withAlignmentColumns(base) if detailed else base

        arr = np.asarray(raw, dtype=np.float64)
        indices = arr[:, 0].astype(np.int64, copy=False)
        costs = arr[:, 1]
        out = self._buildHitsFrame(indices, costs, "cost")
        out.insert(0, "query", normalizedQuery)

        if not detailed:
            return out

        return self._attachDetailedAlignment(
            out,
            maxWorkers=resolvedThreads,
            useMatrix=True,
            maxCost=maxCost,
        )

    def searchForAllWithMatrix(
        self,
        *,
        queries: Sequence[str],
        maxCost: float,
        vGeneFilters: Optional[Sequence[str]] = None,
        jGeneFilters: Optional[Sequence[str]] = None,
        detailed: bool = False,
        numThreads: int = 4,
    ) -> pd.DataFrame:
        normalizedQueries = self._normalizeQueries(queries)
        resolvedThreads = self._resolveNumThreads(numThreads)

        if not normalizedQueries:
            return self._withAlignmentColumns(self._emptyAllMatrix.copy()) if detailed else self._emptyAllMatrix.copy()

        if vGeneFilters is not None and len(vGeneFilters) != len(normalizedQueries):
            raise ValueError("vGeneFilters length must match queries length")
        if jGeneFilters is not None and len(jGeneFilters) != len(normalizedQueries):
            raise ValueError("jGeneFilters length must match queries length")

        raw = self._trie.SearchIndicesForAllWithMatrix(
            queries=normalizedQueries,
            maxCost=maxCost,
            vGeneFilters=self._normalizeOptionalFilters(vGeneFilters),
            jGeneFilters=self._normalizeOptionalFilters(jGeneFilters),
            numThreads=resolvedThreads,
        )

        if not raw:
            return self._withAlignmentColumns(self._emptyAllMatrix.copy()) if detailed else self._emptyAllMatrix.copy()

        frames: list[pd.DataFrame] = []

        for queryIndex, hits in enumerate(raw):
            if not hits:
                continue

            arr = np.asarray(hits, dtype=np.float64)
            indices = arr[:, 0].astype(np.int64, copy=False)
            costs = arr[:, 1]

            out = self._buildHitsFrame(indices, costs, "cost")
            out.insert(0, "query", normalizedQueries[queryIndex])
            out.insert(0, "query_index", queryIndex)
            frames.append(out)

        if not frames:
            return self._withAlignmentColumns(self._emptyAllMatrix.copy()) if detailed else self._emptyAllMatrix.copy()

        out = pd.concat(frames, ignore_index=True, copy=False)

        if not detailed:
            return out

        return self._attachDetailedAlignment(
            out,
            maxWorkers=resolvedThreads,
            useMatrix=True,
            maxCost=maxCost,
        )