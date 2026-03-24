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
    def __init__(self, *, trie: Trie, sqlite_path: Path):
        self._trie = trie
        self._sqlite_path = sqlite_path

        self._df = self._load_table()
        self._columns = list(self._df.columns)
        self._values = self._df.to_numpy(copy=False)

        self._empty = pd.DataFrame(columns=["target_index", "distance", *self._columns])
        self._empty_all = pd.DataFrame(columns=["query_index", "query", "target_index", "distance", *self._columns])
        self._empty_matrix = pd.DataFrame(columns=["target_index", "cost", *self._columns])
        self._empty_all_matrix = pd.DataFrame(columns=["query_index", "query", "target_index", "cost", *self._columns])

        matrix_res = ir.files(__package__).joinpath("data/matrices/blosum62.txt")
        with ir.as_file(matrix_res) as p:
            self._trie.LoadSubstitutionMatrix(str(p), "", 1.5)

    def _load_table(self) -> pd.DataFrame:
        con = sqlite3.connect(self._sqlite_path)
        try:
            return pd.read_sql_query("SELECT * FROM vdjdb", con)
        finally:
            con.close()

    def _normalize_query(self, query: str) -> str:
        if not isinstance(query, str):
            raise TypeError("query must be a string")
        query = query.strip().upper()
        if not query:
            raise ValueError("query must not be empty")
        return query

    def _normalize_queries(self, queries: Sequence[str]) -> list[str]:
        return [self._normalize_query(q) for q in queries]

    def _normalize_optional_filters(self, filters: Optional[Sequence[str]]) -> Optional[list[str]]:
        if filters is None:
            return None
        return [self._normalize_query(f) for f in filters]

    def _resolve_num_threads(self, numThreads: Optional[int]) -> int:
        if numThreads is None:
            return 4
        value = int(numThreads)
        if value < 1:
            raise ValueError("numThreads must be >= 1")
        return value

    def _build_hits_frame(self, idxs: np.ndarray, values: np.ndarray, score_name: str) -> pd.DataFrame:
        selected = self._values[idxs]
        out = pd.DataFrame(selected, columns=self._columns, copy=False)
        out.insert(0, score_name, values)
        out.insert(0, "target_index", idxs)
        return out

    def _with_alignment_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy(deep=False)
        if "alignment_query" not in out.columns:
            out.insert(len(out.columns), "alignment_query", pd.Series(dtype="object"))
        if "alignment_cdr3" not in out.columns:
            out.insert(len(out.columns), "alignment_cdr3", pd.Series(dtype="object"))
        if "alignment_meta" not in out.columns:
            out.insert(len(out.columns), "alignment_meta", pd.Series(dtype="object"))
        return out

    def _extract_int_hits(self, hits: Sequence[tuple[int, int]] | pd.DataFrame) -> list[tuple[int, int]]:
        if isinstance(hits, pd.DataFrame):
            if hits.empty:
                return []
            if "target_index" not in hits.columns:
                raise ValueError("hits DataFrame must contain 'target_index'")
            if "distance" not in hits.columns:
                raise ValueError("hits DataFrame must contain 'distance'")
            idxs = hits["target_index"].to_numpy(dtype=np.int64, copy=False)
            dists = hits["distance"].to_numpy(dtype=np.int64, copy=False)
            return list(zip(idxs.tolist(), dists.tolist()))
        return [(int(idx), int(dist)) for idx, dist in hits]

    def _extract_float_hits(self, hits: Sequence[tuple[int, float]] | pd.DataFrame) -> list[tuple[int, float]]:
        if isinstance(hits, pd.DataFrame):
            if hits.empty:
                return []
            if "target_index" not in hits.columns:
                raise ValueError("hits DataFrame must contain 'target_index'")
            if "cost" not in hits.columns:
                raise ValueError("hits DataFrame must contain 'cost'")
            idxs = hits["target_index"].to_numpy(dtype=np.int64, copy=False)
            costs = hits["cost"].to_numpy(dtype=np.float64, copy=False)
            return list(zip(idxs.tolist(), costs.tolist()))
        return [(int(idx), float(cost)) for idx, cost in hits]

    def _op_type_name(self, op_type: Any) -> str:
        name = getattr(op_type, "name", None)
        if isinstance(name, str) and name:
            return name
        return str(op_type).split(".")[-1]

    def _alignment_to_json(self, alignment: Any) -> str:
        ops = []
        for op in alignment.ops:
            ops.append(
                {
                    "type": self._op_type_name(op.type),
                    "query_pos": int(op.queryPos),
                    "query_char": op.queryChar,
                    "target_char": op.targetChar,
                }
            )
        distance_value = float(alignment.distance)
        if distance_value.is_integer():
            distance_json: int | float = int(distance_value)
        else:
            distance_json = distance_value
        payload = {
            "distance": distance_json,
            "substitutions": int(alignment.substitutions),
            "insertions": int(alignment.insertions),
            "deletions": int(alignment.deletions),
            "ops": ops,
        }
        return json.dumps(payload, ensure_ascii=False)

    def _attach_detailed_alignment(
            self,
            df: pd.DataFrame,
            *,
            max_workers: int,
            use_matrix: bool,
            maxSubstitution: Optional[int] = None,
            maxInsertion: Optional[int] = None,
            maxDeletion: Optional[int] = None,
            maxEdits: Optional[int] = None,
            maxCost: Optional[float] = None,
    ) -> pd.DataFrame:
        out = self._with_alignment_columns(df)
        if out.empty:
            return out

        alignment_query: list[Optional[str]] = [None] * len(out)
        alignment_cdr3: list[Optional[str]] = [None] * len(out)
        alignment_meta: list[Optional[str]] = [None] * len(out)

        group_keys = ["query_index", "query"] if "query_index" in out.columns else ["query"]

        for _, group in out.groupby(group_keys, sort=False):
            query_value = str(group.iloc[0]["query"])
            row_positions = [int(i) for i in group.index]

            if use_matrix:
                hits = list(zip(group["target_index"].astype(int).tolist(), group["cost"].astype(float).tolist()))
                alignments = self._trie.AlignIndexHitsWithMatrix(
                    query=query_value,
                    hits=hits,
                    maxCost=maxCost,
                    numThreads=max_workers,
                )
            else:
                hits = list(zip(group["target_index"].astype(int).tolist(), group["distance"].astype(int).tolist()))
                alignments = self._trie.AlignIndexHits(
                    query=query_value,
                    hits=hits,
                    maxSubstitution=maxSubstitution,
                    maxInsertion=maxInsertion,
                    maxDeletion=maxDeletion,
                    maxEdits=maxEdits,
                    numThreads=max_workers,
                )

            for row_pos, alignment in zip(row_positions, alignments):
                if alignment is None:
                    continue
                alignment_query[row_pos] = alignment.queryAligned
                alignment_cdr3[row_pos] = alignment.targetAligned
                alignment_meta[row_pos] = self._alignment_to_json(alignment)

        out["alignment_query"] = alignment_query
        out["alignment_cdr3"] = alignment_cdr3
        out["alignment_meta"] = alignment_meta
        return out

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
            detailed: bool = False,
            numThreads: int = 4,
    ) -> pd.DataFrame:
        normalized_query = self._normalize_query(query)
        resolved_threads = self._resolve_num_threads(numThreads)

        raw = self._trie.SearchIndices(
            query=normalized_query,
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
            vGeneFilter=None if vGeneFilter is None else self._normalize_query(vGeneFilter),
            jGeneFilter=None if jGeneFilter is None else self._normalize_query(jGeneFilter),
        )

        if not raw:
            base = self._empty.copy()
            base.insert(0, "query", pd.Series(dtype="object"))
            return self._with_alignment_columns(base) if detailed else base

        arr = np.asarray(raw, dtype=np.int64)
        idxs = arr[:, 0]
        dists = arr[:, 1]
        out = self._build_hits_frame(idxs, dists, "distance")
        out.insert(0, "query", normalized_query)

        if not detailed:
            return out

        return self._attach_detailed_alignment(
            out,
            max_workers=resolved_threads,
            use_matrix=False,
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
        )

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
            detailed: bool = False,
            numThreads: int = 4,
    ) -> pd.DataFrame:
        qs = self._normalize_queries(queries)
        resolved_threads = self._resolve_num_threads(numThreads)

        if not qs:
            return self._with_alignment_columns(self._empty_all.copy()) if detailed else self._empty_all.copy()

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
            vGeneFilters=self._normalize_optional_filters(vGeneFilters),
            jGeneFilters=self._normalize_optional_filters(jGeneFilters),
            numThreads=resolved_threads,
        )

        if not raw:
            return self._with_alignment_columns(self._empty_all.copy()) if detailed else self._empty_all.copy()

        frames: list[pd.DataFrame] = []

        for i, hits in enumerate(raw):
            if not hits:
                continue

            arr = np.asarray(hits, dtype=np.int64)
            idxs = arr[:, 0]
            dists = arr[:, 1]

            out = self._build_hits_frame(idxs, dists, "distance")
            out.insert(0, "query", qs[i])
            out.insert(0, "query_index", i)
            frames.append(out)

        if not frames:
            return self._with_alignment_columns(self._empty_all.copy()) if detailed else self._empty_all.copy()

        out = pd.concat(frames, ignore_index=True, copy=False)
        if not detailed:
            return out

        return self._attach_detailed_alignment(
            out,
            max_workers=resolved_threads,
            use_matrix=False,
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
        )

    def search_with_matrix(
            self,
            *,
            query: str,
            maxCost: float,
            vGeneFilter: Optional[str] = None,
            jGeneFilter: Optional[str] = None,
            detailed: bool = False,
            numThreads: int = 4,
    ) -> pd.DataFrame:
        normalized_query = self._normalize_query(query)
        resolved_threads = self._resolve_num_threads(numThreads)

        raw = self._trie.SearchIndicesWithMatrix(
            query=normalized_query,
            maxCost=maxCost,
            vGeneFilter=None if vGeneFilter is None else self._normalize_query(vGeneFilter),
            jGeneFilter=None if jGeneFilter is None else self._normalize_query(jGeneFilter),
        )

        if not raw:
            base = self._empty_matrix.copy()
            base.insert(0, "query", pd.Series(dtype="object"))
            return self._with_alignment_columns(base) if detailed else base

        arr = np.asarray(raw, dtype=np.float64)
        idxs = arr[:, 0].astype(np.int64, copy=False)
        costs = arr[:, 1]
        out = self._build_hits_frame(idxs, costs, "cost")
        out.insert(0, "query", normalized_query)

        if not detailed:
            return out

        return self._attach_detailed_alignment(
            out,
            max_workers=resolved_threads,
            use_matrix=True,
            maxCost=maxCost,
        )

    def search_for_all_with_matrix(
            self,
            *,
            queries: Sequence[str],
            maxCost: float,
            vGeneFilters: Optional[Sequence[str]] = None,
            jGeneFilters: Optional[Sequence[str]] = None,
            detailed: bool = False,
            numThreads: int = 4,
    ) -> pd.DataFrame:
        qs = self._normalize_queries(queries)
        resolved_threads = self._resolve_num_threads(numThreads)

        if not qs:
            return self._with_alignment_columns(self._empty_all_matrix.copy()) if detailed else self._empty_all_matrix.copy()

        if vGeneFilters is not None and len(vGeneFilters) != len(qs):
            raise ValueError("vGeneFilters length must match queries length")
        if jGeneFilters is not None and len(jGeneFilters) != len(qs):
            raise ValueError("jGeneFilters length must match queries length")

        raw = self._trie.SearchIndicesForAllWithMatrix(
            queries=qs,
            maxCost=maxCost,
            vGeneFilters=self._normalize_optional_filters(vGeneFilters),
            jGeneFilters=self._normalize_optional_filters(jGeneFilters),
            numThreads=resolved_threads,
        )

        if not raw:
            return self._with_alignment_columns(self._empty_all_matrix.copy()) if detailed else self._empty_all_matrix.copy()

        frames: list[pd.DataFrame] = []

        for i, hits in enumerate(raw):
            if not hits:
                continue

            arr = np.asarray(hits, dtype=np.float64)
            idxs = arr[:, 0].astype(np.int64, copy=False)
            costs = arr[:, 1]

            out = self._build_hits_frame(idxs, costs, "cost")
            out.insert(0, "query", qs[i])
            out.insert(0, "query_index", i)
            frames.append(out)

        if not frames:
            return self._with_alignment_columns(self._empty_all_matrix.copy()) if detailed else self._empty_all_matrix.copy()

        out = pd.concat(frames, ignore_index=True, copy=False)
        if not detailed:
            return out

        return self._attach_detailed_alignment(
            out,
            max_workers=resolved_threads,
            use_matrix=True,
            maxCost=maxCost,
        )

    def align_query_to_target(
            self,
            *,
            query: str,
            target: str,
            maxSubstitution: Optional[int] = None,
            maxInsertion: Optional[int] = None,
            maxDeletion: Optional[int] = None,
            maxEdits: Optional[int] = None,
    ):
        return self._trie.AlignQueryToTarget(
            query=self._normalize_query(query),
            target=self._normalize_query(target),
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
        )

    def align_hit_by_index(
            self,
            *,
            query: str,
            targetIndex: int,
            maxSubstitution: Optional[int] = None,
            maxInsertion: Optional[int] = None,
            maxDeletion: Optional[int] = None,
            maxEdits: Optional[int] = None,
    ):
        return self._trie.AlignIndexHit(
            query=self._normalize_query(query),
            targetIndex=int(targetIndex),
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
        )

    def align_hits(
            self,
            *,
            query: str,
            hits: Sequence[tuple[int, int]] | pd.DataFrame,
            maxSubstitution: Optional[int] = None,
            maxInsertion: Optional[int] = None,
            maxDeletion: Optional[int] = None,
            maxEdits: Optional[int] = None,
            numThreads: int = 4,
    ):
        normalized_query = self._normalize_query(query)
        extracted_hits = self._extract_int_hits(hits)
        resolved_threads = self._resolve_num_threads(numThreads)
        if not extracted_hits:
            return []
        return self._trie.AlignIndexHits(
            query=normalized_query,
            hits=extracted_hits,
            maxSubstitution=maxSubstitution,
            maxInsertion=maxInsertion,
            maxDeletion=maxDeletion,
            maxEdits=maxEdits,
            numThreads=resolved_threads,
        )

    def align_query_to_target_with_matrix(
            self,
            *,
            query: str,
            target: str,
            maxCost: Optional[float] = None,
    ):
        return self._trie.AlignQueryToTargetWithMatrix(
            query=self._normalize_query(query),
            target=self._normalize_query(target),
            maxCost=maxCost,
        )

    def align_hit_by_index_with_matrix(
            self,
            *,
            query: str,
            targetIndex: int,
            maxCost: Optional[float] = None,
    ):
        return self._trie.AlignIndexHitWithMatrix(
            query=self._normalize_query(query),
            targetIndex=int(targetIndex),
            maxCost=maxCost,
        )

    def align_hits_with_matrix(
            self,
            *,
            query: str,
            hits: Sequence[tuple[int, float]] | pd.DataFrame,
            maxCost: Optional[float] = None,
            numThreads: int = 4,
    ):
        normalized_query = self._normalize_query(query)
        extracted_hits = self._extract_float_hits(hits)
        resolved_threads = self._resolve_num_threads(numThreads)
        if not extracted_hits:
            return []
        return self._trie.AlignIndexHitsWithMatrix(
            query=normalized_query,
            hits=extracted_hits,
            maxCost=maxCost,
            numThreads=resolved_threads,
        )
