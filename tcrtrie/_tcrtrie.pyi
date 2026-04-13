from __future__ import annotations
from typing import Optional, Sequence, overload


class AIRREntity:
    junctionAA: str
    vGene: str
    jGene: str
    distance: int

class AlignmentOpType: ...
class AlignmentOp:
    type: AlignmentOpType
    queryPos: int
    queryChar: str
    targetChar: str

class AlignmentResult:
    queryAligned: str
    targetAligned: str
    substitutions: int
    insertions: int
    deletions: int
    distance: float
    ops: list[AlignmentOp]

class Trie:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, dataPath: str) -> None: ...
    def LoadSubstitutionMatrix(
        self,
        matrixPath: str,
        delimiter: str = "",
        gapFactor: float = 1.5,
    ) -> None: ...
    def PrintMatrix(self) -> None: ...

    @overload
    def Search(self, query: str, maxEdits: int) -> list[str]: ...

    @overload
    def Search(
        self,
        queries: Sequence[str],
        maxEdits: int,
        numThreads: int | None = 4,
    ) -> list[list[str]]: ...

    def SearchIndices(
        self,
        query: str,
        maxSubstitution: int = 0,
        maxInsertion: int = 0,
        maxDeletion: int = 0,
        maxEdits: int | None = None,
        vGeneFilter: str | None = None,
        jGeneFilter: str | None = None,
    ) -> list[tuple[int, int]]: ...

    def SearchIndicesForAll(
        self,
        queries: Sequence[str],
        maxSubstitution: int = 0,
        maxInsertion: int = 0,
        maxDeletion: int = 0,
        maxEdits: int | None = None,
        vGeneFilters: Sequence[str] | None = None,
        jGeneFilters: Sequence[str] | None = None,
        numThreads: int | None = 4,
    ) -> list[list[tuple[int, int]]]: ...