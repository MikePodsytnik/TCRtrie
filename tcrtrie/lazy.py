from __future__ import annotations

from typing import Any, Callable


class LazyObject:
    def __init__(self, factory: Callable[[], Any]):
        self._factory = factory
        self._obj: Any | None = None

    def _get(self) -> Any:
        if self._obj is None:
            self._obj = self._factory()
        return self._obj

    def __getattr__(self, name: str) -> Any:
        return getattr(self._get(), name)

    def __repr__(self) -> str:
        if self._obj is None:
            return "<LazyObject (not initialized)>"
        return repr(self._obj)
