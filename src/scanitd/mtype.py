"""Type of the scanitd2."""

from typing import Protocol
from typing import Any


class LoggerType(Protocol):
    """Logger type."""

    def trace(self, msg: str) -> None:
        """Trace."""

    def debug(self, msg: str) -> None:
        """Debug."""

    def info(self, msg: str) -> None:
        """Info."""

    def warning(self, msg: str) -> None:
        """Warning."""

    def error(self, msg: str) -> None:
        """Error."""

    def critical(self, msr: str) -> None:
        """Critical."""

    def success(self, msg: str) -> None:
        """Success."""

    def complete(self) -> Any:
        """Complete."""
