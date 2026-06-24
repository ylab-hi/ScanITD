"""Protocol type definitions for ScanITD.

Defines :class:`LoggerType`, a structural subtype compatible with loguru logger
instances, enabling type-safe logger injection throughout the codebase.
"""

from typing import Protocol
from typing import Any


class LoggerType(Protocol):
    """Structural Protocol for loguru-compatible logger objects.

    Any object implementing these methods can be used as a logger in ScanITD.
    This is satisfied at runtime by a configured loguru logger.
    """

    def trace(self, msg: str) -> None:
        """Log a TRACE-level message.

        Args:
            msg: The message to log.
        """

    def debug(self, msg: str) -> None:
        """Log a DEBUG-level message.

        Args:
            msg: The message to log.
        """

    def info(self, msg: str) -> None:
        """Log an INFO-level message.

        Args:
            msg: The message to log.
        """

    def warning(self, msg: str) -> None:
        """Log a WARNING-level message.

        Args:
            msg: The message to log.
        """

    def error(self, msg: str) -> None:
        """Log an ERROR-level message.

        Args:
            msg: The message to log.
        """

    def critical(self, msr: str) -> None:
        """Log a CRITICAL-level message.

        Args:
            msr: The message to log.
        """

    def success(self, msg: str) -> None:
        """Log a SUCCESS-level message.

        Args:
            msg: The message to log.
        """

    def complete(self) -> Any:
        """Complete or flush the logger, returning any completion value."""
