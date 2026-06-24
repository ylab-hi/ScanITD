"""Abstract base class for file writers used in ScanITD output generation."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import IO, Any

from loguru import logger


class Writer(ABC):
    """Abstract class for writing object to file."""

    def __init__(self, file_path: str) -> None:
        """Initialize Writer object.

        Args:
            file_path: Path to the output file.  Will be overwritten if it already exists.
        """
        self.file_path = Path(file_path)
        if self.file_path.exists():
            logger.warning(f"{self.file_path} exists, will be overwritten.")
        self.io: IO | None = None

    @abstractmethod
    def write_data(self, data_object: Any, object_id: str):
        """Write data to file.

        Args:
            data_object: The data object to serialize.
            object_id: String identifier for the object (used in IDs/filenames).
        """

    @abstractmethod
    def write_line(self, line: str):
        """Write line to file.

        Args:
            line: Text line to write (should include trailing newline if needed).
        """

    @abstractmethod
    def open(self, mode: str = "w"):
        """Open file.

        Args:
            mode: File open mode (default ``w``).

        Returns:
            An open file-like IO object.
        """

    @abstractmethod
    def close(self):
        """Close file."""
