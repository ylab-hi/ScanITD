"""Writer for VCF file for itd object."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import IO, Any

from loguru import logger


class Writer(ABC):
    """Abstract class for writing object to file."""

    def __init__(self, file_path: str) -> None:
        """Initialize Writer object."""
        self.file_path = Path(file_path)
        if self.file_path.exists():
            logger.warning(f"{self.file_path} exists, will be overwritten.")
        self.io: IO | None = None

    @abstractmethod
    def write_data(self, data_object: Any, object_id: str):
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @abstractmethod
    def write_line(self, line: str):
        """Write line to file.

        :param: line: Line to write to file.
        """

    @abstractmethod
    def open(self, mode: str = "w"):
        """Open file.

        :param: mode: Mode to open file.
        """

    @abstractmethod
    def close(self):
        """Close file."""
