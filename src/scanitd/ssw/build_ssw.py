"""Hatch custom build hook to compile the vendored SSW C shared library."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class CustomBuildHook(BuildHookInterface):
    """Hatch build hook that compiles the native ``ssw.so`` shared library."""

    def initialize(self, version: str, build_data: dict[str, Any]) -> None:
        """Compile ``ssw.so`` during wheel builds and mark the wheel platform-specific."""
        del version
        if self.target_name != "wheel":
            return

        ssw_dir = Path(__file__).resolve().parent
        make_bin = shutil.which("make")
        if make_bin is None:
            msg = "`make` executable not found; cannot build vendored ssw.so."
            raise RuntimeError(msg)

        subprocess.run([make_bin, "clean"], cwd=ssw_dir, capture_output=True, text=True, check=False)  # noqa: S603
        result = subprocess.run([make_bin], cwd=ssw_dir, capture_output=True, text=True, check=False)  # noqa: S603
        if result.returncode != 0:
            sys.stderr.write("Failed to build ssw library:\n")
            sys.stderr.write(result.stderr)
            msg = "ssw build failed"
            raise RuntimeError(msg)

        so_path = ssw_dir / "ssw.so"
        if not so_path.is_file():
            msg = f"Expected compiled library not found at {so_path}"
            raise RuntimeError(msg)

        build_data["pure_python"] = False
        build_data["infer_tag"] = True
        build_data.setdefault("force_include", {})[str(so_path)] = "scanitd/ssw/ssw.so"
        build_data.setdefault("artifacts", []).append("src/scanitd/ssw/ssw.so")
