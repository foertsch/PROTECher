"""Shared test helpers.

Several tests need a stand-in for an external program - DockQ, FoldX - so they can
check what the pipeline passes to it without installing it. The obvious way to build
one is a file holding Python with a `#!/path/to/python` first line, marked
executable. That works on Linux and macOS, where the kernel reads the shebang, and
fails on Windows, which ignores it: `subprocess` hands the file to the OS, the OS
sees something that is not a .exe, and raises OSError [WinError 193].

`write_python_stub` papers over exactly that difference and nothing else.
"""

from __future__ import annotations

import os
import stat
import sys
from pathlib import Path


def write_python_stub(path: Path, body: str) -> Path:
    """Create an executable stand-in that runs `body` with this interpreter.

    Returns the path to invoke, which is NOT always `path`: on Windows the thing
    that can be executed is a .bat next to the script, so use the return value
    rather than the argument.

        exe = write_python_stub(tmp_path / "DockQ", "import sys; print(sys.argv)")
        subprocess.run([str(exe), "--flag"])
    """
    if os.name != "nt":
        path.write_text(f"#!{sys.executable}\n{body}", encoding="utf-8")
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        return path

    script = path.with_suffix(".py")
    script.write_text(body, encoding="utf-8")
    launcher = path.with_suffix(".bat")
    # @ hides the line itself; %* forwards every argument the caller passed.
    launcher.write_text(f'@"{sys.executable}" "{script}" %*\n', encoding="utf-8")
    return launcher
