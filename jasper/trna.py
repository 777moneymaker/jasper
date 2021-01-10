from pathlib import Path

from jasper import database


class tRNAScanner(database.Database):
    def __init__(self, source_dir: Path, name: str) -> None:
        super(tRNAScanner, self).__init__(source_dir, name)
