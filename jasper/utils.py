import json
from shutil import which

LOGO = r"""
      _   _    ____  ____  _____ ____  
     | | / \  / ___||  _ \| ____|  _ \ 
  _  | |/ _ \ \___ \| |_) |  _| | |_) |
 | |_| / ___ \ ___) |  __/| |___|  _ < 
  \___/_/   \_\____/|_|   |_____|_| \_\    

  ========================================

            |-----------------|
             \---------------/
              ~-_---------_-~
                 ~-_---_-~
                    ~-_
                 _-~---~-_
              _-~---------~-_
             /---------------\
            |-----------------|                               
"""

TYPES = ('fasta', 'fa', 'fna')


def chunks(lst: list, n: int):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def parse_config(config_path: str, default_config: dict) -> dict:
    try:
        with open(config_path, 'r') as fh:
            return json.load(fh)
    except (FileNotFoundError, ValueError, TypeError) as e:
        return default_config


def perform_tool_check(name: str):
    """Check whether `name` is on PATH and marked as executable.

    `name` is a tool which JASPER will use, for instance `blastn`.

    Args:
        name (str): Name of tool to perform check on.
    Returns:
        (bool): Value indicating whether tool exists in PATH or not.
    """
    return which(name) is not None
