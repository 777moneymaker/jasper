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


def parse_config(config_path: str, default_config: dict) -> dict:
    """This function parses given config json file and returns default if any error

    Args:
        config_path (str): Path to json file with configuration.
        default_config (dict): Default config to be returned if any error.
    Returns:
        (dict): Parsed config.
    """
    try:
        with open(config_path, 'r') as fh:
            return json.load(fh)
    except (FileNotFoundError, ValueError, TypeError) as e:
        return default_config


def perform_tool_check(names: list):
    """Check whether `name` is on PATH and marked as executable.

    `name` is a tool which JASPER will use, for instance `blastn`.

    Args:
        names (list): Names of tools to perform check on.
    Returns:
        (bool): Value indicating whether all tools exists in PATH or not.
    """
    return all(which(x) is not None for x in names)
