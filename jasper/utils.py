from pathlib import Path
import json


def chunks(lst: list, n: int):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def parse_config(config_path: str, default_config: dict) -> dict:
    try:
        with open(config_path, 'r') as fh:
            return json.loads(fh)
    except (FileNotFoundError, ValueError, TypeError):
        return default_config
