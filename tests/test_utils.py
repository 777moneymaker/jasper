#!/usr/bin/env python3
import unittest
from pathlib import Path

from jasper import utils


class UtilsTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(UtilsTests, self).__init__(*args, **kwargs)
        self.test_dir = Path(__file__).parent.absolute()

    def test_parse_config_ok(self):
        res = utils.parse_config(str(self.test_dir / Path("data/test_config.json")), {"test": 1})
        self.assertEqual(res, {"test": "ok", "value": 1})

    def test_parse_config_default(self):
        res = utils.parse_config(str(self.test_dir / Path("data/not_existing.json")), {"test": 1})
        self.assertEqual(res, {"test": 1})

    def test_perform_tools_check(self):
        self.assertTrue(utils.perform_tool_check(['makeblastdb', 'blastn', 'pilercr', 'WIsH', 'tRNAscan-SE']))


if __name__ == '__main__':
    unittest.main()
