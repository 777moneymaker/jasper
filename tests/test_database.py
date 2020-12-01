import unittest
from unittest.mock import MagicMock

from jasper import database


class DataBaseTests(unittest.TestCase):
    def test_concat_fasta_files(self):
        cmd = database.concat_fasta_files_cmd("data/fasta_test_data")
        final = "find /home/mlchodkowski/Desktop/JASPER/tests/data/fasta_test_data/ \( -iname \"*.fasta\" -o -iname \"*.fa\" -o -iname \"*.fna\" \) -exec cat {} \; > temp.fasta"
        self.assertEqual(cmd, final)

        cmd = database.concat_fasta_files_cmd("data/fasta_test_data", "test")
        final = "find /home/mlchodkowski/Desktop/JASPER/tests/data/fasta_test_data/ \( -iname \"*.fasta\" -o -iname \"*.fa\" -o -iname \"*.fna\" \) -exec cat {} \; > test"
        self.assertEqual(cmd, final)


if __name__ == '__main__':
    unittest.main()
