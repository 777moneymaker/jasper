import unittest
from jasper import io


class IOTests(unittest.TestCase):
    def test_retrieve(self):
        ids = [seq.id for seq in io.retrieve('data/fasta_test_data/')]
        self.assertEqual(ids, ["seq1", "seq2"])

    def test_retrieve_isdir(self):
        with self.assertRaises(IsADirectoryError):
            list(io.retrieve("data/fasta_test_data/seqs.fasta"))

    def test_retrieve_isempty(self):
        with self.assertRaises(OSError):
            list(io.retrieve("data/empty_test_dir"))


if __name__ == '__main__':
    unittest.main()
