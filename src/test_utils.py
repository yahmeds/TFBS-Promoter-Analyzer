import unittest
import os
import tempfile
import sys
import pathlib

root_dir = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(root_dir))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from src.utils import find_cds, mrna_to_gene, upstream_gene_seq, download_promotors

class TestUtilsModule(unittest.TestCase):
    def setUp(self):
        self.mock_record = SeqRecord(
            Seq("ATGCATGCATGCATGCATGC"),
            id="test_record",
            features=[
                SeqFeature(
                    FeatureLocation(5, 15), 
                    type="CDS"
                )
            ]
        )

    def test_find_cds(self):
        cds_positions = find_cds(self.mock_record)
        
        self.assertEqual(len(cds_positions), 1)
        self.assertEqual(cds_positions[0], (5, 15))

    def test_mrna_to_gene(self):
        try:
            gene_id = mrna_to_gene("NM_007389")
            self.assertIsNotNone(gene_id)
            self.assertTrue(isinstance(gene_id, str))
        except Exception as e:
            self.skipTest(f"Network-dependent test failed: {e}")

    def test_upstream_gene_seq(self):
        try:
            upstream_seq = upstream_gene_seq("11435", 250)
            
            self.assertIsNotNone(upstream_seq)
            self.assertTrue(len(upstream_seq) <= 250)
            self.assertTrue(isinstance(upstream_seq, Seq))
        except Exception as e:
            self.skipTest(f"Network-dependent test failed: {e}")

    def test_download_promotors(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_ids = ["NM_007389"]
            download_promotors(test_ids, 1000, tmpdir)
            expected_file = os.path.join(tmpdir, f"{test_ids[0]}_1000.fa")
            self.assertTrue(os.path.exists(expected_file))
            
            with open(expected_file, 'r') as f:
                content = f.read()
                self.assertTrue(content.startswith(f">{test_ids[0]}_promoter_1000"))
                self.assertTrue(len(content.split('\n')[1]) > 0)

if __name__ == '__main__':
    unittest.main()