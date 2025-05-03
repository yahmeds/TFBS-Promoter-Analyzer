import unittest
import os
import json
import sys
import pathlib
import argparse

root_dir = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(root_dir))

from Bio import motifs
from Bio.Seq import Seq

from src.putativeTFBS import json_serializable, get_promoter_sequences, process_tfbs, parse_arguments, format_results

class TestPutativeTFBSModule(unittest.TestCase):
    def setUp(self):
        self.test_mrna_ids = ["NM_000451", "NM_007389"]
        self.test_pfm_file = "./data/MA0083.3.jaspar"
        
        self.mock_args = argparse.Namespace(
            pfm=self.test_pfm_file,
            mrna_ids=self.test_mrna_ids,
            threshold=-20,
            promotor_length=500,
            window_size=40,
            window_threshold=None,
            pseudocount=0.1,
            output=None
        )

    def test_json_serializable(self):
        import numpy as np

        np_int = np.int32(42)
        self.assertEqual(json_serializable(np_int), 42)
        
        np_float = np.float64(3.14)
        self.assertEqual(json_serializable(np_float), 3.14)
        
        np_array = np.array([1, 2, 3])
        self.assertEqual(json_serializable(np_array), [1, 2, 3])

    def test_get_promoter_sequences(self):
        promoter_sequences = get_promoter_sequences(self.test_mrna_ids, 500)
        
        self.assertEqual(len(promoter_sequences), len(self.test_mrna_ids))
        for seq in promoter_sequences:
            self.assertTrue(len(seq) <= 500)
            self.assertTrue(hasattr(seq, 'id'))

    def test_process_tfbs(self):
        promoter_sequences = get_promoter_sequences(self.test_mrna_ids, 500)
        
        results = process_tfbs(self.mock_args, promoter_sequences)
        
        self.assertIn('parameters', results)
        self.assertIn('matches', results)
        self.assertIn('windows', results)
        
        self.assertEqual(results['parameters']['mrna_ids'], self.test_mrna_ids)

    def test_parse_arguments(self):
        import sys
        original_argv = sys.argv
        
        try:
            sys.argv = [
                'script_name', 
                '-m', './data/MA0083.3.jaspar', 
                '-t', '-20', 
                'NM_000451', 'NM_007389'
            ]
            
            args = parse_arguments()
            
            self.assertEqual(args.pfm, './data/MA0083.3.jaspar')
            self.assertEqual(args.threshold, -20)
            self.assertEqual(args.mrna_ids, ['NM_000451', 'NM_007389'])
        finally:
            sys.argv = original_argv

    def test_format_results(self):
        mock_results = {
            'matches': [
                {'mrna_id': 'NM_000451', 'position': -100, 'score': -15, 'motif': 'Test'},
                {'mrna_id': 'NM_007389', 'position': -200, 'score': -18, 'motif': 'Test'}
            ],
            'windows': [
                {
                    'number': 1,
                    'motif': 'Test',
                    'range': '[-200:-160]',
                    'score': -16.5,
                    'match_count': 2,
                    'matches': [
                        {'mrna_id': 'NM_000451', 'position': -100, 'score': -15, 'motif': 'Test'},
                        {'mrna_id': 'NM_007389', 'position': -200, 'score': -18, 'motif': 'Test'}
                    ]
                }
            ]
        }
        
        import io
        import sys
        
        captured_output = io.StringIO()
        sys.stdout = captured_output
        
        format_results(mock_results)
        
        sys.stdout = sys.__stdout__
        
        output = captured_output.getvalue()
        self.assertIn('1 Test [-200:-160]', output)
        self.assertIn('NM_000451', output)
        self.assertIn('NM_007389', output)

if __name__ == '__main__':
    unittest.main()