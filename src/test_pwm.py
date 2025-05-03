import unittest
import sys
import pathlib

root_dir = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(root_dir))

from Bio import motifs
from Bio.Seq import Seq
import numpy as np

from src.pwm import pwm2pssm, scan_sequence, scan_all_sequences, score_window, best_window

class TestPWMModule(unittest.TestCase):
    def setUp(self):
        self.sample_matrix = motifs.Motif()
        self.sample_matrix.name = "TestMotif"
        
        self.sample_matrix.counts = {
            'A': [40, 10, 20, 30],
            'C': [30, 40, 10, 20],
            'G': [20, 30, 40, 10],
            'T': [10, 20, 30, 40]
        }
        
        self.test_sequence = Seq("ACGTACGTACGT")

    def test_pwm2pssm(self):
        pssm = pwm2pssm(motifs.matrix.FrequencyPositionMatrix(
            ['A', 'C', 'G', 'T'], 
            self.sample_matrix.counts
        ))
        self.assertIsNotNone(pssm)
        self.assertEqual(pssm.length, 4) 

    def test_scan_sequence(self):
        pssm = pwm2pssm(motifs.matrix.FrequencyPositionMatrix(
            ['A', 'C', 'G', 'T'], 
            self.sample_matrix.counts
        ))
        
        results = scan_sequence(pssm, self.test_sequence, threshold=-10)
        
        self.assertIsInstance(results, list)
        for result in results:
            self.assertIsInstance(result[0], (int, float))  
            self.assertIsInstance(result[1], float) 

    def test_scan_all_sequences(self):
        pssm = pwm2pssm(motifs.matrix.FrequencyPositionMatrix(
            ['A', 'C', 'G', 'T'], 
            self.sample_matrix.counts
        ))
        sequences = [self.test_sequence, 
                     Seq("TGCATGCATGCA"),
                     Seq("GTACGTACGTAC")]
        
        results = scan_all_sequences(pssm, sequences, threshold=-10)
        
        self.assertEqual(len(results), len(sequences))
        for seq_results in results:
            self.assertIsInstance(seq_results, list)

    def test_score_window(self):
        pssm = pwm2pssm(motifs.matrix.FrequencyPositionMatrix(
            ['A', 'C', 'G', 'T'], 
            self.sample_matrix.counts
        ))
        sequences = [self.test_sequence]
        scan_results = scan_all_sequences(pssm, sequences, threshold=-10)
        
        window_score = score_window(scan_results, -200, -100)
        self.assertIsInstance(window_score, float)

    def test_best_window(self):
        pssm = pwm2pssm(motifs.matrix.FrequencyPositionMatrix(
            ['A', 'C', 'G', 'T'], 
            self.sample_matrix.counts
        ))
        sequences = [self.test_sequence, 
                     Seq("TGCATGCATGCA"),
                     Seq("GTACGTACGTAC")]
        scan_results = scan_all_sequences(pssm, sequences, threshold=-10)
        
        best_start, best_end, best_score = best_window(scan_results, window_size=40)
        
        if best_start is not None:
            self.assertTrue(isinstance(best_start, (int, float)) or best_start is None)
        if best_end is not None:
            self.assertTrue(isinstance(best_end, (int, float)) or best_end is None)
        self.assertIsInstance(best_score, float)

if __name__ == '__main__':
    unittest.main()