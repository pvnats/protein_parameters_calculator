from decimal import *
import unittest
from Tk_ProteinParameters import Protein


class TestCalculations(unittest.TestCase):
    test_sequence = "ARNDCQEGHILKMFPSTWYVARNDCQEGHILKMFPSTWYV"
    
    @classmethod
    def setUpClass(self):
        self.protein = Protein()
        self.protein.calculate(self.test_sequence)
        
    def test_make_sequence(self):
        expected_sequence =['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        self.assertListEqual(expected_sequence, self.protein.sequence, 'The sequence is wrong')
        
    def test_count_aa(self):
        expected_aa_dict = {'A': 2, 'R': 2, 'N': 2, 'D': 2, 'C': 2, 'Q': 2, 'E': 2, 'G': 2, 'H': 2, 'I': 2, 'L': 2, 'K': 2, 'M': 2, 'F': 2, 'P': 2, 'S': 2, 'T': 2, 'W': 2, 'Y': 2, 'V': 2}
        self.assertDictEqual(expected_aa_dict, self.protein.aa_dict, 'The aa_dict is wrong')
        
    def test_aa_number(self):
        expected_aa_number = 40
        self.assertEqual(expected_aa_number, self.protein.aa_number, 'The aa_number is wrong')
    
    def test_count_atoms(self):
        expected_atom_dict = {'C': 214, 'H': 314, 'N': 58, 'O': 58, 'S': 4}
        self.assertDictEqual(expected_atom_dict, self.protein.atom_dict, 'The atom_dict is wrong')
        
    def test_calculate_mol_weight(self):
        expected_mol_weight = Decimal('4773.45410')
        self.assertEqual(expected_mol_weight, self.protein.mol_weight, 'The mol_weight is wrong')
    
    def test_calculate_ext_coef(self):
        expected_ext_coef = {
            "ext_coef_cformed": 14105, 
            "ext_coef_reduced": 13980, 
            "abs_cformed": Decimal('2.954883341184740835781787448'),
            "abs_reduced": Decimal('2.928696852872220977258375649')
            }
        self.assertDictEqual(expected_ext_coef, self.protein.ext_coef, 'The ext_coef dict is wrong')
        
    def test_calculate_pI(self):
        expected_pI = Decimal('6.9541015625')
        self.assertEqual(expected_pI, self.protein.pI, 'The pI is wrong')
        
if __name__ == '__main__':
    unittest.main()