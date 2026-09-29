"""Regression checks for protein sequence preservation and species labels.

Run with: python3 -m unittest discover -s tests -p 'test_protein_headers.py'
"""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


BIN = Path(__file__).resolve().parents[1] / 'bin'


class ProteinHeaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def annotate(self, consequence, hgvsp, amino_acids, position='2/7',
                 translated='NQQQQ', reference='MKKKKKK', source='ANN',
                 organism='Homo sapiens', taxonomy='9606', transcript='TX1'):
        if source == 'ANN':
            fields = ['AT', consequence, 'HIGH', 'GENE', 'G1', 'transcript',
                      'TX1', 'protein_coding', '1/1', 'c.4dup', hgvsp,
                      '4/21', '4/21', position, '', '']
        else:
            fields = ['AT', consequence, 'HIGH', 'GENE', 'G1', 'Transcript',
                      'TX1', 'protein_coding', '1/1', '', 'c.4dup', hgvsp,
                      '4/21', '4/21', position, amino_acids]
        vcf = self.root / 'variants.vcf'
        vcf.write_text('##fileformat=VCFv4.2\n'
                       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                       'chr1\t10\t.\tA\tAT\t.\tPASS\t'
                       + source + '=' + '|'.join(fields) + '\n')
        fasta = self.root / 'translated.fa'
        fasta.write_text('>var_._chr1.10.A.AT_TX1\n' + translated + '\n')
        ref = self.root / 'reference.fa'
        ref.write_text('>' + transcript + '\n' + reference + '\n')
        out = self.root / 'out.fa'
        subprocess.run([sys.executable, str(BIN / 'get_var_aa_change.py'),
                        str(vcf), str(fasta), str(ref), str(out),
                        '--organism', organism, '--taxonomy-id', taxonomy],
                       check=True, capture_output=True, text=True)
        header, *sequence = out.read_text().splitlines()
        return header, ''.join(sequence)

    def test_indels_and_frameshifts_preserve_translation(self):
        cases = [
            ('frameshift_variant', 'p.Lys2fs', 'K/X'),
            ('inframe_deletion', 'p.Lys2del', 'K/-'),
            ('inframe_insertion', 'p.Lys2_Gln3insAla', '-/A'),
            ('inframe_insertion', 'p.Lys2dup', '-/K'),
            ('missense_variant&frameshift_variant', 'p.Lys2Asn', 'K/N'),
        ]
        for source in ['ANN', 'CSQ']:
            for consequence, hgvsp, amino_acids in cases:
                with self.subTest(source=source, consequence=consequence, hgvsp=hgvsp):
                    _, sequence = self.annotate(consequence, hgvsp, amino_acids, source=source)
                    self.assertEqual(sequence, 'NQQQQ')

    def test_verified_missense_correction(self):
        for source in ['ANN', 'CSQ']:
            with self.subTest(source=source):
                _, sequence = self.annotate('missense_variant', 'p.Lys2Asn', 'K/N', source=source)
                self.assertEqual(sequence, 'MNKKKKK')

    def test_unverifiable_missense_preserves_translation(self):
        cases = [
            {'reference': 'MAKKKKK'},
            {'transcript': 'OTHER'},
            {'position': '0/7'},
            {'position': '8/7'},
            {'position': '2-3/7'},
            {'position': ''},
            {'hgvsp': 'p.?', 'amino_acids': '-/-'},
        ]
        for case in cases:
            with self.subTest(case=case):
                args = dict(consequence='missense_variant', hgvsp='p.Lys2Asn', amino_acids='K/N')
                args.update(case)
                _, sequence = self.annotate(**args)
                self.assertEqual(sequence, 'NQQQQ')

    def test_complete_translation_is_preserved(self):
        _, sequence = self.annotate('missense_variant', 'p.Lys2Asn', 'K/N', translated='MNQQQQ')
        self.assertEqual(sequence, 'MNQQQQ')

    def test_species_labels_in_both_scripts(self):
        for organism, taxonomy in [('Homo sapiens', '9606'),
                                   ('Mus musculus', '10090'),
                                   ('Macaca fascicularis', '9541')]:
            with self.subTest(organism=organism):
                header, _ = self.annotate('missense_variant', 'p.Lys2Asn', 'K/N',
                                          organism=organism, taxonomy=taxonomy)
                self.assertIn(f'OS={organism} OX={taxonomy}', header)
                inp, out = self.root / 'novel.fa', self.root / 'novel.out.fa'
                inp.write_text('>STRG.1.1.p1\nMQQQQ\n')
                subprocess.run([sys.executable, str(BIN / 'reformat_novel_protein_headers.py'),
                                '--input', str(inp), '--output', str(out), '--sample-id', 'sample',
                                '--organism', organism, '--taxonomy-id', taxonomy],
                               check=True, capture_output=True, text=True)
                self.assertIn(f'OS={organism} OX={taxonomy}', out.read_text())
                self.assertEqual(out.read_text().splitlines()[1], 'MQQQQ')

    def test_standalone_scripts_require_species(self):
        commands = [
            ['get_var_aa_change.py', 'a.vcf', 'b.fa', 'c.fa', 'out.fa'],
            ['reformat_novel_protein_headers.py', '--input', 'a.fa', '--output', 'out.fa',
             '--sample-id', 'sample'],
        ]
        for script, *args in commands:
            with self.subTest(script=script):
                result = subprocess.run([sys.executable, str(BIN / script), *args],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 2)
                self.assertIn('--organism', result.stderr)
                self.assertIn('--taxonomy-id', result.stderr)


if __name__ == '__main__':
    unittest.main()
