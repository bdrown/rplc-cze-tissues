#!/user/bin/env python
"""calc_protein_props.py

This script allow the user to calculate interesting properties of proteins and
proteoforms that were discovered by TDportal.

One command line argument is required:
  - path to dump of hits from TDportal search
  - path to output file
  
Usage:
   python calc_protein_props.py path/to/hits.csv path/to/output.csv
   
"""

__author__ = "Bryon Drown"
__version__ = "1.0.0"
__maintainer__ = "Bryon Drown"
__email__ = "bryon.drown@gmail.com"

# Imports
import sys
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import progressbar

# Grab commandline arguments
# Decompress gzipped XML
source = sys.argv[1]
output = open(sys.argv[2], 'w')

print(source)
df = pd.read_csv(source, sep=',', dtype=str)
header = ['HitId','PFR','Uniprot_Id','Accession','Description','Sequence','Start_Index','End_Index','Sequence_Length',
          'Modifications','Modification_Codes','N_Terminal_Modification_Code','C_Terminal_Modification_Code',
          'Monoisotopic_Mass','Average_Mass','File_Name','Result_Set','Observed_Precursor_Mass','Precursor_Scans',
          'Fragment_Scans','Average_Retention_Time','Global_Q-value','P-score','E-value','C-score','Percent_Cleavages','Proteoform_Level']
df.columns = header

props = pd.DataFrame(columns=['GRAVY', 'pI', 'Aromaticity', 'Instability'])
for seq in progressbar.progressbar(df.Sequence):
    try:
        X = ProteinAnalysis(seq)
        props = props.append({'GRAVY':X.gravy(), 'pI':X.isoelectric_point(), 'ChargeAt24':X.charge_at_pH(2.4), 'Aromaticity':X.aromaticity(), 'Instability':X.instability_index()}, ignore_index=True)
    except KeyError:
        print("Unknown amino acid found")
        print("Sequence: %s" % seq)
        props = props.append({'GRAVY':'NA', 'pI':'NA', 'ChargeAt24':'NA', 'Aromaticity':'NA', 'Instability':'NA'}, ignore_index=True)
        continue

df = pd.concat([df, props], axis=1)
df.to_csv(output, sep=',', header=True, index=False)