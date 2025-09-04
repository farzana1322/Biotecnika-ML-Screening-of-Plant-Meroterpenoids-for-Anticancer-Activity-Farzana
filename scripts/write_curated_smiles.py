import pandas as pd

# Define your curated SMILES
data = {
    'compound_id': ['CMPD001', 'CMPD002', 'CMPD003', 'CMPD004', 'CMPD005', 'CMPD006'],
    'SMILES': [
        'CC(C)C1=CC=CC=C1C(=O)O',
        'COC1=CC=CC=C1C(=O)O',
        'CCC(=O)OC1=CC=CC=C1',
        'CC1=CC=CC=C1C(=O)OC',
        'CC(C)OC(=O)C1=CC=CC=C1',
        'C[C@H](O)C(=O)OC1=CC=CC=C1'
    ]
}

df = pd.DataFrame(data)
df.to_csv('../data/curated_smiles.csv', index=False)
print("✅ curated_smiles.csv updated with all 6 compounds")
