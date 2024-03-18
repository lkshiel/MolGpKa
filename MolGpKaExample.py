#IMPORTANT: this only works if you are in the MolGpKa/src directory
#import MolGpKa
from src.predict_pka import *

#define MolGpKa function
def RunMolGpKa(x):
	mol = Chem.MolFromSmiles(x)
	base_dict, acid_dict = predict(mol)
	pkas = list(base_dict.values()) + list(acid_dict.values())
	sites=len(pkas)
	yield sites, pkas

#run test with list of examples
examples=['CC(O)=O','CC(C)C(N)C(O)=O','C(O)1=CC=C(N)C=C1','NC(CCS)C(O)=O','NC(CC1=CN=CN1)C(O)=O']
for smi in examples:
	data=RunMolGpKa(smi)
	for n,p in data:
		pkaSites=n
		pkaList=p
		print('SMILES:',smi,'\n','\t','# of sites:',n,'\n','\t','pKa:',p)

# #expected output
# SMILES: CC(O)=O 
# 	# of sites: 1 
# 	pKa: [8.337572]
# SMILES: CC(C)C(N)C(O)=O 
# 	# of sites: 2 
# 	pKa: [9.808557, 8.171646]
# SMILES: C(O)1=CC=C(N)C=C1 
# 	# of sites: 3 
# 	pKa: [5.06671, 14.116651, 12.194197]
# SMILES: NC(CCS)C(O)=O 
# 	# of sites: 3 
# 	pKa: [8.882975, 10.914888, 7.571865]
# SMILES: NC(CC1=CN=CN1)C(O)=O 
# 	# of sites: 4 
# 	pKa: [8.179463, 6.96088, 7.530961, 14.086538]