import logging

from src.predict_pka import *


class CTSMolgpka:

	def __init__(self):
		pass

	def convert_floats(self, pka_list):
		"""
		Python's json serializer doesn't like np.float32 types, so
		converting them into python floats.
		"""
		return [float(i) for i in pka_list]

	def run_molgpka(self, smiles):
		mol = Chem.MolFromSmiles(smiles)
		base_dict, acid_dict = predict(mol)
		pkas = list(base_dict.values()) + list(acid_dict.values())
		sites=len(pkas)
		yield sites, pkas

	def main(self, smiles):
		"""
		Main function for returning pkas and/or microspecies.
		Examples=['CC(O)=O','CC(C)C(N)C(O)=O','C(O)1=CC=C(N)C=C1','NC(CCS)C(O)=O','NC(CC1=CN=CN1)C(O)=O']
		"""
		pka_sites, pka_list = None, None

		data = self.run_molgpka(smiles)
		
		for n,p in data:
			pka_sites=n
			pka_list=p
			print('SMILES:',smiles,'\n','\t','# of sites:',n,'\n','\t','pKa:',p)

		pka_list = self.convert_floats(pka_list)

		return smiles, pka_sites, pka_list



if __name__ == "__main__":
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

	examples = ['CC(O)=O','CC(C)C(N)C(O)=O','C(O)1=CC=C(N)C=C1','NC(CCS)C(O)=O','NC(CC1=CN=CN1)C(O)=O']

	molgpka = CTSMolgpka()

	for smi in examples:
		data = molgpka.main(smi)
