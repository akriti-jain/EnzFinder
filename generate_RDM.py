'''
## Copyright Notice
EnzFinder code repository is a TCS proprietary resource and should be used for academic purposes only. The contents of this repository should not be used for any commercial purpose without the consent of ALL the authors involved. By downloading and utilizing the scripts, the user consents that any and all Intellectual Property derived from the EnzFinder code repository is fully owned by TCS in the associated jurisdictions. EnzFinder code repository usage without citation will be considered illegal.
'''
import re
import kcfconvoy as kcf
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.rdmolops import GetDistanceMatrix
RDLogger.DisableLog('rdApp.*')
from round1 import sortRDM

def generate_RDM(reactant, product):
	'''
	RDM is generated for query reaction
	Input : mapped reactant SMILES, mapped Product SMILES
	
	A reaction can be multi-substrate or multi-product reaction where a reactant pair are not similar. Therefore, in the interest of time, 
	we first check if there is any common mapped position between reactant and product. If yes, reaction center is searched.

	mapped reactant SMILES and mapped product SMILES are converted into canonical SMILES. 
	Using canonical SMILES, KCF are looked for each atom in a molecule.
	Also, dictionary of index to mapped-position, mapped-position-KCF, mapped-position-neighbors, etc are created.
	
	identify_difference function is called to identify reaction center and generate RDM

	Output : RDM, reaction-center-mapped-position, farthest-atom of both reactant and product, reactant pair

	'''
	reaction_RDM, reaction_KCF_mapped_position, reaction_farthest_atom, reactant_product_pair = [],{},{},{}
	for react_one in reactant.split(".") :	
		for product_one in product.split(".") :
			reactant_mappedNumber = find_mapped_number_in_molecule(react_one)
			product_mappedNumber = find_mapped_number_in_molecule(product_one)

			if set(reactant_mappedNumber).intersection(set(product_mappedNumber)) :
						reactant_FN_KCF, reactant_details = smiesTocanonical(react_one)		
						product_FN_KCF, product_details = smiesTocanonical(product_one)
						
						result_RDM_to_be_saved,KCF_mapped_position,farthest_atom = identify_difference(reactant_details, product_details, reactant_FN_KCF,product_FN_KCF,react_one,product_one)
						if len(result_RDM_to_be_saved)>0 :
							reaction_RDM.append(result_RDM_to_be_saved)
							for i in range(len(result_RDM_to_be_saved)) :
								temp_rdm = result_RDM_to_be_saved[i]
								temp_sorted_rdm = sortRDM(temp_rdm)
								if temp_rdm.split(":")[0]==temp_sorted_rdm.split(":")[0] :
									reactant_product_pair[result_RDM_to_be_saved[i]]=react_one+">>"+product_one
								else :
									reactant_product_pair[result_RDM_to_be_saved[i]]=product_one+">>"+react_one
								reaction_KCF_mapped_position[result_RDM_to_be_saved[i]]=KCF_mapped_position[i]
								reaction_farthest_atom[result_RDM_to_be_saved[i]]=farthest_atom[i]

	return (reaction_RDM, reaction_KCF_mapped_position, reaction_farthest_atom,reactant_product_pair)

def find_mapped_number_in_molecule (Smiles) :
	'''
	List of mapped number in a SMILES
	Input : SMILES
	Return : list of mapped number of each atom in a molecule

	'''
	mol = Chem.MolFromSmiles(Smiles)
	mappedNumber = []
	for i,atom in enumerate(mol.GetAtoms()):
		mappedNumber.append(atom.GetAtomMapNum())
	return mappedNumber

def smiesTocanonical(Input_Smiles):
	'''
	SMILES are converted to canonical SMILES
	Input : mapped SMILES
	canonical mapped SMILES is used to call kcfv() function
	Return : results from kcfv() function
	'''
	can = Chem.MolToSmiles(Chem.MolFromSmiles(Input_Smiles))
	return kcfv(can)

def kcfv(smiles) :
	'''
	KCF is generated for each atom of a molecule using kcfconvoy
	Input : canonical mapped SMILES
	kcfextension() function and find_neighbors() function are called
	Return : result from kcfextension() function and find_neighbors() function 
	'''

	k = kcf.KCFvec()
	k.input_smiles(smiles)
	k.convert_kcf_vec()
	KCF_index_mapList = {key: value["kegg_atom"] for key, value in k.kegg_atom_label.items()} #atom_index:atom_kcf
	return kcfextension(smiles,KCF_index_mapList),find_neighbors(smiles,KCF_index_mapList)

def kcfextension(SMILES,KCF_index_mapList) :
	'''
	to identify reaction center that have same KCF but neighbor atoms are change, KCF is added with atom symbol.
	Input : reaction SMILES (either reactant or product) and list of KCF atom at every index
	Return : KCF_mapped_subs
	KCF_mapped_subs ={atom-mapped-position:KCF_neighbor-atoms}
	'''
	mol = Chem.MolFromSmiles(SMILES)
	KCF_mapped_subs ={}
	for i,atom in enumerate(mol.GetAtoms()): # take an atom in molecule
		val=[]
		val = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
		val.sort()
		KCF_mapped_subs[atom.GetAtomMapNum()]=KCF_index_mapList[i]+"_"+''.join(val)
	return KCF_mapped_subs



def find_neighbors(Smiles,KCFmapping) :
	'''
	Find all neighbours of all atoms in a molecule.
	Input : SMILES, dict of KCF and mapped number
	Output : list of mappedNumber_to_index,index_to_mappedNUmber,mappedNumber_to_atom,index_to_atom,index_to_neighbors,mappedNumber_to_KCF
	
	mappedNumber_to_index = {mapped-number:atom-index}
	index_to_mappedNUmber = (atom-index:mapped-number)
	mappedNumber_to_atom = {mapped-number:atom-symbol}
	index_to_atom = {atom-index:atom-symbol}
	index_to_neighbors = {atom-index:(neighbor-atom-index,neighbor-atom-symbol)}
	mappedNumber_to_KCF = {mapped-number:KCF}
	
	'''
	mol = Chem.MolFromSmiles(Smiles)
	mol = Chem.AddHs(mol)
	mappedNumber_to_index ={}
	index_to_mappedNUmber = {}
	mappedNumber_to_atom = {}
	index_to_atom = {}
	index_to_neighbors = {}
	mappedNumber_to_KCF = {}
	for i,atom in enumerate(mol.GetAtoms()):
		mappedNumber_to_index[atom.GetAtomMapNum()]= atom.GetIdx()
		index_to_mappedNUmber[atom.GetIdx()] = atom.GetAtomMapNum()
		mappedNumber_to_atom[atom.GetAtomMapNum()] = atom.GetSymbol()
		index_to_atom[atom.GetIdx()] = atom.GetSymbol()
		index_to_neighbors[atom.GetIdx()]= [(nbr.GetIdx(),nbr.GetSymbol()) for nbr in atom.GetNeighbors()]
		if atom.GetSymbol() != 'H' :
			mappedNumber_to_KCF[atom.GetAtomMapNum()]= KCFmapping[atom.GetIdx()]
	
	return ([mappedNumber_to_index,index_to_mappedNUmber,mappedNumber_to_atom,index_to_atom,index_to_neighbors,mappedNumber_to_KCF])



def calculate_farthest_atom_distance_from_reaction_centre(mol_smiles, kcfPos_index) :
	'''
	Calculate the radius of a molecule considering reaction center as the center
	Input : SMILES
	Return : radius of the molecule 
	'''
	mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(mol_smiles)))	
	return round(max(GetDistanceMatrix(mol)[kcfPos_index]))


def identify_difference(reactant_details, product_details,react_fn_kcf,prod_fn_kcf,r_smiles,p_smiles):
	'''
	Identify the reaction center (R), find different atom (D) and matched atom (M)
	For same mapped number, if KCF is different, then it is R
	If R is not false positive KCF, then find D and M for each neighbor atom position.
	Once all RDMs are generated for a reactant pair, then check if any RDM is overlapping RDM. If yes, then remove the RDM

	Input :  reactant_details, product_details,react_fn_kcf,prod_fn_kcf,r_smiles,p_smiles
	Output : list of RDMs, list of R position, list of radius for each R 

	'''
	temp_kcf_mapped_pos = {}
	temp_farthest_atom_react_prod = {}
	finalRDMlist = []
	reactant_mappedNumber_to_index,reactant_index_to_mappedNUmber,reactant_mappedNumber_to_atom,reactant_index_to_atom,reactant_index_to_neighbors,reactant_mappedNumber_to_KCF = reactant_details[0],reactant_details[1], reactant_details[2], reactant_details[3], reactant_details[4],reactant_details[5]
	product_mappedNumber_to_index,product_index_to_mappedNUmber,product_mappedNumber_to_atom,product_index_to_atom,product_index_to_neighbors,product_mappedNumber_to_KCF = product_details[0],product_details[1], product_details[2], product_details[3], product_details[4],product_details[5]
	
	for kcfMappedpos in set(react_fn_kcf).intersection(set(prod_fn_kcf)) :
			reaction_center, dissimilar, match ="","",""
			temp_reactant_sub_str, temp_product_sub_str =[],[]
			if (react_fn_kcf[kcfMappedpos] != prod_fn_kcf[kcfMappedpos]) :
				if (falsePositiveCheck((reactant_mappedNumber_to_KCF[kcfMappedpos],product_mappedNumber_to_KCF[kcfMappedpos]))) :
					continue
				elif (falsePositiveCheck((product_mappedNumber_to_KCF[kcfMappedpos],reactant_mappedNumber_to_KCF[kcfMappedpos]))) :
					continue
				else :
					reactant_neighbors = reactant_index_to_neighbors[reactant_mappedNumber_to_index[kcfMappedpos]]
					product_neighbors = product_index_to_neighbors[product_mappedNumber_to_index[kcfMappedpos]]
					
					reaction_center = tuple((reactant_mappedNumber_to_KCF[kcfMappedpos], product_mappedNumber_to_KCF[kcfMappedpos]))
					dissimilar, match = find_d_m(reactant_neighbors, product_neighbors, reactant_index_to_mappedNUmber, product_index_to_mappedNUmber,reactant_mappedNumber_to_KCF,product_mappedNumber_to_KCF)
					if len(reaction_center)>0 :
						reaction_center, dissimilar, match = completeRDM(reaction_center, dissimilar, match)
						# if dissimilar=='*-*' :
						# 	print ("No atom in D", (reaction_center, dissimilar, match))
						
						if (reaction_center+":"+dissimilar+":"+match) not in finalRDMlist:
							# sorted_rdm = sortRDM(reaction_center+":"+dissimilar+":"+match)
							sorted_rdm = reaction_center+":"+dissimilar+":"+match
							finalRDMlist.append(sorted_rdm)
							temp_kcf_mapped_pos[sorted_rdm]=kcfMappedpos
							reactant_farthest_atom = calculate_farthest_atom_distance_from_reaction_centre(r_smiles,reactant_mappedNumber_to_index[kcfMappedpos])
							product_farthest_atom = calculate_farthest_atom_distance_from_reaction_centre(p_smiles,product_mappedNumber_to_index[kcfMappedpos])
							temp_farthest_atom_react_prod[sorted_rdm]=','.join([str(reactant_farthest_atom),str(product_farthest_atom)])
							
	
	# This is for checking Overlapping RDM
	finalRDMlist2 = tuple(finalRDMlist)
	for rdm_ind in list(range(len(finalRDMlist))) :
		RDM_rc = finalRDMlist[rdm_ind].split(":")[0].split("-")
		
		if RDM_rc[0][0]=='O' :
			if overlappingRDM(RDM_rc[0]) :
				RDM_d = finalRDMlist[rdm_ind].split(":")[1].split("-")[0]
				RDM_m = finalRDMlist[rdm_ind].split(":")[2].split("-")[0]
				temp = RDM_rc[0].replace('O','C')
				if temp in RDM_d :
					for i in finalRDMlist :
						j=i.split(":")[0].split("-")[0]
						if temp in j :
							finalRDMlist2 = list(finalRDMlist2)
							finalRDMlist2.remove(finalRDMlist[rdm_ind])
							finalRDMlist2 = tuple(finalRDMlist2)
							break
				elif temp in RDM_m :
					for i in finalRDMlist :
						j=i.split(":")[0].split("-")[0]
						
						if temp in j :
							finalRDMlist2 = list(finalRDMlist2)
							finalRDMlist2.remove(finalRDMlist[rdm_ind])
							finalRDMlist2 = tuple(finalRDMlist2)
							break
			
			elif RDM_rc[1][0]=='O' :
				if overlappingRDM(RDM_rc[1]) :
					RDM_d = finalRDMlist[rdm_ind].split(":")[1].split("-")[1]
					RDM_m = finalRDMlist[rdm_ind].split(":")[2].split("-")[1]
					temp = RDM_rc[1].replace('O','C')
					
					if temp in RDM_d :
						for i in finalRDMlist :
							j=i.split(":")[0].split("-")[1]
							
							if temp in j :
								finalRDMlist2 = list(finalRDMlist2)
								finalRDMlist2.remove(finalRDMlist[rdm_ind])
								finalRDMlist2 = tuple(finalRDMlist2)
								break
					elif temp in RDM_m :
						for i in finalRDMlist :
							j=i.split(":")[0].split("-")[1]
							if temp in j :
								finalRDMlist2 = list(finalRDMlist2)
								finalRDMlist2.remove(finalRDMlist[rdm_ind])
								finalRDMlist2 = tuple(finalRDMlist2)
								break

	pos, dist = [],[]
	for var in finalRDMlist2 :
		pos.append(temp_kcf_mapped_pos[var])
		dist.append(temp_farthest_atom_react_prod[var])
	return finalRDMlist2,pos,dist

def falsePositiveCheck(inputPair) :
	'''
	Check if input KCF is false positive KCF,
	Input : KCF
	Output : True if KCF is false positive KCF, else return False
	'''
	FP = [('C1b','C1x'),('C1c','C1y'),('C1d','C1z'),('C2b','C2x'),('C2c','C2z'),('C5a','C5x'),('C8x','C2x'),('C8x','C2b'),('C2y','C8y'),('C2c','C8y'),('N1b','N1x'),('N1c','N1y'),('N2b','N2x'),('N1b','N4x'),('N1x','N4x'),('N1c','N4y'),('N1y','N4y'),('N2y','N5y'),('O2a','O2x'),('O5a','O5x'),('O7a','O7x'), ('S2a','S2x'),('S3a','S3x')]
	if inputPair in FP :
		return True
	else :
		return False


def overlappingRDM(inputO) :
	'''
	Check if input KCF is overlapping KCF,
	Input : KCF
	Output : True if KCF is overlapping KCF, else return False
	'''
	OverR_cList = ['C4a','C5a','C5x','C6a','C7a''C7x']
	OverR_OList = ['O4a','O5a','O5x','O6a','O7a','O7x']
	if inputO in OverR_cList :
		return True
	elif inputO in OverR_OList :
		return True
	else :
		return False

def delete_neighbors(posc,temp) :
	'''
	This function is used in find_d_m() function
	'''
	for i in posc :
		del temp[i]

	return temp

def neighborsToDict(N,ItoM,MtoKCF) :
	'''
	This function is used in find_d_m() function to convert neighbors list to dict
	'''
	temp={}
	for i in N :
		if 'H'== i[1] :
			hIndex = max(map(lambda x: x[0],N))+1000
			if hIndex not in temp.keys() :
				temp[hIndex] = 'H'
			else :
				temp[hIndex+1] = 'H'
		else :
			temp[ItoM[i[0]]] = MtoKCF[ItoM[i[0]]]
	return temp

def find_d_m(r_neighbors, p_neighbors, r_index_to_mappedNUmber, p_index_to_mappedNUmber,r_mappedNumber_to_KCF,p_mappedNumber_to_KCF):
	'''
	Map position wise KCF for D and M part and save it
	Input : reactant-neighbors, product-neighbors, dict of reactant-index-to-mapped-position, dict of product-index-to-mapped-position,
			dict of reactant-mapped-position-to-KCF, dict of product-mapped-position-to-KCF
	For each mapped number of reactant D part, product KCF of same mapped number is noted. If same mapped number is missing, then '*' is added.
	Same is repeated with product D part if not checked and '*' is added in reactant.
	Similar approach is taken to find KCF for M part.

	Return : d and m
	d = list of tuple KCF for reactant and product at D part eg [('C1a','C2a')]
	m =  list of tuple KCF for reactant and product at M part


	'''
	reactant_KCF_neighbors = neighborsToDict(r_neighbors,r_index_to_mappedNUmber,r_mappedNumber_to_KCF)
	product_KCF_neighbors = neighborsToDict(p_neighbors,p_index_to_mappedNUmber,p_mappedNumber_to_KCF)
	d,m={},{}
	posCheck =[]
	# if the position matches between reactant and product
	for p in set(reactant_KCF_neighbors.keys()).intersection(set(product_KCF_neighbors.keys())) :
		if reactant_KCF_neighbors[p]==product_KCF_neighbors[p] :
			m[p]=(reactant_KCF_neighbors[p],product_KCF_neighbors[p])
		else :
			d[p]=(reactant_KCF_neighbors[p],product_KCF_neighbors[p])		
		posCheck.append(p)
	# removed matched position KCF from reactant and product 
	if len(posCheck)>0 :
		reactant_KCF_neighbors = delete_neighbors(posCheck,reactant_KCF_neighbors)
		product_KCF_neighbors = delete_neighbors(posCheck,product_KCF_neighbors)		

	# if the neigbor KCF matched but position does not matches
	rposcheck,pposcheck = [],[]
	for i,j in reactant_KCF_neighbors.items() :
		tempq = [k for k, v in product_KCF_neighbors.items() if v == j]
		pflag=0
		for q in tempq :
			if q not in pposcheck :
				pposcheck.append(q)
				pflag=1
				break
		if pflag==1 :
			m[i]=(j,j)
			rposcheck.append(i)

	# removed matched KCF positions but different position in reactant and product
	if len(pposcheck)>0 :
		reactant_KCF_neighbors = delete_neighbors(rposcheck,reactant_KCF_neighbors)
		product_KCF_neighbors = delete_neighbors(pposcheck,product_KCF_neighbors)
		
	# check for the remaining KCF 
	count= 0 # to check track that while loop does not become infinity loop

	while reactant_KCF_neighbors or product_KCF_neighbors :
		# if one KCF is left in reactant and product, combine them
		if (len(reactant_KCF_neighbors)==1) and (len(product_KCF_neighbors)==1) :
			p,q = list(reactant_KCF_neighbors.keys())[0],list(product_KCF_neighbors.keys())[0]
			if reactant_KCF_neighbors.values()==product_KCF_neighbors.values() :
				m[p]=(reactant_KCF_neighbors[p],product_KCF_neighbors[q])
			else :
				d[p]=(reactant_KCF_neighbors[p],product_KCF_neighbors[q])
			del reactant_KCF_neighbors[p]
			del product_KCF_neighbors[q]
			
		#if one kCF is left in reactant
		elif (len(reactant_KCF_neighbors)>=1) and (len(product_KCF_neighbors)==0) :
			p = list(reactant_KCF_neighbors.keys())[0]
			d[p]=(reactant_KCF_neighbors[p],'*')
			del reactant_KCF_neighbors[p]
			
		# if one KCF is left in product
		elif (len(reactant_KCF_neighbors)==0) and (len(product_KCF_neighbors)>=1) :
			q =list(product_KCF_neighbors.keys())[0]
			d[q]=('*',product_KCF_neighbors[q])
			del product_KCF_neighbors[q]
			

		# if more than 1 KCF in reactant and product
		else :
			# check false positive case, if true combine them as D
			rposcheck,pposcheck = [],[]
			found=0
			for i in reactant_KCF_neighbors.values() :
				for j in product_KCF_neighbors.values() :
					if (falsePositiveCheck((i,j))) :
						p = list(reactant_KCF_neighbors.keys())[list(reactant_KCF_neighbors.values()).index(i)]
						q = list(product_KCF_neighbors.keys())[list(product_KCF_neighbors.values()).index(j)]
						d[p]=(i,j)
						rposcheck.append(p)
						pposcheck.append(q)
					elif (falsePositiveCheck((j,i))) :
						p = list(reactant_KCF_neighbors.keys())[list(reactant_KCF_neighbors.values()).index(i)]
						q = list(product_KCF_neighbors.keys())[list(product_KCF_neighbors.values()).index(j)]
						d[p]=(i,j)
						rposcheck.append(p)
						pposcheck.append(q)

					elif (len(i)>1) and (len(j)>1) : 
						if (i[0:2]==j[0:2])  :
							p = list(reactant_KCF_neighbors.keys())[list(reactant_KCF_neighbors.values()).index(i)]
							q = list(product_KCF_neighbors.keys())[list(product_KCF_neighbors.values()).index(j)]
							d[p]=(i,j)
							rposcheck.append(p)
							pposcheck.append(q)
						elif (len(i)>2) and (len(j)>2) : 
							if ((i[0]+i[2])==(j[0]+j[2]))  :
								p = list(reactant_KCF_neighbors.keys())[list(reactant_KCF_neighbors.values()).index(i)]
								q = list(product_KCF_neighbors.keys())[list(product_KCF_neighbors.values()).index(j)]
								d[p]=(i,j)
								rposcheck.append(p)
								pposcheck.append(q)

					elif i[0]==j[0] :
						p = list(reactant_KCF_neighbors.keys())[list(reactant_KCF_neighbors.values()).index(i)]
						q = list(product_KCF_neighbors.keys())[list(product_KCF_neighbors.values()).index(j)]
						d[p]=(i,j)
						rposcheck.append(p)
						pposcheck.append(q)
					if len(rposcheck)>0 :
						break
				if len(rposcheck)>0 :
						break


				
			if len(rposcheck)==0 :
				# nothing matches, randomly assign any KCF of reactant to any KCF of product.
				p = list(reactant_KCF_neighbors.keys())[list(reactant_KCF_neighbors.values()).index(i)]
				q = list(product_KCF_neighbors.keys())[list(product_KCF_neighbors.values()).index(j)]
				d[p]=(reactant_KCF_neighbors[p],product_KCF_neighbors[q])
				reactant_KCF_neighbors = delete_neighbors([p],reactant_KCF_neighbors)
				product_KCF_neighbors = delete_neighbors([q],product_KCF_neighbors)
				
			# removed matched KCF but different position from reactant and product 
			elif len(rposcheck) >0 :
				reactant_KCF_neighbors = delete_neighbors(rposcheck,reactant_KCF_neighbors)
				product_KCF_neighbors = delete_neighbors(pposcheck,product_KCF_neighbors)
				
				continue

		if (len(reactant_KCF_neighbors)==0) and (len(product_KCF_neighbors)==0) :
			break
		count += 1
		if count==10 :
			return {0:'error'},{0:'error'}
	
	if len(d)==0 :
		d[0]=('*','*')
	if len(m)==0 :
		m[0]=('*','*')
	return d,m


def completeRDM(reaction_center, dissimilar, match) :
	'''
	Write RMD as reactant_R-Product_R:reactant_D-product_D:reactant_M-product_M
	Input : KCF of reaction-center, KCF of D part and KCF of M part
	Return : R:D:M
	'''
	R = reaction_center[0]+"-"+reaction_center[1]
	D=''
	t=0
	for i in (dissimilar.keys()) :
		if t==0 :
			D+=dissimilar[i][0]
			t=1
		else :
			D=D+"+"+dissimilar[i][0]
	D+='-'
	t=0
	for i in (dissimilar.keys()) :
		if t==0 :
			D+=dissimilar[i][1]
			t=1
		else :
			D=D+"+"+dissimilar[i][1]

	M=''
	t=0
	
	for i in (match.keys()) :
		if t==0 :
			M+=match[i][0]
			t=1
		else :
			M=M+"+"+match[i][0]
	M+='-'
	t=0
	for i in (match.keys()) :
		if t==0 :
			M+=match[i][1]
			t=1
		else :
			M=M+"+"+match[i][1]

	return R,D,M
