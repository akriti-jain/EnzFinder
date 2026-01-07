import re
import kcfconvoy as kcf
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.rdmolops import GetDistanceMatrix
RDLogger.DisableLog('rdApp.*')

def smiesTocanonical(Input_Smiles):
	can = Chem.MolToSmiles(Chem.MolFromSmiles(Input_Smiles))
	# print ('initial smile', can)
	return kcfv(can)

def kcfv(smiles) :
	k = kcf.KCFvec()
	k.input_smiles(smiles)
	k.convert_kcf_vec()
	KCF_index_mapList = {key: value["kegg_atom"] for key, value in k.kegg_atom_label.items()} #atom_index:atom_kcf
	# print ("KCF index malist", KCF_index_mapList)
	return kcfextension(smiles,KCF_index_mapList),find_neighbors(smiles,KCF_index_mapList)

def kcfextension(SMILES,KCF_index_mapList) :
	mol = Chem.MolFromSmiles(SMILES)
	KCF_mapped_subs ={}
	for i,atom in enumerate(mol.GetAtoms()): # take an atom in molecule
		val=''
		for nbr in atom.GetNeighbors() : # take each neighbor and check
			val+=nbr.GetSymbol()
		KCF_mapped_subs[atom.GetAtomMapNum()]=KCF_index_mapList[i]+"_"+''.join(sorted(val))
	# print ("KCF mapped subs", KCF_mapped_subs)
	return KCF_mapped_subs

def find_mapped_number_in_molecule (Smiles) :
	mol = Chem.MolFromSmiles(Smiles)
	mappedNumber = []
	for i,atom in enumerate(mol.GetAtoms()):
		mappedNumber.append(atom.GetAtomMapNum())
	return mappedNumber

def find_neighbors(Smiles,KCFmapping) :
	mol = Chem.MolFromSmiles(Smiles)
	mol = Chem.AddHs(mol)
	mappedNumber_to_index ={}
	index_to_mappedNUmber = {}
	mappedNumber_to_atom = {}
	index_to_atom = {}
	index_to_neighbors = {}
	mappedNumber_to_KCF = {}
	for i,atom in enumerate(mol.GetAtoms()):
		# print (atom.GetAtomMapNum(),atom.GetIdx())
		mappedNumber_to_index[atom.GetAtomMapNum()]= atom.GetIdx()
		index_to_mappedNUmber[atom.GetIdx()] = atom.GetAtomMapNum()
		mappedNumber_to_atom[atom.GetAtomMapNum()] = atom.GetSymbol()
		index_to_atom[atom.GetIdx()] = atom.GetSymbol()
		index_to_neighbors[atom.GetIdx()]= [(nbr.GetIdx(),nbr.GetSymbol()) for nbr in atom.GetNeighbors()]
		if atom.GetSymbol() != 'H' :
			mappedNumber_to_KCF[atom.GetAtomMapNum()]= KCFmapping[atom.GetIdx()]
			# print (atom.GetIdx(),KCFmapping[atom.GetIdx()],atom.GetAtomMapNum(), atom.GetSymbol(),[(nbr.GetIdx(),nbr.GetSymbol()) for nbr in atom.GetNeighbors()])
	# print ('map to kcf',mappedNumber_to_KCF)
	return ([mappedNumber_to_index,index_to_mappedNUmber,mappedNumber_to_atom,index_to_atom,index_to_neighbors,mappedNumber_to_KCF])



def calculate_farthest_atom_distance_from_reaction_centre(mol_smiles, kcfPos_index) :
	mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(mol_smiles)))	
	return round(max(GetDistanceMatrix(mol)[kcfPos_index]))

def sortRDM(dbRDM) :
	temp_rdmList = list(([sorted(re.split("\+",x)) if "+" in x else x for x in y ]) for y in list(map(lambda x : re.split("-",x),re.split(":",dbRDM))))
	temp_R = sorted(temp_rdmList[0])
	# print (temp_rdmList)
	Rpos = temp_rdmList[0].index(temp_R[0])
	Ppos = temp_rdmList[0].index(temp_R[1])
	rdm = ''
	for i in range(3) :
		if i==0:
			rdm=temp_R[0]+"-"+temp_R[1]
		else :
			if temp_R[0]==temp_R[1] :
				Ppos=Rpos+1
			if str(type(temp_rdmList[i][0])) != "<class 'list'>" :
				rdm+=":"+temp_rdmList[i][Rpos]+"-"+temp_rdmList[i][Ppos]
			else :	
				rdm+=":"+"+".join(temp_rdmList[i][Rpos])+"-"+"+".join(temp_rdmList[i][Ppos])
	return rdm


def identify_difference(reactant_details, product_details,react_fn_kcf,prod_fn_kcf,r_smiles,p_smiles):
	temp_kcf_mapped_pos = {}
	temp_farthest_atom_react_prod = {}
	#print (r_smiles,p_smiles)
	reactant_mappedNumber_to_index,reactant_index_to_mappedNUmber,reactant_mappedNumber_to_atom,reactant_index_to_atom,reactant_index_to_neighbors,reactant_mappedNumber_to_KCF = reactant_details[0],reactant_details[1], reactant_details[2], reactant_details[3], reactant_details[4],reactant_details[5]
	product_mappedNumber_to_index,product_index_to_mappedNUmber,product_mappedNumber_to_atom,product_index_to_atom,product_index_to_neighbors,product_mappedNumber_to_KCF = product_details[0],product_details[1], product_details[2], product_details[3], product_details[4],product_details[5]
	
	finalRDMlist = []
	RCformed =[]
	overkcfList=[]
	reactant_sub_str={}
	product_sub_str = {}
	unsorted_sorted_hash = {}
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
						temp = reaction_center+":"+dissimilar+":"+match#np
						sorted_temp = sortRDM(reaction_center+":"+dissimilar+":"+match)#np
						unsorted_sorted_hash[sorted_temp] = temp#np
						if dissimilar=='*-*' :
							print ("No atom in D", (reaction_center, dissimilar, match))
						
						if (reaction_center+":"+dissimilar+":"+match) not in finalRDMlist:
							sorted_rdm = sortRDM(reaction_center+":"+dissimilar+":"+match)
							#print (sorted_rdm)
							finalRDMlist.append(sorted_rdm)
							temp_kcf_mapped_pos[sorted_rdm]=kcfMappedpos
							reactant_farthest_atom = calculate_farthest_atom_distance_from_reaction_centre(r_smiles,reactant_mappedNumber_to_index[kcfMappedpos])
							product_farthest_atom = calculate_farthest_atom_distance_from_reaction_centre(p_smiles,product_mappedNumber_to_index[kcfMappedpos])
							temp_farthest_atom_react_prod[sorted_rdm]=','.join([str(reactant_farthest_atom),str(product_farthest_atom)])
						#print (unsorted_sorted_hash)
							
	
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

	pos, dist, unsortedRDM = [],[],[]

	#print (finalRDMlist2[0])
	
	for var in finalRDMlist2 :
		#print (var)
		pos.append(temp_kcf_mapped_pos[var])
		dist.append(temp_farthest_atom_react_prod[var])
		unsortedRDM.append(unsorted_sorted_hash[var])
	
	return finalRDMlist2,unsortedRDM,pos,dist

def falsePositiveCheck(inputPair) :
	FP = [('C1b','C1x'),('C1c','C1y'),('C1d','C1z'),('C2b','C2x'),('C2c','C2z'),('C5a','C5x'),('C8x','C2x'),('C8x','C2b'),('C2y','C8y'),('C2c','C8y'),('N1b','N1x'),('N1c','N1y'),('N2b','N2x'),('N1b','N4x'),('N1x','N4x'),('N1c','N4y'),('N1y','N4y'),('N2y','N5y'),('O2a','O2x'),('O5a','O5x'),('O7a','O7x'), ('S2a','S2x'),('S3a','S3x')]
	if inputPair in FP :
		return True
	else :
		return False

def neighborsToDict(N,ItoM,MtoKCF) :
	temp={}
	for i in N :
		if 'H'== i[1] :
			hIndex = max(map(lambda x: x[0],N))+1000
			# print (hIndex,temp)
			if hIndex not in temp.keys() :
				temp[hIndex] = 'H'
			else :
				temp[hIndex+1] = 'H'
		else :
			temp[ItoM[i[0]]] = MtoKCF[ItoM[i[0]]]
	# print ('temp',temp)
	return temp

def delete_neighbors(posc,temp) :
	# print ('del neigh',posc, temp)
	for i in posc :
		del temp[i]

	return temp

def overlappingRDM(inputO) :
	# print ('in',inputO)
	OverR_cList = ['C4a','C5a','C5x','C6a','C7a''C7x']
	OverR_OList = ['O4a','O5a','O5x','O6a','O7a','O7x']
	if inputO in OverR_cList :
		return True
	elif inputO in OverR_OList :
		return True
	else :
		return False

def checkHpair(inputPairH) :
	Hpair = [('*','H'),('H','*')]
	if inputPairH in Hpair :
		return True
	else :
		return False



def find_d_m(r_neighbors, p_neighbors, r_index_to_mappedNUmber, p_index_to_mappedNUmber,r_mappedNumber_to_KCF,p_mappedNumber_to_KCF):
	# print ('here',r_neighbors,p_neighbors)
	reactant_KCF_neighbors = neighborsToDict(r_neighbors,r_index_to_mappedNUmber,r_mappedNumber_to_KCF)
	product_KCF_neighbors = neighborsToDict(p_neighbors,p_index_to_mappedNUmber,p_mappedNumber_to_KCF)
	# print (reactant_KCF_neighbors,product_KCF_neighbors)
	d,m={},{}
	posCheck =[]
	# if the position matches between reactant and product
	for p in set(reactant_KCF_neighbors.keys()).intersection(set(product_KCF_neighbors.keys())) :
		if reactant_KCF_neighbors[p]==product_KCF_neighbors[p] :
			m[p]=(reactant_KCF_neighbors[p],product_KCF_neighbors[p])
		else :
			d[p]=(reactant_KCF_neighbors[p],product_KCF_neighbors[p])		
		posCheck.append(p)
	# print ('before del 1st level', reactant_KCF_neighbont in [8,9,10] :
	# removed matched position KCF from reactant and product 
	if len(posCheck)>0 :
		reactant_KCF_neighbors = delete_neighbors(posCheck,reactant_KCF_neighbors)
		product_KCF_neighbors = delete_neighbors(posCheck,product_KCF_neighbors)		

	# if the neigbor KCF matched but position does not matches
	rposcheck,pposcheck = [],[]
	# print ("reached here", reactant_KCF_neighbors, product_KCF_neighbors)
	for i,j in reactant_KCF_neighbors.items() :
		tempq = [k for k, v in product_KCF_neighbors.items() if v == j]
		# print ("Yes, it is",i,j, tempq,pposcheck)
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
		# print ('levell 2',reactant_KCF_neighbors, product_KCF_neighbors, rposcheck,pposcheck) 
		reactant_KCF_neighbors = delete_neighbors(rposcheck,reactant_KCF_neighbors)
		product_KCF_neighbors = delete_neighbors(pposcheck,product_KCF_neighbors)
		# print ('2nd level', reactant_KCF_neighbors, product_KCF_neighbors)

	# check for the remaining KCF 
	count= 0 # to check track that while loop does not become infinity loop

	while reactant_KCF_neighbors or product_KCF_neighbors :
		# print ('while', reactant_KCF_neighbors, product_KCF_neighbors)
		# if one KCF is left in reactant and product, combine them
		if (len(reactant_KCF_neighbors)==1) and (len(product_KCF_neighbors)==1) :
			# print ('while started',reactant_KCF_neighbors, product_KCF_neighbors)
			p,q = list(reactant_KCF_neighbors.keys())[0],list(product_KCF_neighbors.keys())[0]
			if reactant_KCF_neighbors.values()==product_KCF_neighbors.values() :
				m[p]=(reactant_KCF_neighbors[p],product_KCF_neighbors[q])
			else :
				d[p]=(reactant_KCF_neighbors[p],product_KCF_neighbors[q])
			del reactant_KCF_neighbors[p]
			del product_KCF_neighbors[q]
			# print ('while 1st', reactant_KCF_neighbors, product_KCF_neighbors)

		#if one kCF is left in reactant
		elif (len(reactant_KCF_neighbors)>=1) and (len(product_KCF_neighbors)==0) :
			p = list(reactant_KCF_neighbors.keys())[0]
			d[p]=(reactant_KCF_neighbors[p],'*')
			del reactant_KCF_neighbors[p]
			# print ('while 2nd', reactant_KCF_neighbors, product_KCF_neighbors,d)

		# if one KCF is left in product
		elif (len(reactant_KCF_neighbors)==0) and (len(product_KCF_neighbors)>=1) :
			q =list(product_KCF_neighbors.keys())[0]
			d[q]=('*',product_KCF_neighbors[q])
			del product_KCF_neighbors[q]
			# print ('while 3rd', reactant_KCF_neighbors, product_KCF_neighbors)		 


		# if more than 1 KCF in reactant and product
		else :
			# check false positive case, if true combine them as D
			rposcheck,pposcheck = [],[]
			found=0
			for i in reactant_KCF_neighbors.values() :
				for j in product_KCF_neighbors.values() :
					if (falsePositiveCheck((i,j))) :
						# print ('FP',(i,j), reactant_KCF_neighbors,list(reactant_KCF_neighbors.values()).index(i),list(reactant_KCF_neighbors.keys())[0])
						p = list(reactant_KCF_neighbors.keys())[list(reactant_KCF_neighbors.values()).index(i)]
						q = list(product_KCF_neighbors.keys())[list(product_KCF_neighbors.values()).index(j)]
						d[p]=(i,j)
						rposcheck.append(p)
						pposcheck.append(q)
					elif (falsePositiveCheck((j,i))) :
						# print ('elseFP',(i,j), reactant_KCF_neighbors,list(reactant_KCF_neighbors.values()).index(i),list(reactant_KCF_neighbors.keys())[0])
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
							# print ('here',i,j)
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
				# print ('nothing matches, randonly')
				p = list(reactant_KCF_neighbors.keys())[list(reactant_KCF_neighbors.values()).index(i)]
				q = list(product_KCF_neighbors.keys())[list(product_KCF_neighbors.values()).index(j)]
				d[p]=(reactant_KCF_neighbors[p],product_KCF_neighbors[q])
				# print ('check here')
				reactant_KCF_neighbors = delete_neighbors([p],reactant_KCF_neighbors)
				product_KCF_neighbors = delete_neighbors([q],product_KCF_neighbors)
				
			# removed matched KCF but different position from reactant and product 
			elif len(rposcheck) >0 :
				# print ('issue here')
				reactant_KCF_neighbors = delete_neighbors(rposcheck,reactant_KCF_neighbors)
				product_KCF_neighbors = delete_neighbors(pposcheck,product_KCF_neighbors)
				
				# print ('while else', reactant_KCF_neighbors, product_KCF_neighbors)
				continue
			# print ("end of else statem")

		if (len(reactant_KCF_neighbors)==0) and (len(product_KCF_neighbors)==0) :
			break
		count += 1
		if count==10 :
			# print ("Becoming infinite loop. Check KCF neighbors and exceptions.", reactant_KCF_neighbors,product_KCF_neighbors)
			return {0:'error'},{0:'error'}
	
	# print ('d and m' ,d,m)
	if len(d)==0 :
		d[0]=('*','*')
	if len(m)==0 :
		m[0]=('*','*')
	# print ('check this also' ,d,m)
	return d,m


def completeRDM(reaction_center, dissimilar, match) :
	# print ('joining',reaction_center, dissimilar,match)
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
	# print (match)
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

def generate_RDM(reactant, product):
	#  If only reactant and product molecule are similar, then find reaction centre 
	# Convert the reaction into canonical form 
	# Also, Convert the atom of a molecule into KCF 
	# and find neighbours with index and mapped position
	
	reaction_RDM, reaction_unsortedRDM, reaction_KCF_mapped_position, reaction_farthest_atom, reactant_product_pair = [],[],[],[],[]
	for react_one in reactant.split(".") :	
		for product_one in product.split(".") :
			reactant_mappedNumber = find_mapped_number_in_molecule(react_one)
			product_mappedNumber = find_mapped_number_in_molecule(product_one)

			if set(reactant_mappedNumber).intersection(set(product_mappedNumber)) :
						reactant_FN_KCF, reactant_details = smiesTocanonical(react_one)		
						product_FN_KCF, product_details = smiesTocanonical(product_one)
						#print (react_one,product_one)
						result_RDM_to_be_saved,result_unsorted_RDM_to_be_saved,KCF_mapped_position,farthest_atom = identify_difference(reactant_details, product_details, reactant_FN_KCF,product_FN_KCF,react_one,product_one)
						#print(result_RDM_to_be_saved,KCF_mapped_position,farthest_atom)
						if len(result_RDM_to_be_saved)>0 :
							reaction_RDM.append(result_RDM_to_be_saved)
							reaction_unsortedRDM.append(result_unsorted_RDM_to_be_saved)
							reaction_KCF_mapped_position.append(KCF_mapped_position)
							reaction_farthest_atom.append(farthest_atom)
							reactant_product_pair.append(react_one+">>"+product_one)
							#row = "\t".join([result_unsorted_RDM_to_be_saved,result_RDM_to_be_saved,KCF_mapped_position,farthest_atom,react_one+">>"+product_one])
							#row = "\t".join([result_unsorted_RDM_to_be_saved,result_RDM_to_be_saved,str(KCF_mapped_position)])
							#print (row)
			

	

	return (reaction_unsortedRDM, reaction_RDM, reaction_KCF_mapped_position, reaction_farthest_atom,reactant_product_pair)
				
# if __name__ == "__main__":
#     main()