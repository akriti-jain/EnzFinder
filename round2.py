'''
## Copyright Notice
EnzFinder code repository is a TCS proprietary resource and should be used for academic purposes only. The contents of this repository should not be used for any commercial purpose without the consent of ALL the authors involved. By downloading and utilizing the scripts, the user consents that any and all Intellectual Property derived from the EnzFinder code repository is fully owned by TCS in the associated jurisdictions. EnzFinder code repository usage without citation will be considered illegal.
'''
from rdkit import Chem
from rdkit import DataStructs
import re
import os
from rdkit.Chem.rdmolops import PatternFingerprint
from itertools import chain
import subprocess
from round1 import sortRDM

def prioritize_EC(queryReactionID, Top10_EC_round1, reaction_SMILES, RDM_reaction_centre_pos, radius_from_reaction_centre,db_RDM_mapped_pos,db_reaction,cofactor_ec) :
	'''
			Final score of top 10 EC level 3 and prioritize EC level 4

			Input : Query-reaction-name, EC-level3-from-round1, query-reaction-SMARTS, RDM-reaction-center-mapped-position, size-of-the-molecule, 
					database-RDM-reaction-center-mapped-position, database-reaction=db_reaction, cofactor-RDM-EC-list
			For every RDM in an EC, check if reactant and product molecule size of query reaction is less than or equal to 3, if yes, ignore that RDM
			Also, RDM should not be cofactor RDM. RDM is matched with RDM given in cofactor-RDM-EC-list (cofactor_ec)

			Openbabel is used to generate molecule fragments.
			Each database reaction mapped to RDM of an EC number is checked.
			Using the SMARTS of database reaction and query reactions, fragments are generated for different radius depending on the size of the query molecule.
			Center point in a molecule is reaction center mapped position.
			For each pair of fragments of same radius, tanimoto score is calculated and average is calculated.
			EC number are ranked based on final score.

			Output : EC level 3 with initial screening (round 1) and final score (radius= radius 2) in ./result/{query-reaction-name}/{query-reaction-name}__predicted_EC_level3.csv
			and 
				ranked EC level 4 with final score (radius = dynamic radius) in ./result/{query-reaction-name}/{query-reaction-name}__predicted_EC_level4.csv
	'''
	try :
		reaction_EC_3digit_predict,EC_3digit_predict = {},{}
		EC_4digit_predict, reaction_EC_4digit_predict = {},{}
		for topEC in Top10_EC_round1.keys() :
			queryrdm_all = Top10_EC_round1[topEC][0].split("||") # Query RDM
			Top1_EC_initial_score = Top10_EC_round1[topEC][3]
			rdm_list, EC_picked = {},0 
			for i in range(len(queryrdm_all)) :
				qrdm_check = queryrdm_all[i]
				# Condition1 : Check the size of query molecule, it should be more than 3.
				molSize = check_molecule_size(reaction_SMILES[qrdm_check])
				if molSize==0 :
					# Check if query RDM is a cofactor or not
					if (sortRDM(qrdm_check) in cofactor_ec.keys()) and (topEC in cofactor_ec[sortRDM(qrdm_check)]) :
						pass
					else :
						rdm_list[qrdm_check] = [Top10_EC_round1[topEC][1].split("||")[i],Top10_EC_round1[topEC][2].split("||")[i],reaction_SMILES[qrdm_check],RDM_reaction_centre_pos[qrdm_check],radius_from_reaction_centre[qrdm_check]]
			
			# This is the list of all selected RDM that passed the above condition
			fragment_tanimoto_score_3digit, fragment_tanimoto_score_4digit, database_reaction_id_3digit,database_reaction_id_4digit=[],[],[],[]
			if len(rdm_list) >= 1 :
				# Once RDM is selected, pick the database reaction one by one and calculate the tanimoto score.
				for original_rdm in rdm_list.keys() : 
					dbreactionID = list(set(rdm_list[original_rdm][0].split(",")))
					levelOfMatch = rdm_list[original_rdm][1]
					query_reactantPair = rdm_list[original_rdm][2]
					query_mapped_position = rdm_list[original_rdm][3]
					dynamic_radius_reactant = int(rdm_list[original_rdm][4].split(",")[0])
					dynamic_radius_product = int(rdm_list[original_rdm][4].split(",")[1])
					for metacycID in dbreactionID :
						for dbrdm in db_RDM_mapped_pos[metacycID] :
							if levelOfMatch=='RpDM' :
								if len(set(identify_rdm_levelOfmatch(sortRDM(original_rdm), levelOfMatch)).intersection(set(identify_rdm_levelOfmatch(sortRDM(dbrdm),levelOfMatch))))>0 :
									db_mapped_position = db_RDM_mapped_pos[metacycID][dbrdm]
									db_reactantPair = db_reaction[metacycID][0][dbrdm]
									predicted_EC = db_reaction[metacycID][1]
									break
							else : 
								if identify_rdm_levelOfmatch(sortRDM(original_rdm), levelOfMatch) == identify_rdm_levelOfmatch(sortRDM(dbrdm),levelOfMatch) :
									db_mapped_position = db_RDM_mapped_pos[metacycID][dbrdm]
									db_reactantPair = db_reaction[metacycID][0][dbrdm]
									predicted_EC = db_reaction[metacycID][1]
									break

						db_molsize = check_molecule_size(db_reactantPair)
						if db_molsize == 0 :
							predicted_EC_upto_3digit = ".".join(predicted_EC.split(".")[:3])
							simScore_rad1, simScore_rad_all = 0,0
							if EC_picked == 0 :
								EC_picked = 1

							# Find database reaction's fragments
							db_reactant_fragments =  findSubStructure(2,db_reactantPair.split(">>")[0],db_mapped_position,'ec3')
							db_product_fragments = findSubStructure(2,db_reactantPair.split(">>")[1],db_mapped_position,'ec3')

							query_reactant_fragments = findSubStructure(2,query_reactantPair.split(">>")[0],query_mapped_position,'ec3')
							query_product_fragments = findSubStructure(2,query_reactantPair.split(">>")[1],query_mapped_position,'ec3')

							simScore_pairwise_avg = calculateTanimotoScore(db_reactant_fragments,db_product_fragments,query_reactant_fragments,query_product_fragments)
								# simScore_pairwise_avg has reactant score, product score and average.
							simScore_rad1 = simScore_pairwise_avg[0]
							simScore_rad_all = simScore_pairwise_avg[1]
							fragment_tanimoto_score_3digit.append(simScore_rad_all)
							database_reaction_id_3digit.append(metacycID)
								
							if len(predicted_EC.split("."))==4 :
								db_reactant_fragments = findSubStructure(dynamic_radius_reactant,db_reactantPair.split(">>")[0],db_mapped_position,'ec4')
								db_product_fragments = findSubStructure(dynamic_radius_product,db_reactantPair.split(">>")[1],db_mapped_position,'ec4')

								query_reactant_fragments = findSubStructure(dynamic_radius_reactant,query_reactantPair.split(">>")[0],query_mapped_position,'ec4')
								query_product_fragments = findSubStructure(dynamic_radius_product,query_reactantPair.split(">>")[1],query_mapped_position,'ec4')

								simScore_pairwise_avg = calculateTanimotoScore(db_reactant_fragments,db_product_fragments,query_reactant_fragments,query_product_fragments)
								# simScore_pairwise_avg has reactant score, product score and average.
								simScore_rad1 = simScore_pairwise_avg[0]
								simScore_rad_all = simScore_pairwise_avg[1]
								
							
								if predicted_EC_upto_3digit not in EC_4digit_predict :
									EC_4digit_predict[predicted_EC_upto_3digit] = {predicted_EC : simScore_rad_all}
									reaction_EC_4digit_predict[predicted_EC_upto_3digit] = {predicted_EC : [simScore_rad_all, metacycID]}
								else :
									if predicted_EC in EC_4digit_predict[predicted_EC_upto_3digit] :
										if simScore_rad_all> EC_4digit_predict[predicted_EC_upto_3digit][predicted_EC] :
											EC_4digit_predict[predicted_EC_upto_3digit][predicted_EC] = simScore_rad_all
											reaction_EC_4digit_predict[predicted_EC_upto_3digit][predicted_EC] =[simScore_rad_all,metacycID]
										elif simScore_rad_all == EC_4digit_predict[predicted_EC_upto_3digit][predicted_EC] :
											reaction_EC_4digit_predict[predicted_EC_upto_3digit][predicted_EC][1] +=","+metacycID
									else :
										EC_4digit_predict[predicted_EC_upto_3digit][predicted_EC]=simScore_rad_all
										reaction_EC_4digit_predict[predicted_EC_upto_3digit][predicted_EC] =[simScore_rad_all,metacycID]

							
			if len(fragment_tanimoto_score_3digit)>0 :
				maxScore =  max(fragment_tanimoto_score_3digit)
				metacycID = [database_reaction_id_3digit[i] for i, val in enumerate(fragment_tanimoto_score_3digit) if val == maxScore]
				EC_3digit_predict[predicted_EC_upto_3digit] = maxScore
				reaction_EC_3digit_predict[predicted_EC_upto_3digit] =[maxScore,",".join(metacycID),Top1_EC_initial_score]
						
			else :
				reaction_EC_3digit_predict[topEC.split("-")[1]] =[0,Top10_EC_round1[topEC][1],Top1_EC_initial_score]
				EC_3digit_predict[topEC.split("-")[1]] = 0
				EC_picked = 1

			if EC_picked == 0 :
				reaction_EC_3digit_predict[topEC.split("-")[1]] =[0,Top10_EC_round1[topEC][1],Top1_EC_initial_score]
				EC_3digit_predict[topEC.split("-")[1]] = 0

		
		# All the results are saved

		if len(reaction_EC_3digit_predict)>0 :
			if not os.path.exists('./result'): 
				os.makedirs('./result')
			if not os.path.exists('./result/'+queryReactionID): 
				os.makedirs('./result/'+queryReactionID)  
			with open('./result/'+queryReactionID+'/'+queryReactionID+"_predicted_EC_level4.csv","w") as  out,open('./result/'+queryReactionID+'/'+queryReactionID+"_predicted_EC_level3.csv","w") as  out1 :
				out.write('\t'.join(['Reaction-Name',"MetaCyc-Reaction-ID","Predicted-EC-number","Final-Score","Rank","\n"]))
				out1.write('\t'.join(['Reaction-Name',"MetaCyc-Reaction-ID","Predicted-EC-number",'Initial-Screening-Weighted-Score',"Final-Score","\n"]))
					
				sorted_simscore_rad1 = sorted(EC_3digit_predict.items(), key=lambda item: item[1], reverse=True)[:10]
				final_EC_4digit_ranked,final_EC_4digit_zero = {},[]
				for val in sorted_simscore_rad1 :
					k = val[0]
					out1.write("\t".join([queryReactionID,reaction_EC_3digit_predict[k][1],k,str(reaction_EC_3digit_predict[k][2]),str(reaction_EC_3digit_predict[k][0]),'\n'])) 
					
					if k in EC_4digit_predict :
						sorted_simscore_radall_EC = sorted(EC_4digit_predict[k].items(), key=lambda item: item[1], reverse=True)[:10]
						for val1 in sorted_simscore_radall_EC :
							k1 = val1[0]
							if val1[1] not in final_EC_4digit_ranked :
								final_EC_4digit_ranked[val1[1]] = [[queryReactionID,reaction_EC_4digit_predict[k][k1][1],k1,str(reaction_EC_4digit_predict[k][k1][0])]]
							else :
								final_EC_4digit_ranked[val1[1]].append([queryReactionID,reaction_EC_4digit_predict[k][k1][1],k1,str(reaction_EC_4digit_predict[k][k1][0])])
					# else :
					# 	final_EC_4digit_zero.append([queryReactionID,reaction_EC_3digit_predict[k][1],k+".-",'0'])
					

				final_EC_4digit_ranked = sorted(final_EC_4digit_ranked.items(), key=lambda item: item[0], reverse=True)
				rank=1

				for k,v in final_EC_4digit_ranked :
					for v1 in v :
						out.write("\t".join(v1+[str(rank),'\n']))
					rank+=1
				for v in final_EC_4digit_zero :
					out.write("\t".join(v+[str(rank),'\n']))

	
	except :
			return "Error"
		
	return "Done"


def check_molecule_size(RPpair) :
	'''
	Size of reactant and product is checked
	Input : SMARTS reaction
	Output : 0 if size of molecule is greater than 3, else 1
	'''
	reactantProdpair = RPpair.strip().split(">>")
	molLessThan3Atom = 0
	for i in reactantProdpair :
		if i=='NA' :
			continue
		mol = Chem.MolFromSmiles(i)
		# the following condition is specifically for NH3
		if (mol.GetNumAtoms()==1) :
			molLessThan3Atom = 1
			continue
		mol = Chem.AddHs(mol)
		atomCount = mol.GetNumAtoms()
		if atomCount <=3 :
			molLessThan3Atom = 1
			continue
	return molLessThan3Atom

def identify_rdm_levelOfmatch(rdm, level) :
	'''
	Divide the RDM into different RDM patterns

	Input : RDM
	Output : All RDM patterns
	'''
	if level=='RDM' :
		return rdm
	elif level=='RpDM' :
		return (rdm.split(":")[0]+":"+rdm.split(":")[1].split("-")[0]+":"+rdm.split(":")[2],rdm.split(":")[0]+":"+rdm.split(":")[1].split("-")[1]+":"+rdm.split(":")[2])
	elif level=='DM' :
		return rdm.split(":")[1]+":"+rdm.split(":")[2]
	elif level=='RM' :
		return rdm.split(":")[0]+":"+rdm.split(":")[2]
	elif level=='R' :
		return rdm.split(":")[0]
	elif level=='D' :
		return rdm.split(":")[1]
	elif level=='RD' :
		return rdm.split(":")[0]+":"+rdm.split(":")[1]
	else :
		return rdm

def fragment_generate(mSMILE,matomIdx,rad) :
	'''
	This function is called by findSubStructure() to find fragments using openbabel
	Return fragment of specified radius 
	'''
	m = Chem.MolFromSmiles(mSMILE)
	for iatom in m.GetAtoms() :
			if iatom.GetAtomMapNum() == matomIdx :	
				atomIdx = iatom.GetIdx()			
				neighbours_dict = {}
				neigh_list = [atomIdx]
				fcount=0
				indx = [atomIdx]
				while fcount!=rad: # take an atom in molecule
					temp =[]
					for i in indx :
						if i not in neighbours_dict:
							val=[]
							for x in m.GetAtomWithIdx(i).GetNeighbors() :
								val.append(x.GetIdx())
								if x.GetIdx() not in neigh_list :
									neigh_list.append(x.GetIdx())
							neighbours_dict[i] = val
							temp+=val
					fcount+=1
					indx=temp
		

	frag_SMILE = "-:"+Chem.MolFragmentToSmiles(m, neigh_list)
	command = 'obabel -:\"'+frag_SMILE+'\" -osmi -xk'

	result = subprocess.run(['obabel', frag_SMILE, '-osmi', '-xk'], capture_output=True, text=True, check=True)
	return result.stdout.strip().replace("\t","")

def findSubStructure(dradius, mSMILE, matomIdx,eclevel) :
	'''
	Fragments of different radius is calculated
	Input : radius, SMILES, mapped position, EC level 
	Neighbour atoms are saved using rdkit as a list which is used as input in openbabel
	radius may vary. Radius for EC level 3 is 2. Therefore, fragments are generated for radius 1 and 2.
	Radius for EC level 4 is dynamic which depends on the size of the molecule. Fragments are generated from radius 1 to dynamic radius.
	EC level is added because for EC level 3 only radius 2 fragments are generated and 
	for EC level 4 fragments from 1 to dynamic radius fragments are generated.
	output:fragments of all radius
	'''
	
	var = []
	if eclevel == 'ec4' :
		for radius in range(1, dradius+1) :
			fr = fragment_generate(mSMILE,matomIdx,radius)
			var.append(fr)
	else :
		fr = fragment_generate(mSMILE,matomIdx,dradius)
		var.append(fr)

	return var


def calculateTanimotoScore(reactant_qfragments, product_qfragments, reactant_dbfragments, product_dbfragments) :
	'''
	Tanimoto score of each pair of fragments is calculated
	Input : query-reactant-fragments-of-all-radius,query-product-fragments-of-all-radius,database-reactant-fragments-of-all-radius,database-product-fragments-of-all-radius
	Output : [average tanimoto score of radius 1, average tanimoto score of all radius]
	'''
	var = []
	temp_reatant, temp_product = [],[]
	for rad in range(len(reactant_qfragments)) :
		# qReactant1,dbReactant1
		fp1 = Chem.PatternFingerprint(Chem.MolFromSmiles(reactant_qfragments[rad]))
		fp2 = Chem.PatternFingerprint(Chem.MolFromSmiles(reactant_dbfragments[rad]))
		temp_reatant.append(DataStructs.TanimotoSimilarity(fp1,fp2))
	scoreR = sum(temp_reatant)/len(temp_reatant)
	
	
	for rad in range(len(product_qfragments)) :
		fp1 = Chem.PatternFingerprint(Chem.MolFromSmiles(product_qfragments[rad]))
		fp2 = Chem.PatternFingerprint(Chem.MolFromSmiles(product_dbfragments[rad]))
		temp_product.append(DataStructs.TanimotoSimilarity(fp1,fp2))		
	scoreP = sum(temp_product)/len(temp_product)
	
	var.append(round(((temp_reatant[0]+temp_product[0])/2),3))
	var.append(round(((scoreR+scoreP)/2),3))

	return var


						
