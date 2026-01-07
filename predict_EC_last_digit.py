# pattern fingerprint, took all reactions mapped to an EC
from rdkit import Chem
from rdkit import DataStructs
import re
import os
from rdkit.Chem.rdmolops import PatternFingerprint
from itertools import chain
import subprocess


def findSubStructure(dradius, mSMILE, matomIdx) :
	m = Chem.MolFromSmiles(mSMILE)
	var = []
	for radius in range(1, dradius+1) :
		for iatom in m.GetAtoms() :
			if iatom.GetAtomMapNum() == matomIdx :	
				atomIdx = iatom.GetIdx()			
				neighbours_dict = {}
				neigh_list = [atomIdx]
				fcount=0
				indx = [atomIdx]
				while fcount!=radius: # take an atom in molecule
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
		
		
		frag_SMILE = Chem.MolFragmentToSmiles(m, neigh_list)
		command = 'obabel -:\"'+frag_SMILE+'\" -osmi -xk'
		
		var.append(subprocess.check_output(command, shell=True, text=True).split("\t")[0].strip())

	return var

def sortRDM(dbRDM) :
	temp_rdmList = list(([sorted(re.split("\+",x)) if "+" in x else x for x in y ]) for y in list(map(lambda x : re.split("-",x),re.split(":",dbRDM))))
	temp_R = sorted(temp_rdmList[0])
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

def database_fragment(input_filename) :
	# The details of database will only be saved if EC is given up to 4th digit.
	temp = {}
	flag = 0
	rdmForSmallMol = []
	for line in open(input_filename,"r") :
		if flag == 0 :
			flag = 1
			continue
		line = line.strip().split("\t")
		# if molecule size is 3 or less, then that RDM deatils are not saved
		reactantProdpair = line[3].strip().split(">>")
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
				# [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
				# rdmForSmallMol.append(Chem.MolToSmiles(mol))
				continue
		# temp[reaction ID] = [rdm, sorted_Rdm, level of match,kcf mapped position,EC number, reactant pair]
		if (molLessThan3Atom==0) and (line[2].strip()!='NA') :
			rid = line[0].strip()
			ec = line[1].strip()
			rdm = line[4].strip()
			srdm = line[5].strip()
			mappedPosition = int(line[6])
			reactantPair = line[3].strip()
			if rid not in temp :
				temp[rid] = [[rdm], [srdm],[], [mappedPosition],ec,[reactantPair]]
				var=[srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1]]
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1].split("-")[0]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1].split("-")[1]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0])
				var.append(srdm.strip().split(":")[1]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[1])
				var.append(srdm.strip())
				temp[rid][2].append(var)
			else :
				temp[rid][0].append(rdm)
				temp[rid][1].append(srdm)
				temp[rid][3].append(mappedPosition)
				var=[srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1]]
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1].split("-")[0]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1].split("-")[1]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0])
				var.append(srdm.strip().split(":")[1]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[1])
				var.append(srdm.strip())
				temp[rid][2].append(var)
				temp[rid][5].append(reactantPair)
	return temp

def query_reaction_fragments(input_filename) :
	temp = {}
	flag = 0
	for line in open(input_filename,"r") :
		if flag == 0 :
			flag = 1
			continue
		line = line.strip().split("\t")
		rid = line[0].strip()
		reactantProdpair = line[3].strip().split(">>")
		molLessThan3Atom = 0
		for i in reactantProdpair :
			if i=='NA' :
				continue
			mol = Chem.MolFromSmiles(i)
			if (mol.GetNumAtoms()==1) :
				molLessThan3Atom = 1
				continue
			mol = Chem.AddHs(mol)
			atomCount = mol.GetNumAtoms()
			if atomCount <=3 :
				molLessThan3Atom=1
				continue
		if molLessThan3Atom == 0 :
			rid = line[0].strip()
			rdm = line[4].strip()
			srdm = line[5].strip()
			mappedPosition = int(line[6])
			#radius = [int(line[7].strip()),int(line[8].strip())]
			radius = [int(line[7].split(",")[0]),int(line[7].split(",")[1])]
			reactantPair = line[3].strip()
			# temp[reaction ID] = [rdm, sorted_Rdm, level of match, kcf mapped position, distance,reactantPair]
			if rid not in temp :
				temp[rid] = [[rdm], [srdm],[], [mappedPosition],[radius],[reactantPair]]
				var=[srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1]]
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1].split("-")[0]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1].split("-")[1]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0])
				var.append(srdm.strip().split(":")[1]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[1])
				var.append(srdm.strip())
				temp[rid][2].append(var)

			else :
				temp[rid][0].append(rdm)
				temp[rid][1].append(srdm)
				temp[rid][3].append(mappedPosition)
				temp[rid][4].append(radius)
				var=[srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1]]
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1].split("-")[0]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[1].split("-")[1]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[0])
				var.append(srdm.strip().split(":")[1]+":"+srdm.strip().split(":")[2])
				var.append(srdm.strip().split(":")[1])
				var.append(srdm.strip())
				temp[rid][2].append(var)
				temp[rid][5].append(reactantPair)
			
	return temp

def calculateTanimotoScore(reactant_qfragments, product_qfragments, reactant_dbfragments, product_dbfragments) :
	var = []
	temp_reatant, temp_product = [],[]
	#change range to 2 for re-ranking of EC numbers at sub-subclass level/Level 3 after round 1 screening
	for rad in range(len(reactant_qfragments)) :
		fp1 = Chem.PatternFingerprint(Chem.MolFromSmiles(reactant_qfragments[rad]))
		fp2 = Chem.PatternFingerprint(Chem.MolFromSmiles(reactant_dbfragments[rad]))
		temp_reatant.append(DataStructs.TanimotoSimilarity(fp1,fp2))
	scoreR = sum(temp_reatant)/len(temp_reatant)
	
	#change range to 2 for re-ranking of EC numbers at sub-subclass level/Level 3 after round 1 screening
	for rad in range(len(product_qfragments)) :
		fp1 = Chem.PatternFingerprint(Chem.MolFromSmiles(product_qfragments[rad]))
		fp2 = Chem.PatternFingerprint(Chem.MolFromSmiles(product_dbfragments[rad]))
		temp_product.append(DataStructs.TanimotoSimilarity(fp1,fp2))		
	scoreP = sum(temp_product)/len(temp_product)
	
	var.append(round(((temp_reatant[0]+temp_product[0])/2),3))
	var.append(round(((scoreR+scoreP)/2),3))

	return var

def list_of_cofactors(input_filename) :
	cofactors = {}
	f=0
	for line in open(input_filename,"r") :
		if f==0 :
			f=1
			continue
		line=line.strip().split("\t")
		rdm = sortRDM(line[1].strip())
		cofactors[rdm]=line[3].strip()
		
	return cofactors


def main( ) :
	db_frag = database_fragment("/home/nishtha/Documents/TCS_Research/python_codes/EnzFinder_part1/EnzFinder/data/metacycv04_RDM_result_21KCF_v06.csv")
	query_frag = query_reaction_fragments("*_RDM.csv file genearted as output with enzfinder.py")
	cofactor_ec = list_of_cofactors("./data/cofactor/cofactor_EC_v02.csv")
	path_to_ec_score = "Path to the *round1_ECwise_score.csv files genearted as output with enzfinder.py"
	path_to_output_file = "User_defined"
	path_to_output_file_3digit = "User_defined"
			
	count = 0
	for filename in os.listdir(path_to_ec_score) :
			score_ec, flag, ec_details ={},0,{}
			print (filename)
			count=0
			for line in open(path_to_ec_score+"/"+filename,"r") :
				if flag==0 :
					flag = 1
					continue
				line = line.strip().split("\t")
				ec = '.'.join(line[6].strip().split(".")[:3]) # predicted EC number
				score = float(line[5].strip()) # calculated score in round1
				rxnid = line[0].strip()
				# if rxnid not in high_confidence_rxn :
				# 	break 
				if ec not in score_ec :
					score_ec[ec] = score 
					ec_details[ec] = line

				else :
					if score > score_ec[ec] :
						score_ec[ec] = score
						ec_details[ec] = line 

			
			if len(score_ec) == 0:
				print ("There is no prediction at 3rd level in round1,", rxnid)
				continue

			# To pick the top EC from 3rd level prediction, sort the dict based on score	
			score_ec = list(score_ec.items())
			score_ec_sorted = sorted(score_ec, key=lambda x: x[1], reverse=True)

			EC_3digit_predict, EC_4digit_predict, reaction_EC_3digit_predict,reaction_EC_4digit_predict={},{}, {},{}
			if len(score_ec_sorted)>10 :
				len_top5i = 10
			else :
				len_top5i = len(score_ec_sorted)
			
			for top5i in range(len_top5i) :
					# pick the top 10 EC number based on score from 3rd digit
				topEC = score_ec_sorted[top5i][0] # score_ec_sorted[EC-number]=score, dict is passed as list, [(EC-number1,score),()
					
				queryReactionID = ec_details[topEC][0].strip()
				if queryReactionID not in query_frag.keys() :
				# 	out.write("\t".join([queryReactionID,ec_details[topEC][7],topEC,'0','\n']))
					continue
				

				queryrdm_all = ec_details[topEC][1].strip().split("||") # Query RDM
				rdm_list = {} # rdm_list[rdm] = [part_matched, matched_rdm,metacyc_reactions]
				rdm_selected_round1 = []
				metacyc_rdm_selected_round1 =[]
				
				
				for i in range(len(queryrdm_all)) :
						
					qrdm_check = queryrdm_all[i]
					query_all_level_of_match_RDM = list(chain.from_iterable(query_frag[queryReactionID][2]))
					if ec_details[topEC][2].strip().split("||")[i] not in query_all_level_of_match_RDM : # check if molecule is less than equal to 3
						pass
					elif (qrdm_check in cofactor_ec.keys()) and (topEC in cofactor_ec[qrdm_check]) : # check if RDM is cofactor
						pass
					else :
						metacyc_reaction_list = ec_details[topEC][7].strip().split("||")[i].split(",")
						if (len(set(metacyc_reaction_list).intersection(set(db_frag.keys()))) == 0) : # check if 4th digit for that EC is given in our database
							pass
						else :
							rdm_list[ec_details[topEC][1].strip().split("||")[i]]=[ec_details[topEC][3].strip().split("||")[i],ec_details[topEC][2].strip().split("||")[i],list(set(metacyc_reaction_list).intersection(set(db_frag.keys())))]
								
				# This is the list of all selected RDM from round1 analysis
				queryrdm_selected_count = len(rdm_list)
				if queryrdm_selected_count >= 1 :
					# for the selected database reaction, use RDM to map the molecule and calculate tanimoto score
					# Select RDM one by one for an EC number and pick db reaction ID to compare structure based on radius
					for rdm_name in rdm_list.keys() : 
						# Select RDM, check if selected database reaction is mapped to this rdm, if yes, calculate tanimoto score
						if len(set(rdm_list[rdm_name][2])) >0 :
							rdm_dbreactionID = list(set(rdm_list[rdm_name][2]))
							rdm = rdm_name
							matched_rdm = rdm_list[rdm_name][1]
							part_matched = rdm_list[rdm_name][0]
							for i in range(len(query_frag[queryReactionID][2])) :
								if matched_rdm in query_frag[queryReactionID][2][i] :
									dynamic_radius_reactant = query_frag[queryReactionID][4][i][0]
									dynamic_radius_product = query_frag[queryReactionID][4][i][1]
									query_unsorted_original_RDM = query_frag[queryReactionID][0][i]
									query_reactantPair = query_frag[queryReactionID][5][i].split(">>")  # SMILES
									query_mapped_position = query_frag[queryReactionID][3][i]
									break
							# Selected one db reaction at a time corresponding to RDM
							for metacycID in rdm_dbreactionID :
								db_level_of_match_RDM = list(chain.from_iterable(db_frag[metacycID][2]))
								if matched_rdm in db_level_of_match_RDM :
									predicted_EC = db_frag[metacycID][4]
									predicted_EC_upto_3digit = ".".join(predicted_EC.split(".")[:3])
										# Identify the pair of substratre and product from query and db
										# Use this information to calculate the tanimoto score for substrate and product
											
									for i in range(len(db_frag[metacycID][2])) :
										if matched_rdm in db_frag[metacycID][2][i] :
											db_unsorted_original_RDM = db_frag[metacycID][0][i]
											db_reactantPair = db_frag[metacycID][5][i].split(">>")
											db_mapped_position = db_frag[metacycID][3][i]
											break
									
									query_unsorted_original_RDM_reaction_centre = query_unsorted_original_RDM.split(":")[0].split("-")
									db_unsorted_original_RDM_reaction_centre = db_unsorted_original_RDM.split(":")[0].split("-")

									sorted_RDM_reaction_centre = rdm.split(":")[0].split("-") # because query and db reaction center are same
									query_molecule = 0 # the reactant reaction center of RDM is taken from reactant
									db_molecule = 0
													
									if query_unsorted_original_RDM_reaction_centre[0]==sorted_RDM_reaction_centre[1] :
										query_molecule = 1 # the reactant reaction center of RDM is taken from product

									if db_unsorted_original_RDM_reaction_centre[0]==sorted_RDM_reaction_centre[1] :
										db_molecule = 1 # the reactant reaction center of RDM is taken from product


									if db_molecule == 0 :
										db_reactant_fragments = findSubStructure(dynamic_radius_reactant,db_reactantPair[0].strip(),db_mapped_position)
										db_product_fragments = findSubStructure(dynamic_radius_product,db_reactantPair[1].strip(),db_mapped_position)
										
									else :
										db_reactant_fragments = findSubStructure(dynamic_radius_reactant,db_reactantPair[1].strip(),db_mapped_position)
										db_product_fragments = findSubStructure(dynamic_radius_product,db_reactantPair[0].strip(),db_mapped_position)
												
									

									if query_molecule == 0 :
										query_reactant_fragments = findSubStructure(dynamic_radius_reactant,query_reactantPair[0].strip(),query_mapped_position)
										query_product_fragments = findSubStructure(dynamic_radius_product,query_reactantPair[1].strip(),query_mapped_position)
									else :
										query_reactant_fragments = findSubStructure(dynamic_radius_reactant,query_reactantPair[1].strip(),query_mapped_position)
										query_product_fragments = findSubStructure(dynamic_radius_product,query_reactantPair[0].strip(),query_mapped_position)
									
									# print ('reactant fragments :',db_reactant_fragments,query_reactant_fragments)
									# print ('product fragments :',db_product_fragments,query_product_fragments )
									
									simScore_pairwise_avg = calculateTanimotoScore(db_reactant_fragments,db_product_fragments,query_reactant_fragments,query_product_fragments)
										# simScore_pairwise_avg has reactant score, product score and average.
									simScore_rad1 = simScore_pairwise_avg[0]
									simScore_rad_all = simScore_pairwise_avg[1]
									# print (simScore_rad1,simScore_rad_all, predicted_EC)
									
									if predicted_EC_upto_3digit not in reaction_EC_3digit_predict :
										EC_3digit_predict[predicted_EC_upto_3digit] = simScore_rad1
										reaction_EC_3digit_predict[predicted_EC_upto_3digit] =[simScore_rad1,metacycID]
									else :
										# print (simScore_rad1,EC_3digit_predict[predicted_EC_upto_3digit])
										if simScore_rad1> EC_3digit_predict[predicted_EC_upto_3digit] :
											EC_3digit_predict[predicted_EC_upto_3digit] = simScore_rad1
											reaction_EC_3digit_predict[predicted_EC_upto_3digit] =[simScore_rad1,metacycID]
										elif simScore_rad1 == EC_3digit_predict[predicted_EC_upto_3digit] :
											reaction_EC_3digit_predict[predicted_EC_upto_3digit][1] +=","+metacycID
									
									if len(predicted_EC.split("."))==4 :
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
										
			if len(reaction_EC_3digit_predict)>0 :
				with open(path_to_output_file+queryReactionID+"_round2_tanimoto_score_dynamic_radius.csv","w") as  out,open(path_to_output_file_3digit+queryReactionID+"_round2_tanimoto_score_radius1_3digit.csv","w") as  out1 :
					out.write('\t'.join(['Reaction-Name',"MetaCyc-Reaction-ID","Predicted-EC-number","Tanimoto_score_average_radius_dynamic_radius","\n"]))
					out1.write('\t'.join(['Reaction-Name',"MetaCyc-Reaction-ID","Predicted-EC-number","Tanimoto_score_radius1","\n"]))
				
					sorted_simscore_rad1 = sorted(EC_3digit_predict.items(), key=lambda item: item[1], reverse=True)[:10]

					for val in sorted_simscore_rad1 :
						k = val[0]
						out1.write("\t".join([queryReactionID,reaction_EC_3digit_predict[k][1],k,str(reaction_EC_3digit_predict[k][0]),'\n'])) 
						
						if k in EC_4digit_predict :
							sorted_simscore_radall_EC = sorted(EC_4digit_predict[k].items(), key=lambda item: item[1], reverse=True)[:10]

							for val1 in sorted_simscore_radall_EC :
								k1 = val1[0]
								out.write("\t".join([queryReactionID,reaction_EC_4digit_predict[k][k1][1],k1,str(reaction_EC_4digit_predict[k][k1][0]),'\n']))

	return
						
					
			
if __name__ == "__main__":
    main()
