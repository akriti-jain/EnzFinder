'''
## Copyright Notice
EnzFinder code repository is a TCS proprietary resource and should be used for academic purposes only. The contents of this repository should not be used for any commercial purpose without the consent of ALL the authors involved. By downloading and utilizing the scripts, the user consents that any and all Intellectual Property derived from the EnzFinder code repository is fully owned by TCS in the associated jurisdictions. EnzFinder code repository usage without citation will be considered illegal.
'''

import re
from rdkit import Chem
import math
			
def score_ec_number(all_RDM, reaction,single_molecule_cofactor,databaseFreqEC) :
	'''
	Round 1 : Top 10 EC level 3 are selected

	Input : All-RDM-generated-for-a reaction, query-reaction, byproduct-cofactor-SMILES-list, unique-RDM-database

	If the query reaction has any byproduct from the list, then the corresponsing RDM is considered as cofactor-RDM and the cofactor-weightage 
	is used.
	query-RDM is splitted into RDM pattern to match different RDM pattern with database
	RDM pattern - RDM, RD, R, DM, D, RpDM, RM
	score of each RDM pattern matched is calculated and saved in rdm_score_EC

	rdm_score_EC = {EC:{unsorted-RDM : [score, database-reaction-name, RDM-pattern-match]}}
	If a reaction is multi-substrate where stoichiometry of the molecule is more than 1, then it will generate same RDM for same EC.
	Therefore, if score changes for any reactant pair, then only score is updated.

	Final score of round 1 is calculated by taking sum of all unsorted-RDM score for a specific EC.
	Prioritize EC based on final score.

	Output : Top10_round1_ECnumber
	Top10_round1_ECnumber = Extracted top 10 EC level 3 from rdm_score_EC

	return {EC level 3:RDM1|RDM2||RDM3, db_reaction1,db_reaction2||db_reaction3, db_Reaction2||db_reaction1,db_reaction3,pattern_match_for_RDM1||pattern_match_for_RDM2||pattern_match_for_RDM3,Initial-Screening-Weighted-Score}
	'''
	try :
		rdm_score_EC = {}
		for var in range(len(all_RDM)) :
			for original_rdm in all_RDM[var] :
				rxn = reaction[original_rdm]
				reactant_pair = rxn.split(">>")
				water_cofactor_weightage = 0

				for mol in reactant_pair : 
					mol = Chem.MolFromSmiles(mol)
					[a.SetAtomMapNum(0) for a in mol.GetAtoms()]
					if Chem.MolToSmiles(mol) in single_molecule_cofactor :
						water_cofactor_weightage = 0.5
						break 

				rdm = sortRDM(original_rdm)
				rd = rdm.strip().split(":")[0]+":"+rdm.strip().split(":")[1]
				rrpdm = rdm.strip().split(":")[0]+":"+rdm.strip().split(":")[1].split("-")[0]+":"+rdm.strip().split(":")[2]
				prpdm = rdm.strip().split(":")[0]+":"+rdm.strip().split(":")[1].split("-")[1]+":"+rdm.strip().split(":")[2]
				rm = rdm.strip().split(":")[0]+":"+rdm.strip().split(":")[2]
				r = rdm.strip().split(":")[0]
				dm = rdm.strip().split(":")[1]+":"+rdm.strip().split(":")[2]
				d = rdm.strip().split(":")[1]

				if rdm in databaseFreqEC['RDM'].keys() :
					for ifdata in range(len(databaseFreqEC['RDM'][rdm][1])) :
						ftotal = databaseFreqEC['RDM'][rdm][0]
						freq = int(databaseFreqEC['RDM'][rdm][1][ifdata].strip().split(":")[1])
						ec = databaseFreqEC['RDM'][rdm][1][ifdata].strip().split(":")[0]
						dbreaction = databaseFreqEC['RDM'][rdm][2][ifdata].strip().split(":")[1]
						cofactor_weightage = databaseFreqEC['RDM'][rdm][3][ifdata]
						if water_cofactor_weightage==0.5 :
							cofactor_weightage=0.5
						part_weight = databaseFreqEC['RDM'][rdm][4]
						score = scoreCalculate(freq,ftotal)*cofactor_weightage*part_weight

						# There is a possibility that there are more than one RDM for a reaction
						# Select reaction that has highest score in different level of match mapped to same EC
						if ec not in rdm_score_EC :
							rdm_score_EC[ec]={original_rdm:[score,dbreaction,'RDM']}
						else :
							if original_rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][original_rdm]=[score,dbreaction,'RDM']
							else :
								if score > rdm_score_EC[ec][original_rdm][0] :
									rdm_score_EC[ec][original_rdm]=[score,dbreaction,'RDM']

				if rd in databaseFreqEC['RD'].keys() :
					for ifdata in range(len(databaseFreqEC['RD'][rd][1])) :
						ftotal = databaseFreqEC['RD'][rd][0]
						freq = int(databaseFreqEC['RD'][rd][1][ifdata].strip().split(":")[1])
						ec = databaseFreqEC['RD'][rd][1][ifdata].strip().split(":")[0]
						dbreaction = databaseFreqEC['RD'][rd][2][ifdata].strip().split(":")[1]
						cofactor_weightage = databaseFreqEC['RD'][rd][3][ifdata]
						if water_cofactor_weightage==0.5 :
							cofactor_weightage=0.5
						part_weight = databaseFreqEC['RD'][rd][4]
						score = scoreCalculate(freq,ftotal)*cofactor_weightage*part_weight
						if ec not in rdm_score_EC :
							rdm_score_EC[ec]={original_rdm:[score,dbreaction,'RD']}
						else :
							if original_rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][original_rdm]=[score,dbreaction,'RD']
							else :
								if score > rdm_score_EC[ec][original_rdm][0] :
									rdm_score_EC[ec][original_rdm]=[score,dbreaction,'RD']

				if prpdm in databaseFreqEC['RpDM'].keys() :
					for ifdata in range(len(databaseFreqEC['RpDM'][prpdm][1])) :
						ftotal = databaseFreqEC['RpDM'][prpdm][0]
						freq = int(databaseFreqEC['RpDM'][prpdm][1][ifdata].strip().split(":")[1])
						ec = databaseFreqEC['RpDM'][prpdm][1][ifdata].strip().split(":")[0]
						dbreaction = databaseFreqEC['RpDM'][prpdm][2][ifdata].strip().split(":")[1]
						cofactor_weightage = databaseFreqEC['RpDM'][prpdm][3][ifdata]
						if water_cofactor_weightage==0.5 :
							cofactor_weightage=0.5
						part_weight = databaseFreqEC['RpDM'][prpdm][4]
						score = scoreCalculate(freq,ftotal)*cofactor_weightage*part_weight
						if ec not in rdm_score_EC :
							rdm_score_EC[ec]={original_rdm:[score,dbreaction,'RpDM']}
						else :
							if original_rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][original_rdm]=[score,dbreaction,'RpDM']
							else :
								if score > rdm_score_EC[ec][original_rdm][0] :
									rdm_score_EC[ec][original_rdm]=[score,dbreaction,'RpDM']

				if rrpdm in databaseFreqEC['RpDM'].keys() :
					for ifdata in range(len(databaseFreqEC['RpDM'][rrpdm][1])) :
						ftotal = databaseFreqEC['RpDM'][rrpdm][0]
						freq = int(databaseFreqEC['RpDM'][rrpdm][1][ifdata].strip().split(":")[1])
						ec = databaseFreqEC['RpDM'][rrpdm][1][ifdata].strip().split(":")[0]
						dbreaction = databaseFreqEC['RpDM'][rrpdm][2][ifdata].strip().split(":")[1]
						cofactor_weightage = databaseFreqEC['RpDM'][rrpdm][3][ifdata]
						if water_cofactor_weightage==0.5 :
							cofactor_weightage=0.5
						part_weight = databaseFreqEC['RpDM'][rrpdm][4]
						score = scoreCalculate(freq,ftotal)*cofactor_weightage*part_weight
						if ec not in rdm_score_EC :
							rdm_score_EC[ec]={original_rdm:[score,dbreaction,'RpDM']}
						else :
							if original_rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][original_rdm]=[score,dbreaction,'RpDM']
							else :
								if score > rdm_score_EC[ec][original_rdm][0] :
									rdm_score_EC[ec][original_rdm]=[score,dbreaction,'RpDM']

				if dm in databaseFreqEC['DM'].keys() :
					for ifdata in range(len(databaseFreqEC['DM'][dm][1])) :
						ftotal = databaseFreqEC['DM'][dm][0]
						freq = int(databaseFreqEC['DM'][dm][1][ifdata].strip().split(":")[1])
						ec = databaseFreqEC['DM'][dm][1][ifdata].strip().split(":")[0]
						dbreaction = databaseFreqEC['DM'][dm][2][ifdata].strip().split(":")[1]
						cofactor_weightage = databaseFreqEC['DM'][dm][3][ifdata]
						if water_cofactor_weightage==0.5 :
							cofactor_weightage=0.5
						part_weight = databaseFreqEC['DM'][dm][4]
						score = scoreCalculate(freq,ftotal)*cofactor_weightage*part_weight
						if ec not in rdm_score_EC :
							rdm_score_EC[ec]={original_rdm:[score,dbreaction,'DM']}
						else :
							if original_rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][original_rdm]=[score,dbreaction,'DM']
							else :
								if score > rdm_score_EC[ec][original_rdm][0] :
									rdm_score_EC[ec][original_rdm]=[score,dbreaction,'DM']

				if d in databaseFreqEC['D'].keys() :
					for ifdata in range(len(databaseFreqEC['D'][d][1])) :
						ftotal = databaseFreqEC['D'][d][0]
						freq = int(databaseFreqEC['D'][d][1][ifdata].strip().split(":")[1])
						ec = databaseFreqEC['D'][d][1][ifdata].strip().split(":")[0]
						dbreaction = databaseFreqEC['D'][d][2][ifdata].strip().split(":")[1]
						cofactor_weightage = databaseFreqEC['D'][d][3][ifdata]
						if water_cofactor_weightage==0.5 :
							cofactor_weightage=0.5
						part_weight = databaseFreqEC['D'][d][4]
						score = scoreCalculate(freq,ftotal)*cofactor_weightage*part_weight
						if ec not in rdm_score_EC :
							rdm_score_EC[ec]={original_rdm:[score,dbreaction,'D']}
						else :
							if original_rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][original_rdm]=[score,dbreaction,'D']
							else :
								if score > rdm_score_EC[ec][original_rdm][0] :
									rdm_score_EC[ec][original_rdm]=[score,dbreaction,'D']

				if rm in databaseFreqEC['RM'].keys() :
					for ifdata in range(len(databaseFreqEC['RM'][rm][1])) :
						ftotal = databaseFreqEC['RM'][rm][0]
						freq = int(databaseFreqEC['RM'][rm][1][ifdata].strip().split(":")[1])
						ec = databaseFreqEC['RM'][rm][1][ifdata].strip().split(":")[0]
						dbreaction = databaseFreqEC['RM'][rm][2][ifdata].strip().split(":")[1]
						cofactor_weightage = databaseFreqEC['RM'][rm][3][ifdata]
						if water_cofactor_weightage==0.5 :
							cofactor_weightage=0.5
						part_weight = databaseFreqEC['RM'][rm][4]
						score = scoreCalculate(freq,ftotal)*cofactor_weightage*part_weight
						if ec not in rdm_score_EC :
							rdm_score_EC[ec]={original_rdm:[score,dbreaction,'RM']}
						else :
							if original_rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][original_rdm]=[score,dbreaction,'RM']
							else :
								if score > rdm_score_EC[ec][original_rdm][0] :
									rdm_score_EC[ec][original_rdm]=[score,dbreaction,'RM']

				if r in databaseFreqEC['R'].keys() :
					for ifdata in range(len(databaseFreqEC['R'][r][1])) :
						ftotal = databaseFreqEC['R'][r][0]
						freq = int(databaseFreqEC['R'][r][1][ifdata].strip().split(":")[1])
						ec = databaseFreqEC['R'][r][1][ifdata].strip().split(":")[0]
						dbreaction = databaseFreqEC['R'][r][2][ifdata].strip().split(":")[1]
						cofactor_weightage = databaseFreqEC['R'][r][3][ifdata]
						if water_cofactor_weightage==0.5 :
							cofactor_weightage=0.5
						part_weight = databaseFreqEC['R'][r][4]
						score = scoreCalculate(freq,ftotal)*cofactor_weightage*part_weight
						if ec not in rdm_score_EC :
							rdm_score_EC[ec]={original_rdm:[score,dbreaction,'R']}
						else :
							if original_rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][original_rdm]=[score,dbreaction,'R']
							else :
								if score > rdm_score_EC[ec][original_rdm][0] :
									rdm_score_EC[ec][original_rdm]=[score,dbreaction,'R']

		round1_ECnumber_score, round1_ECnumber_RDM_reactions = round1_final_score(rdm_score_EC)
		round1_ECnumber_score_sorted = dict(sorted(round1_ECnumber_score.items(), key=lambda item: item[1], reverse=True)[:10])

		Top10_round1_ECnumber = identify_round1_top10(round1_ECnumber_score_sorted,round1_ECnumber_RDM_reactions)

		return Top10_round1_ECnumber

	except :
		return "Error"

def scoreCalculate(f,t) :
	'''
	Input : number of reaction per EC, total number of reaction for respective RDM
	Output : Ratio of number-of-reaction-per-EC and total-number-of-reaction-for-respective-RDM
	'''
	return (f/t)
	

def round1_final_score(EC_RDM_score_all) :
	'''
	Final score of each EC level 3 is calculated

	Input : EC_RDM_score_all = rdm_score_EC

	Add score of all RDM for specific EC

	Output : EC_final_Score_RDM_reactions
	'''
	EC_final_Score = {}
	EC_final_Score_RDM_reactions = {}
	for ec in EC_RDM_score_all.keys() :
		score = []
		all_dbReactions = []
		all_ecRDM = []
		level_of_match = []
		for ECrdm in EC_RDM_score_all[ec].keys() :
			all_ecRDM.append(ECrdm)
			score.append(EC_RDM_score_all[ec][ECrdm][0])
			all_dbReactions.append(EC_RDM_score_all[ec][ECrdm][1])
			level_of_match.append(EC_RDM_score_all[ec][ECrdm][2])

		fscore=round(sum(score),3)
		EC_final_Score[ec] = fscore
		EC_final_Score_RDM_reactions[ec] = ['||'.join(all_ecRDM),'||'.join(all_dbReactions),'||'.join(level_of_match), fscore]
	
		 
	
	return EC_final_Score,EC_final_Score_RDM_reactions

def identify_round1_top10(sorted_ECnumber, EC_rdm_rxn_dict) :
	'''
	Based on the final score, select top 10 EC

	Input : sorted EC number with score, EC_rdm_rxn_dict

	output : top 10 EC level 3 result
	'''
	temp_top_EC = {}
	for k,v in  EC_rdm_rxn_dict.items() :
		if k in sorted_ECnumber.keys() :
			temp_top_EC[k]=v
	return temp_top_EC

def sortRDM(dbRDM) :
	'''
	To avoid double count of similar RDM, RDMs are sorted.
	Input : RDM
	Output : sorted RDM
	'''
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
