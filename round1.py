import re
from rdkit import Chem
import math
			
def score_ec_number(sorted_RDM, reaction, reaction_name, filemarker) :
	# Import unique RDM database and single molecule cofactor list
	try :
		databaseFreqEC = freqData("./data/uniqueRDM/uniqueRDM_db.csv")
		single_molecule_cofactor = singleCofactor("./data/cofactor/single_cofactor.csv")
	except :
		return "Error: Please check the uniqueRDM file and single cofactor file path in round1.py."
	
	try :
		rdm_score_EC = {}
		#2D list, one reactant product pair counted once even if the number of RDMs is more than 1 for that pair
		for var in range(len(sorted_RDM)) :
		 	# If single molecule such as water is present in a reaction as cofactor
			# Then, assign the cofactor weightage as 0.5 for that reactant pair
			rxn = reaction[var]
			reactant_pair = rxn.split(">>")
			water_cofactor_weightage = 0
			for mol in reactant_pair :
				mol = Chem.MolFromSmiles(mol)
				[a.SetAtomMapNum(0) for a in mol.GetAtoms()]
				if Chem.MolToSmiles(mol) in single_molecule_cofactor :
					water_cofactor_weightage = 0.5
					break 
							
			rdm_all = sorted_RDM[var]
			#list of all RDMs for a reactant product pair 
			for rdm in rdm_all :
				rd = rdm.strip().split(":")[0]+":"+rdm.strip().split(":")[1]
				rrpdm = rdm.strip().split(":")[0]+":"+rdm.strip().split(":")[1].split("-")[0]+":"+rdm.strip().split(":")[2]
				prpdm = rdm.strip().split(":")[0]+":"+rdm.strip().split(":")[1].split("-")[1]+":"+rdm.strip().split(":")[2]
				rm = rdm.strip().split(":")[0]+":"+rdm.strip().split(":")[2]
				r = rdm.strip().split(":")[0]
				dm = rdm.strip().split(":")[1]+":"+rdm.strip().split(":")[2]
				d = rdm.strip().split(":")[1]

				if rdm in databaseFreqEC['RDM'].keys() :
					for ifdata in range(len(databaseFreqEC['RDM'][rdm][1])) :
						ftotal = databaseFreqEC['RDM'][rdm][0]#total RDM count all EC all reaction
						freq = int(databaseFreqEC['RDM'][rdm][1][ifdata].strip().split(":")[1])#one by one EC wise (ifdata - for iteration through each EC) frequency e.g., EC-4.2.1:3
						ec = databaseFreqEC['RDM'][rdm][1][ifdata].strip().split(":")[0]#EC in the format EC-4.2.1 
						dbreaction = databaseFreqEC['RDM'][rdm][2][ifdata].strip().split(":")[1]#Again ifdata would iterate through EC, this would fetch the reactions
						cofactor_weightage = databaseFreqEC['RDM'][rdm][3][ifdata]#cofactor RDM for the EC or not
						if water_cofactor_weightage==0.5 :
							cofactor_weightage=0.5
						part_weight = databaseFreqEC['RDM'][rdm][4]
						score = scoreCalculate(freq,ftotal)*cofactor_weightage*part_weight#calculate score as per formula

						# There is a possibility that there are more than one RDM for a reaction
						# Select reaction that has highest score in different level of match mapped to same EC
						if ec not in rdm_score_EC :
							rdm_score_EC[ec]={rdm:[score,dbreaction,rdm,'RDM']}#if EC was not defined, first time EC
						else :
							if rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][rdm]=[score,dbreaction,rdm,'RDM']#if EC was defined, but not for that RDM then add score for that RDM
							else :
								if score > rdm_score_EC[ec][rdm][0] :
									rdm_score_EC[ec][rdm]=[score,dbreaction,rdm,'RDM']#if EC was defined for that RDM but now a better score is found
										
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
							rdm_score_EC[ec]={rdm:[score,dbreaction,rd,'RD']}
							
						else :
							if rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][rdm]=[score,dbreaction,rd,'RD']
								
							else :
								if score > rdm_score_EC[ec][rdm][0] :
									rdm_score_EC[ec][rdm]=[score,dbreaction,rd,'RD']
									
				
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
							rdm_score_EC[ec]={rdm:[score,dbreaction,prpdm,'RpDM']}
						else :
							if rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][rdm]=[score,dbreaction,prpdm,'RpDM']
							else :
								if score > rdm_score_EC[ec][rdm][0] :
									rdm_score_EC[ec][rdm]=[score,dbreaction,prpdm,'RpDM']
							
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
							rdm_score_EC[ec]={rdm:[score,dbreaction,rrpdm,'RpDM']}
						else :
							if rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][rdm]=[score,dbreaction,rrpdm,'RpDM']
							else :
								if score > rdm_score_EC[ec][rdm][0] :
									rdm_score_EC[ec][rdm]=[score,dbreaction,rrpdm,'RpDM']
								
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
							rdm_score_EC[ec]={rdm:[score,dbreaction,dm,'DM']}
						else :
							if rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][rdm]=[score,dbreaction,dm,'DM']
							else :
								if score > rdm_score_EC[ec][rdm][0] :
									rdm_score_EC[ec][rdm]=[score,dbreaction,dm,'DM']
				
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
							rdm_score_EC[ec]={rdm:[score,dbreaction,d,'D']}
						else :
							if rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][rdm]=[score,dbreaction,d,'D']
							else :
								if score > rdm_score_EC[ec][rdm][0] :
									rdm_score_EC[ec][rdm]=[score,dbreaction,d,'D']
										
						
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
							rdm_score_EC[ec]={rdm:[score,dbreaction,rm,'RM']}
						else :
							if rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][rdm]=[score,dbreaction,rm,'RM']
							else :
								if score > rdm_score_EC[ec][rdm][0] :
									rdm_score_EC[ec][rdm]=[score,dbreaction,rm,'RM']
							
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
							rdm_score_EC[ec]={rdm:[score,dbreaction,r,'R']}
						else :
							if rdm not in rdm_score_EC[ec] :
								rdm_score_EC[ec][rdm]=[score,dbreaction,r,'R']
							else :
								if score > rdm_score_EC[ec][rdm][0] :
									rdm_score_EC[ec][rdm]=[score,dbreaction,r,'R']

		round1_ECnumber_score, round1_ECnumber_RDM_reactions = round1_final_score(rdm_score_EC, reaction_name, filemarker)
		round1_ECnumber_score_sorted = dict(sorted(round1_ECnumber_score.items(), key=lambda item: item[1], reverse=True)[:10])
		Top10_round1_ECnumber = identify_round1_top10(round1_ECnumber_score_sorted,round1_ECnumber_RDM_reactions)
		return Top10_round1_ECnumber

	except :
		return "Error: Error in predicting EC number in round1."

def scoreCalculate(f,t) :
	return (f/t)

def freqData(fname) :
	flag,rdm_ec_freq_rxn = 0,{}
	for line in open(fname,"r") :
		if flag== 0 :
			flag=1
			continue
		line = line.strip().split("\t")
		rdm=line[0].strip()#C1a-C1b:.. etc
		matchedPart = line[1].strip()#RDM/RD/RM etc
		total = int(line[2].strip())
		freq = line[5].strip().split("||")
		reaction = line[6].strip().split("||")
		cofactor = list(map(lambda x : float(x),line[3].strip().split("||")))
		part_weightage = float(line[4].strip())

		#---------------------------------------------------------------------------------------------------------------
		#########  EC-RDM freq and reaction ########
		# rdm_ec_freq_rxn[matchedPart]={rdm1:[total, freq, reaction],rdm2:[total,freq, reaction,cofactor,part_weightage]}
		if matchedPart not in rdm_ec_freq_rxn.keys() :
			rdm_ec_freq_rxn[matchedPart]={}
		else :
			if rdm not in rdm_ec_freq_rxn[matchedPart] :
				rdm_ec_freq_rxn[matchedPart][rdm]=[total,freq,reaction,cofactor,part_weightage]
				#processes cofactor_db data contains rdm which denotes the pattern followed by total count, EC wise frequency, EC wise reaction, is it cofactor and corresponding weight
			else :
				rdm_ec_freq_rxn[matchedPart][rdm][0]+=total#returns total count for a pattern with type RDM.RD etc
				rdm_ec_freq_rxn[matchedPart][rdm][1]+=freq#EC wise count
				rdm_ec_freq_rxn[matchedPart][rdm][2]+=reaction#EC wise reaction
				rdm_ec_freq_rxn[matchedPart][rdm][3]+=cofactor#cofactor weight
				rdm_ec_freq_rxn[matchedPart][rdm][4]=part_weightage#pattern weight
				#Each RpDM has 2 entries, for either side! this part of the code only applies to RpDM for addition#processes cofactor_db data contains rdm which denotes the pattern followed by total count, EC wise frequency, EC wise reaction, is it cofactor and corresponding weight
		#---------------------------------------------------------------------------------------------------------------
	return rdm_ec_freq_rxn#rdm - C1b+C1b+H+O1a-*+C1b+C8x+C8y matched part - D rdm_ec_freq_rxn - [3, ['EC-4.2.1:3'], ['EC-4.2.1:RXN-14874,RXN1A0-6305,RXN-17987'], [1.0], 0.5]


def singleCofactor(fname) :
	temp=[]
	for line in open(fname,"r") :
		line = line.strip().split("\t")
		smi=line[1].strip()
		mol = Chem.MolFromSmiles(smi)
		temp.append(Chem.CanonSmiles(smi))
	return temp	

def round1_final_score(EC_RDM_score_all, reaction_name, filemarker) :
	EC_final_Score = {}
	EC_final_Score_RDM_reactions = {}
	output_folder = "user_defined/"+filemarker+"/"
	for ec in EC_RDM_score_all.keys() :
		score = []
		all_dbReactions = []
		all_ecRDM = []
		all_patterns = []
		all_patterns_type = []
		for ECrdm in EC_RDM_score_all[ec].keys() :
			all_ecRDM.append(ECrdm)
			score.append(EC_RDM_score_all[ec][ECrdm][0])
			all_dbReactions.append(EC_RDM_score_all[ec][ECrdm][1])
			all_patterns.append(EC_RDM_score_all[ec][ECrdm][2])
			all_patterns_type.append(EC_RDM_score_all[ec][ECrdm][3])

		fscore=round(sum(score),3)
		EC_final_Score[ec] = fscore
		EC_final_Score_RDM_reactions[ec] = ['||'.join(all_ecRDM),'||'.join(all_dbReactions)]
		query_rdms = ''
		#for i in all_ecRDM :
		query_rdms = '||'.join(all_ecRDM)
		database_rdms = ''
		database_rdms = '||'.join(all_patterns)
		database_rdm_types = ''
		database_rdm_types = '||'.join(all_patterns_type)
		database_scores = ''
		reactions = '||'.join(all_dbReactions)
		database_scores = '||'.join(map(str, score))
		with open (output_folder+reaction_name+'round1_ECwise_score.csv', 'a') as out :
			print (reaction_name, "\t", query_rdms, "\t", database_rdms, "\t", database_rdm_types, "\t", database_scores, "\t", fscore, "\t", ec, "\t", reactions, file = out)
			
	
		 
	
	return EC_final_Score,EC_final_Score_RDM_reactions

def identify_round1_top10(sorted_ECnumber, EC_rdm_rxn_dict) :
	temp_top_EC = {}
	for k,v in  EC_rdm_rxn_dict.items() :
		if k in sorted_ECnumber.keys() :
			temp_top_EC[k]=v
	return temp_top_EC
