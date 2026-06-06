'''
## Copyright Notice
EnzFinder code repository is a TCS proprietary resource and should be used for academic purposes only. The contents of this repository should not be used for any commercial purpose without the consent of ALL the authors involved. By downloading and utilizing the scripts, the user consents that any and all Intellectual Property derived from the EnzFinder code repository is fully owned by TCS in the associated jurisdictions. EnzFinder code repository usage without citation will be considered illegal.
'''
import sys
import argparse
from map_reaction_graphormer import map_reaction
from generate_RDM import generate_RDM 
from round1 import score_ec_number,sortRDM
from round2 import prioritize_EC
from rdkit import Chem

def main(mapped, input_fileName) :
	'''
	Database from ./data folder are saved. 
	Database : unique RDM, cofactor, single cofactor, database reaction with reaction center mapped position
	If there is any error is uploading the database, then it will show error and exit.

	map_reaction is called to map the reaction
	generate_RDM is called to generate RDM of a given mapped reaction
	score_ec_number is called to calculate the select top 10 EC level 3 based on Initial-Screening-Weighted-Score 
	prioritize_EC is called to calculate final score based on tanimoto score for EC level 3 and EC level 4. Save the results

	Input : mapped, input_fileName
	mapped - int, if reactions are atom-to-atom not mapped, then 0, otherise it should be 1. Default value is 1
	input_fileName - Path to input file, expected delimiter is tab

	Output : It will display success or error at every step. 
		If there is any error in uploading database, then program will show error and exit.
		If there is any error at any step, then it will show error message at that step and move to next reaction.
		If there is no error, then results will be saved in ./results/{reaction-name} folder

	'''
	try :
		print ("Loading data. ...")
		databaseFreqEC = freqData("./data/uniqueRDM/uniqueRDM_db.csv")
		single_molecule_cofactor = singleCofactor("./data/cofactor/single_cofactor.csv")		
		db_RDM_mapped_pos, db_reaction = database_upload("./data/metacyc_db_input.csv")
		cofactor_ec = list_of_cofactors('./data/cofactor/cofactor_pair_with_EC.csv')
	except :
		print ("Error: There is an error in uploading the data.")
		return 
	
	print ('EnzFinder has started predicting enzyme for given query reaction(s).')
	flag = 0 
	for line in open(input_fileName,"r") :
		line = line.strip().split("\t")
		if flag==0 :
			flag=1
			header = line
			continue

		reaction_ID = line[0].strip()
		rxn = line[1].strip().replace("\"","")

		print (f'started with reaction name : {reaction_ID}')
		
		# Map the reaction using GraphormerMapper
		if mapped == 0 :
			rxn = map_reaction(rxn)

		# If there is any error in mapping, it will show error and will not predict EC.
		if rxn =='Error' :
			print (f'Error : Could not map the {reaction_ID} query reaction.')
		
		else :
			# Break the reaction into reactant and product
			reactant, product = rxn.strip().split(">>")[0],rxn.strip().split(">>")[1]
			
			#  Call the function to generate the RDM
			reaction_RDM, reaction_KCF_mapped_position, reaction_farthest_atom,reactant_product_pair = generate_RDM(reactant, product)
			if len(reaction_RDM) > 0 :
				print (f'{reaction_ID} : RDM generated.')
				# call the function to predict EC at Step-1 analysis
				Top10_round1_ECnumber = score_ec_number(reaction_RDM,reactant_product_pair,single_molecule_cofactor,databaseFreqEC)
				if 'Error' in Top10_round1_ECnumber :
					print (f'Error : There is an error in Step-1 analysis for {reaction_ID} query reaction.')
				else :
					print (f'{reaction_ID} : Round-1 completed.')
					# Call the function to prioritize EC number up to 4th digit in round2 analysis
					ECnumber = prioritize_EC(reaction_ID,Top10_round1_ECnumber,reactant_product_pair,reaction_KCF_mapped_position,reaction_farthest_atom,db_RDM_mapped_pos,db_reaction,cofactor_ec)
					if 'Error' in ECnumber :
						print (f'Error : There is an error in prioritizing EC number up to level-4 analysis for {reaction_ID} query reaction.')
					else :
						print (f'Successful : Prioritizing EC level-3 and level-4 for {reaction_ID} is completed.')
						print (f'Result EC level-3 and level-4 for {reaction_ID} is saved in ./result/{reaction_ID} folder.')
			else :
				print (f'Error : Could not generate RDM for {reaction_ID} query reaction.')	
	
	print ('All results are saved in ./result/{query reaction name} folder')
	return

def list_of_cofactors(input__cofactor_filename) :
	'''
	Input : path of cofactor list 

	cofactor file = cofactor-name, RDM, EC level 3

	output : cofactor RDM with EC number are saved 
	'''
	cofactors = {}
	flag=0
	for line in open(input__cofactor_filename,"r") :
		if flag==0 :
			flag=1
			continue
		line=line.strip().split("\t")
		rdm = sortRDM(line[1].strip())
		cofactors[rdm]=line[2].strip()
		
	return cofactors

def freqData(fname) :
	'''
	Input : path of unique RDM database, default delimiter tab

	unique RDM database = unique RDM pattern, which-RDM-pattern, total-number-of-reaction,cofactor-weightage-per-EC, RDM-pattern-weighatge,
							count-of-reaction-per-EC, database-reaction-name-mapped-to-RDM
	which-RDM-pattern : R/RD/RM/RpDM/DM/D/RDM
	rdm_ec_freq_rxn = {which-RDM-pattern :{RDM : [total-number-of-reaction,count-of-reaction-per-EC,database-reaction-name-mapped-to-RDM,
				cofactor-weightage-per-EC,RDM-pattern-weighatge]}}

	output : rdm_ec_freq_rxn 
	'''

	flag,rdm_ec_freq_rxn = 0,{}
	for line in open(fname,"r") :
		if flag== 0 :
			flag=1
			continue
		line = line.strip().split("\t")
		rdm=line[0].strip()
		matchedPart = line[1].strip()
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
			else :
				rdm_ec_freq_rxn[matchedPart][rdm][0]+=total
				rdm_ec_freq_rxn[matchedPart][rdm][1]+=freq
				rdm_ec_freq_rxn[matchedPart][rdm][2]+=reaction
				rdm_ec_freq_rxn[matchedPart][rdm][3]+=cofactor
				rdm_ec_freq_rxn[matchedPart][rdm][4]=part_weightage
		#---------------------------------------------------------------------------------------------------------------
	return rdm_ec_freq_rxn

def database_upload(input_filename) :
	'''
	Input : path of reaction database, default delimiter tab

	temp_map_pos = {reaction-name : {RDM : reaction-center-mapped-position}}
	temp_reaction = {reaction-name : [{RDM : reactant-pair-SMARTS}, EC]}

	output : temp_map_pos,  temp_reaction
	'''
	
	temp_map_pos, temp_reaction = {},{}
	flag = 0
	rdmForSmallMol = []
	for line in open(input_filename,"r") :
		if flag == 0 :
			flag = 1
			continue
		line = line.strip().split("\t")
		if line[2].strip() == 'NA' :
			continue
		rid = line[0].strip()
		ec = line[1].strip()
		rdm = line[3].strip()
		srdm = line[4].strip()
		mappedPosition = int(line[5])
		reactantPair = line[2].strip().split(">>")
		if rdm.split(":")[0] == srdm.split(":")[0] :
			reactantPair = reactantPair[0]+">>"+reactantPair[1]
		else :
			reactantPair = reactantPair[1]+">>"+reactantPair[0]
		if rid not in temp_map_pos :
			temp_map_pos[rid] = {rdm:mappedPosition}
		else :
			temp_map_pos[rid][rdm] = mappedPosition

		if rid not in temp_reaction :
			temp_reaction[rid] = [{rdm:reactantPair},ec]
		else :
			temp_reaction[rid][0][rdm] = reactantPair
				
	return temp_map_pos, temp_reaction
		
def singleCofactor(fname) :
	'''
	Input : path of byproduct cofactor list 

	cofactor file = cofactor-name, SMILES
	convert SMILES into canonical SMILES

	output : list of canonical SMILES
	'''
	temp=[]
	for line in open(fname,"r") :
		line = line.strip().split("\t")
		smi=line[1].strip()
		mol = Chem.MolFromSmiles(smi)
		temp.append(Chem.CanonSmiles(smi))
		
	return temp

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Predict EC number for a chemical reaction")
	parser.add_argument("--mapped", required=True, type=int, default=1, help='If reaction is atom-atom mapped - 1, If reaction is not atom-atom mapped - 0, Default-1')
	parser.add_argument("--i", required=True, type=str, help='Enter input file path with header. File should have be in following format. column A - Reaction ID, column B - Atom-atom mapped or unmapped SMARTS of a reaction')
	# parser.add_argument("--o", required=True, type=str, help='Path to Output file')
	args = parser.parse_args()

	mapped = args.mapped
	input_fileName = args.i
	# output_fileName = args.o

	# main(mapped,input_fileName,output_fileName)
	main(mapped,input_fileName)
