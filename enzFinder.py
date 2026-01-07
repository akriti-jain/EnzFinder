import sys
import argparse
from map_reaction_graphormer import map_reaction
from generate_query_RDM_v012 import generate_RDM
from round1 import score_ec_number

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

def main(mapped, input_fileName) :
    flag = 0
    filemarker = input_fileName.replace(".csv","")
    output_folder2 = "user_specified"
    with open (output_folder2+filemarker+'_RDM.csv', 'a') as outfile2 :
        print ("KEGG_ID", "\t", "EC_number", "\t", "new_Graphormer_Mapper_MAPPING", "\t", "Reactant-Product-Used", "\t", "RDM", "\t", "sorted_RDM", "\t", "KCF_mapped_pos", "\t", "Farthest_atom_distance", file = outfile2)
        for line in open(input_fileName,"r") :
            line = line.strip().split(",")
            if flag==0 :
                flag=1
                header = line
                continue
                
            reaction_ID = line[0].strip()
            Enz = line[1].strip()
            rxn = line[2].strip().replace("\"","")

            if mapped == 0 :
                rxn = map_reaction(rxn)
            
            if rxn =='Error' :
                print ("\t".join([reaction_ID,'Could not map the reaction.',"\n"]))
            else :
                reactant, product = rxn.strip().split(">>")[0],rxn.strip().split(">>")[1]
                reaction_unsortedRDM, reaction_RDM, reaction_KCF_mapped_position, reaction_farthest_atom,reactant_product_pair = generate_RDM(reactant, product)
                for i in range(len(reaction_unsortedRDM)) :
                    for j in range(len(reaction_unsortedRDM[i])) :
                        row = "\t".join([reaction_ID, Enz, rxn, str(reactant_product_pair[i]), str(reaction_unsortedRDM[i][j]), str(reaction_RDM[i][j]), str(reaction_KCF_mapped_position[i][j]), str(reaction_farthest_atom[i][j])])
                        print(row, file = outfile2)
                
                name = str(reaction_ID)
                
                output_folder = "user_specified"+filemarker+"/"
                with open (output_folder+name+'round1_ECwise_score.csv', 'a') as out :
                    print ("Reaction-Name", "\t", "query-RDM", "\t", "Matched-RDM-part", "\t", "Matched-part", "\t", "RDMwise_score", "\t", "score", "\t", "predicted-EC-number", "\t", "Associated-Reaction-MetaCycID", file = out)
                
                Top10_round1_ECnumber = score_ec_number(reaction_RDM,reactant_product_pair,name,filemarker)
                print (line[0], "\nTop 10 EC sub-subclass", list(Top10_round1_ECnumber.keys()), "\n")
                

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict EC number for a chemical reaction")
    parser.add_argument("--mapped", required=True, type=int, default=1, help='If reaction is mapped - 1, If reaction is unmapped - 0, Default-1')
    parser.add_argument("--i", required=True, type=str, help='Enter input file path with header. File should have be in following format. column A - Reaction ID, column B - Mapped or unmapped SMARTS of a reaction')
    args = parser.parse_args()
    mapped = args.mapped
    input_fileName = args.i
    main(mapped,input_fileName)
    