from chython import smiles
'''
	GraphormerMapper is used to map the reactions.

	Input : SMARTS
	Output : atom-to-atom mapped reaction

	If GraphormerMapper is unable to map the reaction, then error is returned.
'''
def map_reaction(reaction) :
	try :
		r=smiles(reaction)
		r.reset_mapping(keep_reactants_numbering=True) # this is modified line as per author suggestion to keep numbering intact
		return format(r, 'm')
	except  :
		return 'Error'