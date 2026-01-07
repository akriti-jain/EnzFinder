from chython import smiles
def map_reaction(reaction) :
	try :
		r=smiles(reaction)
		r.reset_mapping(keep_reactants_numbering=True) # this is modified line as per author suggestion to keep numbering intact
		return format(r, 'm')
	except  :
		return 'Error'