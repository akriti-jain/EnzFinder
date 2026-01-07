Requirement :
RDKit
GraphormerMapper
KCFconvoy
re
math
subprocess

###########################################################################################################

Install chython (GraphormerMapper) python package.

Requirement : Only python3.8+
To install chython, use the following command :
>> pip install chython

Github link : https://github.com/chython/chython

###########################################################################################################

Install KCFconvoy python package

To install KCFconvoy, use the following command :
>> pip install kcfconvoy

For more details to install this package, follow 
https://github.com/KCF-Convoy/kcfconvoy/wiki/Installation-method

now, replace ~/python3.10/site-packages/kcfconvoy-0.0.5-py3.10.egg/kcfconvoy/KCFvec.py with EnzFinder/data/extendedKCF/KCFvec.py


###########################################################################################################

To install RDKit, use the following command :
>> pip install rdkit

###########################################################################################################

To run EnzFinder :

>> python enzFinder.py --mapped 1 --i KEGG_test.csv
>> python predict_EC_last_digit.py

For help :
>> python enzFinder.py -h

usage: enzFinder.py [-h] --mapped MAPPED --i I

Predict EC number for a chemical reaction

options:
  -h, --help       show this help message and exit
  --mapped MAPPED  If reaction is mapped - 1, If reaction is unmapped - 0,
                   Default-1
  --i I            Enter input file path with header. File should have be in
                   following format. column A - Reaction ID, column B - Mapped
                   or unmapped SMARTS of a reaction
