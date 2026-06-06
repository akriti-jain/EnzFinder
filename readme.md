
# [*EnzFinder: a sustainable alternative to chemical synthesis*](https://www.biorxiv.org/content/10.64898/2026.02.12.705490v1)
## Authors : Akriti Jain, Nishtha Pandey, and Arijit Roy

## Usage disclaimer
This is a EnzFinder prediction and inference source code necessary to reproduce the results shown in the manuscript. Wherever possible, appropriate sample input and output files are provided for user reference. Any changes made to the source code (except paths to stand-alone programs) are done at your own risk. The authors will not be liable to any discrepancies observed in the results due to changes made to the source code.

## 1. Requirement
python 3.10+
RDKit
chython
obabel


## 2. Data

Following files are required for enzFinder.py
1. data/metacyc_db_input.csv
2. data/cofactor/cofactor_pair_with_EC.csv

round1.py is used to prioritize and select EC level 3. Following files are required for round1.py.
1. data/uniqueRDM/uniqueRDM_db.csv
2. data/cofactor/single_cofactor.csv

Test input file
1. test/sample_input.csv

*kcfconvoy folder should be in same folder as enzFinder.py*


## 3. Setup

```
pip install chython
pip install chytorch-rxnmap
pip install rdkit
sudo apt-get install openbabel
```

## 4. How to use

Run the command for help
```
python enzFinder.py -h
```
Run the following command to predict EC level-3 and EC level-4 for query reaction.
If reaction is atom-atom mapped, then use `--mapped 1`, otherwise `--mapped 0`

```
python enzFinder.py --mapped 1 --i sample_input.csv
```

### Code usage
Detailed instructions on how to use the codes are provided. For any queries related to code usage, contact the corresponding author for more information.

### Copyright Notice
EnzFinder code repository is a TCS proprietary resource and should be used for academic purposes only. The contents of this repository should not be used for any commercial purpose without the consent of ALL the authors involved. By downloading and utilizing the scripts, the user consents that any and all Intellectual Property derived from the EnzFinder code repository is fully owned by TCS in the associated jurisdictions. EnzFinder code repository usage without citation will be considered illegal.

### Contact Us
For further queries related to code usage, please write to us: akriti.j@tcs.com

### Citation
Please cite this article if you use the codes in this repository for your research:

### License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International


                   
