# EukDetect

# Installation

**Install conda**

If you do not have the conda package manager installed already, you will need to install it. Follow the instructions for that here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

**Download this repository**

Download this repository & download the Fishare repo
```
git clone https://github.com/allind/EukDetect.git
cd EukDetect
```
Use whatever method you prefer to move the taxonomy_db folder to the project folder

**Create conda environment and install EukDetect**
```
conda env create -n eukdetect -f environment.yml
conda activate eukdetect
python setup.py install
```
