HOW TO USE

1. Create a Python Virtual Environment using venv

python3 -m venv <name>

And activate using:

source <name>/bin/activate

-------------------------------------------------

2. Run the following commands:

pip install -r requirements.txt
pysmt-install --msat

-------------------------------------------------

3. Convert any directory into projected AllSAT by using the script:

python3 parse_to_projected.py <directory>

It will generate a new directory ending in _DIMACS, containing all files trasformed in CNF and projected to the original variables