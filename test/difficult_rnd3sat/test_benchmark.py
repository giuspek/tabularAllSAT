import time
import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description='Run benchmark tests and return CPU time')

parser.add_argument('--directory', type=str, default='rnd3sat', help='Directory of the benchmark files')
parser.add_argument("--output", type=str, default='result.csv', help='Name outputfile')

args = parser.parse_args()

directory = args.directory + "/" if not args.directory.endswith("/") else args.directory
entries = os.listdir(directory)
entries.sort()
outputfilename = directory + args.output
outputfile = open(outputfilename, "w")
outputfile.write("filename,n_solutions,BC,NBC,BDD,TABULARALLSAT,TABULARALLSAT_models\n")
outputfile.close()

timeout_bc = True
timeout_nbc = False
timeout_bdd = False
timeout_tabularallsat = False

timeout_time = 1200
counter = 0
print(entries)
for file in entries[:]:
	sol = -1
	print("FILE: {}".format(file))
	if file[-3:] != "cnf":
		continue
	
	if(timeout_nbc):
		start_nbc = 0
		end_nbc = timeout_time
	else:
	#output = check_output("./mc_tabularasat_nbc/cdcl-vsids/solver binary_clauses/{} -c -q".format(file), stderr=STDOUT, timeout=1)
		try:
			cmd = ["./../nbc_minisat_all-1.0.2/nbc_minisat_all_release", "{}{}".format(directory, file)]
			start_nbc = 0
			result = subprocess.run(cmd, capture_output=True, timeout=timeout_time)
			result = result.stdout.decode().splitlines()
			model_count = result[-1].split(" ")[-1]
			end_nbc = float(result[-6].split(" ")[-2])
			if sol == -1:
				sol = int(model_count)
			elif not sol == int(model_count):
				end_nbc = -1
				start_nbc = 0
		except subprocess.TimeoutExpired:
			start_nbc = 0
			end_nbc = timeout_time
			#timeout_nbc = True

	
	if(timeout_bc):
		start_bc = 0
		end_bc = timeout_time
	else:
	#output = check_output("./mc_tabularasat_bc/cdcl-vsids/solver binary_clauses/{} -c -q".format(file), stderr=STDOUT, timeout=1)
		try:
			cmd = ["./../bc_minisat_all-1.1.2/bc_minisat_all_release", "{}{}".format(directory, file)]
			start_bc = 0
			result = subprocess.run(cmd, capture_output=True, timeout=timeout_time)
			result = result.stdout.decode().splitlines()
			model_count = result[-1].split(" ")[-1]
			end_bc = float(result[-7].split(" ")[-2])
			if sol == -1:
				sol = int(model_count)
			elif not sol == int(model_count):
				end_bc = -1
				start_bc = 0
		except subprocess.TimeoutExpired:
			start_bc = 0
			end_bc = timeout_time
			#timeout_bc = True
	
		

	if(timeout_bdd):
		start_bdd = 0
		end_bdd = timeout_time
	else:
	#output = check_output("./mc_tabularasat_bdd/cdcl-vsids/solver binary_clauses/{} -c -q".format(file), stderr=STDOUT, timeout=1)
		try:
			cmd = ["./../bdd_minisat_all-1.0.2/bdd_minisat_all_release", "{}{}".format(directory, file)]
			start_bdd = 0
			result = subprocess.run(cmd, capture_output=True, timeout=timeout_time)
			result = result.stdout.decode().splitlines()
			
			model_count = result[-1].split(" ")[-1]
			end_bdd = float(result[-12].split(" ")[-2])
			if sol == -1:
				sol = int(model_count)
			elif not sol == int(model_count):
				end_bdd = -1
				start_bdd = 0
		except subprocess.TimeoutExpired:
			start_bdd = 0
			end_bdd = timeout_time
			#timeout_bdd = True
		except IndexError:
			start_bdd = 0
			end_bdd = timeout_time + 1		
	
	tabularallsat_models = -1
	if(timeout_tabularallsat):
		start_tabularallsat = 0
		end_tabularallsat = timeout_time
	else:
		try:
			cmd = ["./../cdcl-vsids/solver", "{}{}".format(directory, file), "-c"]
			start_tabularallsat = 0
			result = subprocess.run(cmd, capture_output=True, timeout=timeout_time)
			result = result.stdout.decode().splitlines()
			model_count = result[8].split(" ")[-1]
			end_tabularallsat = float(result[-4].split(" ")[-2])
			tabularallsat_models = int(result[-18].split(" ")[-1])
			if sol == -1:
				sol = int(model_count)
			elif not sol == int(model_count):
				end_tabularallsat = -1
				start_tabularallsat = 0
		except subprocess.TimeoutExpired:
			start_tabularallsat = 0
			end_tabularallsat = timeout_time
			# timeout_tabularallsat = True


	outputfile = open(outputfilename, "a")
	outputfile.write("{},{},{},{},{},{},{}\n".format(file, sol, end_bc - start_bc, end_nbc - start_nbc, end_bdd - start_bdd,end_tabularallsat - start_tabularallsat, tabularallsat_models))
	outputfile.close()

	
