
'''
#input: l1, l2 [instance1, instance2, instance3, .....]
#output: output.csv
'''

import argparse
import time
import numpy
import multiprocessing
import subprocess
import shlex
import re


# Config path
app_path = "/home/luan/project/IVMP/"
instances_path = "/home/luan/project/IVMP/instances_new/"
output_path="/home/luan/project/IVMP/output/gupta_conf4/"

input_inst = ""

verbose = False
write_flag = True
combination = False

def read_args():
	parser = argparse.ArgumentParser()

	parser.add_argument('--inst', nargs='+')
	parser.add_argument('--input_inst', default='instances.txt')
	parser.add_argument('--repeat', nargs='?',  type=int)
	parser.add_argument('--gupta', help='Print more data',  action='store_true')
	parser.add_argument('--gupta_comb', help='Print more data',  action='store_true')
	parser.add_argument('--ils', help='Print more data',  action='store_true')
	parser.add_argument('--write', help='write output in a file', action='store_true')
	parser.add_argument('--verbose', help='Print more data', action='store_true')
	parser.add_argument('--output_path', default=".")

	global verbose, write_flag, output_path, combination, input_inst
	
	instances = parser.parse_args().inst
	repeat = parser.parse_args().repeat
	gupta_flag = parser.parse_args().gupta
	ils_flag = parser.parse_args().ils

	write_flag =  parser.parse_args().write

	verbose =  parser.parse_args().verbose
	output_path = parser.parse_args().output_path
	combination = parser.parse_args().gupta_comb

	input_inst = parser.parse_args().input_inst

	# print verbose, write_flag, output_path
	
	if repeat is None:
		repeat = 1

	return instances, repeat, gupta_flag, ils_flag

def write_file(output, path, file_name):
	#print path + file_name
	with open(path + file_name, "w") as f:
		f.write(output)

def prepare_cmd(instance, parameters = [], algorithm ='i'):
	return app_path + "interference " +  instance + "  " + algorithm

def exec_command(cmd):
	p = None
	
	try:
		p = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	except:
		print "error"

	return p 

def top_msg(msg):
	out = 20 * "#" 
	out +=  "\n\t"+msg + "\n"
	out += 20 * "#"
	return out

def print_warning(msg):
	if verbose:
		print msg

def exec_gupta(instances, l1 = None, l2 = None):
	instances.sort()
	
	out_file = ""
	for instance in instances:
		output, error = exec_command(prepare_cmd(instance, algorithm="g"))
		
		print_warning(top_msg("Output"))
		print_warning(output)

		if write_flag:
			with open(output_path+instance+"_out", "w") as f:
				f.write(output)

		num_machines = int(output.splitlines()[4].split()[-1])
		interference = float(output.splitlines()[6].split()[-1])

		print_warning( instance + ", " + str(interference) + ", " + str(num_machines))
		out_file += instance + ", " + str(interference) + ", " + str(num_machines) + "\n"

	write_file(out_file, "", "gupta.csv")

def exec_ils(instances, repeat):
	instances.sort()
	out_file = ""
	for instance in instances:
		interference = []
		machines = []
		for i in range(repeat):
			output, error = exec_command(prepare_cmd(instance, algorithm="i"))
			time.sleep(1)
			if output:
				print_warning(top_msg("Output"))
				print_warning(output)
			else :
				print_warning(top_msg("Output"))
				print_warning(error)
			if write_flag:
				with open(output_path+instance+"_out", "a") as f:
					f.write(output + "\n")
			#  Computa a media
			machines.append(int(output.splitlines()[4].split()[-1]))
			interference.append(float(output.splitlines()[6].split()[-1]))
		# statistic
		avg_interference = numpy.mean(interference)
		avg_numMachine = numpy.mean(machines)
		std_interference = numpy.std(interference)
		std_machines = numpy.std(machines)
		if write_flag:
				with open(output_path+instance+"_out", "a") as f:
					f.write("\n" + str(avg_interference) + "  +-" + str(std_interference) + " " + str(avg_numMachine) + "\n")
		print_warning( instance + ", " + str(avg_interference) + ", " + str(avg_numMachine))
		out_file += instance + ", " + str(avg_interference) +  ", " + str(std_interference) + ", " + str(avg_numMachine) + "\n"
	write_file(out_file, "", "ils.csv")

def exec_gupta_combinations(instances):
	pass

def exec_sequential(repeat, i_nfc, i_fwidth) : 
	# Defining the programs to run
	versions = []
	versions.append("sequential-scalar")
	versions.append("sequential-auto")
	versions.append("sequential-auto-optimized")
	versions.append("sequential-manual")
	exec_batch(versions, repeat, nfc = i_nfc, fwidth = i_fwidth) 

def exec_parallel_omp(repeat, i_nfc, i_fwidth, i_nt) :
	# Defining the programs to run
	versions = []
	versions.append("openmp-scalar")
	versions.append("openmp-auto")
	versions.append("openmp-manual")
	exec_batch(versions, repeat, nfc = i_nfc, fwidth = i_fwidth, nt = i_nt)

def exec_parallel_mpi(repeat, i_nfc, i_fwidth,  i_np) :
	# Defining the programs to run
	versions = []
#	versions.append("mpi-scalar")
#	versions.append("mpi-auto")
	versions.append("mpi-auto-optimized")
#	versions.append("mpi-manual")
	exec_batch(versions, repeat, nfc = i_nfc, fwidth = i_fwidth, np = i_np)

def exec_hybrid(repeat, i_nfc, i_fwidth,  i_np, i_nt) :
	# Defining the programs to run
	versions = []
	versions.append("mpi+openmp-scalar")
	versions.append("mpi+openmp-auto")
	versions.append("mpi+openmp-manual")
	exec_batch(versions, repeat, nfc = i_nfc, fwidth = i_fwidth, np = i_np, nt = i_nt)

def exec_batch(versions, repeat,  nfc = 16, fwidth = 5, np = 1, nt = 1) :
	out_file = ""
	for version in versions :
		exec_times = []
		command = "/usr/bin/time -f '%e' ./run-" + version + ".sh"
		command += " " + str(nfc)
		command += " " + str(fwidth)
		if (np != 1) : command += " " + str(np)
		if (nt != 1) : command += " " + str(nt)
		print command
		for i in range(repeat):
			print version + " " + str(i)
			output, error = exec_command(command)
			m = re.search('^[0-9]*[.][0-9]*$', error, re.M)
			exec_times.append(float(m.group(0)))
			#with open(output_path + "/error_output", "w") as f : f.write(error)
			#output, error = exec_command("awk '/^[0-9]*.[0-9]*$/{print $1};' error_output")
			#exec_times.append(float(output.split()[0]))
			if (i == 0) :
				with open("./kirchhoff_executions.csv", "a") as f : 
					f.write(version + " (PxT=" + str(np) + "x" + str(nt) + ")," + str(exec_times[i]))
			else :
				with open("./kirchhoff_executions.csv", "a") as f : f.write("," + str(exec_times[i]))
		with open("./kirchhoff_executions.csv", "a") as f : f.write("\n")
		avg_exec_times = numpy.mean(exec_times)
		std_exec_times = numpy.std(exec_times)
		print "Average of execution times: " + str(avg_exec_times)
		print "Standart deviation: " + str(std_exec_times)
		out_file += version + ", " + str(avg_exec_times) +  ", " + str(std_exec_times) + "\n"
	#write_file(out_file, "", "./kirchhoff.csv")
	with open("./kirchhoff.csv", "a") as f : f.write(out_file)

# --- Main Program ---

versions = []
instances, repeat, gupta_flag, ils_flag = read_args()

#output, error = exec_command("/opt/intel/parallel_studio_xe_2018/bin/psxevars.sh")

# Cleaning files
with open(output_path + "/kirchhoff_executions.csv", "w") as f : f.write("")
with open(output_path + "/kirchhoff.csv", "w") as f : f.write("")

# Execution in batch
with open("./kirchhoff_executions.csv", "a") as f : f.write("NFC=8, FWIDTH=5\n")
with open("./kirchhoff.csv", "a") as f : f.write("NFC=8, FWIDTH=5\n")
#exec_sequential(repeat, 8, 5)
#exec_parallel_omp(repeat, 8, 5, 2)
exec_parallel_mpi(repeat, 8, 5, 2)
#exec_parallel_omp(repeat, 8, 5, 3)
exec_parallel_mpi(repeat, 8, 5, 3)
#exec_parallel_omp(repeat, 8, 5, 4)
exec_parallel_mpi(repeat, 8, 5, 4)
#exec_parallel_omp(repeat, 8, 5, 5)
exec_parallel_mpi(repeat, 8, 5, 5)
#exec_parallel_omp(repeat, 8, 5, 6)
exec_parallel_mpi(repeat, 8, 5, 6)
#exec_hybrid(repeat, 8, 5, 2, 3)
#exec_hybrid(repeat, 8, 5, 3, 2)

with open("./kirchhoff_executions.csv", "a") as f : f.write("NFC=16, FWIDTH=3\n")
with open("./kirchhoff.csv", "a") as f : f.write("NFC=16, FWIDTH=3\n")
#exec_sequential(repeat, 16, 3)
#exec_parallel_omp(repeat, 16, 3, 2)
exec_parallel_mpi(repeat, 16, 3, 2)
#exec_parallel_omp(repeat, 16, 3, 3)
exec_parallel_mpi(repeat, 16, 3, 3)
#exec_parallel_omp(repeat, 16, 3, 4)
exec_parallel_mpi(repeat, 16, 3, 4)
#exec_parallel_omp(repeat, 16, 3, 5)
exec_parallel_mpi(repeat, 16, 3, 5)
#exec_parallel_omp(repeat, 16, 3, 6)
exec_parallel_mpi(repeat, 16, 3, 6)
#exec_hybrid(repeat, 16, 3, 2, 3)
#exec_hybrid(repeat, 16, 3, 3, 2)

with open("./kirchhoff_executions.csv", "a") as f : f.write("NFC=32, FWIDTH=1\n")
with open("./kirchhoff.csv", "a") as f : f.write("NFC=32, FWIDTH=1\n")
#exec_sequential(repeat, 32, 1)
#exec_parallel_omp(repeat, 32, 1, 2)
exec_parallel_mpi(repeat, 32, 1, 2)
#exec_parallel_omp(repeat, 32, 1, 3)
exec_parallel_mpi(repeat, 32, 1, 3)
#exec_parallel_omp(repeat, 32, 1, 4)
exec_parallel_mpi(repeat, 32, 1, 4)
#exec_parallel_omp(repeat, 32, 1, 5)
exec_parallel_mpi(repeat, 32, 1, 5)
#exec_parallel_omp(repeat, 32, 1, 6)
exec_parallel_mpi(repeat, 32, 1, 6)
#exec_hybrid(repeat, 32, 1, 2, 3)
#exec_hybrid(repeat, 32, 1, 3, 2)

#exec_batch(versions, repeat)
