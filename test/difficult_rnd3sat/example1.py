import subprocess

for i in range(10, 51):
    for j in range(10):
        f = open("rnd3sat-diff-{}-{}.cnf".format(i,j), "w")
        bashCommand = "cnfgen randkcnf 3 {} {}".format(i,round(i*1.5))
        process = subprocess.Popen(bashCommand.split(), stdout=f)
        output, error = process.communicate()