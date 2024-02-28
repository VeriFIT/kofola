import time 
import subprocess

with open("bp_programs.txt") as file:
    progs = [line.rstrip() for line in file if line.strip() != ""]

print("| Program |", end="", flush=True)
print("")
print("-----------------------------------------------------------------------")

for p in progs:
    print("|  " + p + "  |", end="", flush=True)
    startTime = time.time()
    try:
        result = subprocess.run(["../build/src/kofola", "--hyper", "./bp/" + p,  "./bp/gni.txt"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)

        endTime = time.time()
        et = endTime - startTime 

        out = result.stdout.decode("utf-8").strip()
        print(out, end=" ")
        if out == "SAT" or out == "UNSAT":
            print("  " + "%.2f" % et + "  |", end="", flush=True)
        else:
            print("   ERR  |", end="", flush=True)
            print(out)
    except subprocess.TimeoutExpired:
        print("   TO  |", end="", flush=True)


        

    print("")
