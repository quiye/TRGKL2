import sys
rgvs = sys.argv
x = rgvs[1]
y = rgvs[2]
print ("#!/bin/bash")
print ("#============ PBS Options ============")
print ("#QSUB -q gr20100a")
print ("#QSUB -ug gr20114")
print ("#QSUB -W 48:00")
print ("#QSUB -A p=1:t=68:c=68:m=90G")
print ("#============ Shel Script ============")
print ("set -x")
print ("module load intel")
print ("ulimit -s unlimited")
print ("export OMP_STACKSIZE=10G")
print ("aprun -n $QSUB_PROCS -d $QSUB_THREADS -N $QSUB_PPN ./fileinput_gkl "+str(y)+" "+str(x)+"0 s l 100 43933 321961 inputmatrix > " +str(x)+ "0_"+str(y)+".txt")
print ("")

