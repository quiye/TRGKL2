#how2use_try  
gcc -o sparse_nonsymm sparse_nonsymm.c  
./sparse_nonsymm 10000 1000 20 > data.txt  
mkdir memo  
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 30 s < data.txt) 
