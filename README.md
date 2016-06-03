#how2use_try  
gcc -o nonsymm nonsymm.c  
./nonsymm 10000 10000 10 > data.txt  
mkdir memo  
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 30 < data.txt)  

