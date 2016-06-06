#how2use_try  

gcc -o sparse_nonsymm sparse_nonsymm.c  
./sparse_nonsymm 10000 1000 20 > data.txt  
make clean  
make  
#下の引数は、得る特異対, s(sparse) or d(dense),accuracy = 10^(-x)
./fileinput_gkl 30 s 7 < sparse_data.txt
