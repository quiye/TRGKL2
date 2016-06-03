./sparse_nonsymm 1000000 10000 10 > sparse_data.txt
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 10 s < sparse_data.txt) > 10-10.txt
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 20 s < sparse_data.txt) > 10-20.txt
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 30 s < sparse_data.txt) > 10-30.txt
./sparse_nonsymm 1000000 10000 20 > sparse_data.txt
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 10 s < sparse_data.txt) > 20-10.txt
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 20 s < sparse_data.txt) > 20-20.txt
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 30 s < sparse_data.txt) > 20-30.txt
./sparse_nonsymm 1000000 10000 30 > sparse_data.txt
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 10 s < sparse_data.txt) > 30-10.txt
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 20 s < sparse_data.txt) > 30-20.txt
make clean > memo/logofclean.txt && make > memo/logofmake.txt && (./fileinput_gkl 30 s < sparse_data.txt) > 30-30.txt
rm sparse_data.txt
