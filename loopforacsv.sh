mkdir $@.dir
for arg in 1 2 3 4 5 6
do
 cat $@ | (awk -v arg=${arg} '{if ($1 == arg) print;}') > $@.dir/${arg}.csv
done
