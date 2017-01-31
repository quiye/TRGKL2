for x in 1 2 3 ; do
for y in 1 2 3 4 ; do
python write.py $x $y > job_$x$y
done
done
ls -l | awk '{print "qsub " $9}'|grep job_ > jobs
