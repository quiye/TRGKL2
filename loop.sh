for x in `ls | grep .txt`; do
echo $x
cat $x | grep " TIME "
done

