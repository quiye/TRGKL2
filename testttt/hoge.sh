mkdir pdfs
for x in $( ls . | grep .txt$ ); do
    d_name=${x%.txt}
    mkdir ${d_name}
    for arg in 1 3 5; do
        cat ${x} | (awk -v arg=${arg} '{if ($1 == arg) print;}') > ${d_name}/${arg}.csv
    done
    cd ${d_name}/
    mv 1.csv QR_1.0
    mv 3.csv OQDS
    mv 5.csv One-Sided_Jacobi
    gnuplot -e "dirname='${d_name}'" ../gpshell.gp
    pdftk gnuplot-pdf.pdf cat 3 output ../pdfs/${d_name}.pdf
    rm gnuplot-pdf.pdf
    cd ..
done

