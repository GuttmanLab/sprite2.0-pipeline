#bin/bash

barcode_arr=()



for fq in *unmapped.1.fastq.gz
do
	IFS=. read name lane barcode map rpe fq gz <<< $fq
	barcode_arr+=("$barcode")
uniq=($(printf "%s\n" "${barcode_arr[@]}" | sort -u))
done

#echo "${uniq[@]}"

#4_HL522DSXX.4.TTGTCTAT_TTAGCCAG.unmapped.1.fastq.gz
for i in "${uniq[@]}"  
do  
	echo "Combining R1:"
	echo *.${i}.unmapped.1.fastq.gz
	cat *.${i}.unmapped.1.fastq.gz > HL522DSXX.${i}_R1.fastq.gz
	echo "Combining R2:"
	echo cat *.${i}.unmapped.2.fastq.gz
	cat *.${i}.unmapped.2.fastq.gz > HL522DSXX.${i}_R2.fastq.gz
done  
