mkdir -p logs/map_failed
grep error slurm-242* | awk -F":" '{print $1}' | while read out;do mv $out logs/map_failed/.;done
mkdir -p logs/map
grep "Finished job" slurm-242* | awk -F":" '{print $1}' | while read out;do mv $out logs/map/.;done

grep output logs/map/*.out | awk -F"/" '{print $5}' > done.txt
grep "Real time" data/mapped_reads/*.bwamem | awk -F"[./]" '{print $3}' | sort > done_bwamem.txt
cat done.txt | awk -F"." '{print $1}' | sort > done_out.txt
diff done_bwamem.txt done_out.txt # AM4426
#390 done_bwamem.txt
#390 done_out.txt

ls -tral data/mapped_reads/*.bam | awk -F"/" '{print $3}' > found.txt
cat done.txt | grep -vFwf - found.txt | awk -F"." '{print $1}' > incomplete.txt
ls -tral data/mapped_reads/AM* | awk '{print $9}' > files.txt
cat incomplete.txt | grep -Fwf - files.txt | while read f;do rm $f;done

