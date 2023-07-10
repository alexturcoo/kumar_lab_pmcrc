#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=300000
#SBATCH --time=08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=blatmerge
#SBATCH --output=blatmerge.o
#SBATCH --error=blatmerge.e
#SBATCH --qos=privileged

source $HOME/virtenv/bin/activate
module load scipy-stack

# starting with the X chromo = 2023/04/18
python adding_Blatmerges_230418.py \
/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/Blat/results/chrX.Blatresults_full.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/SVlen/variants_freeze4_sv_insdel_alt.vcf.SVlength.csv chrX \
>> ./sv/adding_Blatmerges_230418.py.txt
echo "Counting for Chromo X"
rm /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/Blat/results/chrX.Blatresults_full.csv


# Running on Y chromo
python adding_Blatmerges_230418.py \
/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/Blat/results/chrY.Blatresults_full.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/SVlen/variants_freeze4_sv_insdel_alt.vcf.SVlength.csv chrY \
>> ./sv/adding_Blatmerges_230418.py.txt
echo "Counting for Chromo Y"
rm /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/Blat/results/chrY.Blatresults_full.csv

# runnning through autosomes
for i in {1..22};
do
	python adding_Blatmerges_230418.py \
/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/Blat/results/chr${i}.Blatresults_full.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/SVlen/variants_freeze4_sv_insdel_alt.vcf.SVlength.csv chr${i} \
>> ./sv/adding_Blatmerges_230418.py.txt
	echo "Counting for Chromo "$i""
	rm /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/Blat/results/chr${i}.Blatresults_full.csv
done
