#!bin/bash 

env_dir=./conda_env

out_dir=./out_2Feb23

reads_dir=./dados_brutos/fastqs


mkdir -p ${out_dir}/fastqc

source activate ${env_dir}/fastqc

fastqc \
-o ${out_dir}/fastqc \
--casava \
${reads_dir}/*

conda deactivate 


#---------------------------------------------------------------------------
# Detect overrepresented sequences 
#---------------------------------------------------------------------------
# 1. Make a copy of sequences 
# 2. Decompress 
# 3. Combine all R1 sequences togather and all R2 (resulted in 2 files)
# 4. Run cut adapt on each file 

mkdir -p ${out_dir}/seq_data/fastq
mkdir -p ${out_dir}/seq_data/fastq_comb
mkdir -p ${out_dir}/fastqc_comb

cp ${reads_dir}/* ${out_dir}/seq_data/fastq/

cd ${out_dir}/seq_data/fastq/

gzip -d *.gz

cd - 

cat ${out_dir}/seq_data/fastq/*_R1_001.fastq > ${out_dir}/seq_data/fastq_comb/comb_R1.fastq

cat ${out_dir}/seq_data/fastq/*_R2_001.fastq > ${out_dir}/seq_data/fastq_comb/comb_R2.fastq

source activate ${env_dir}/fastqc

fastqc \
-o ${out_dir}/fastqc_comb \
${out_dir}/seq_data/fastq_comb/*

conda deactivate 


#---------------------------------------------------------------------------------------------
# Remove primers from forward and reverse reads
#---------------------------------------------------------------------------------------------
source activate ${env_dir}/cutadapt

mkdir -p /home/jovyan/work/pvc/Brazil/out_2Feb23/seq_data/no_prim_fastq_gz

for i in ${reads_dir}/*_R1_001.fastq.gz; do 

name=$(basename ${i%%"_L001_R"*})
name=${name##*$"-1_"}

echo $name

echo $path
cutadapt \
-g CCTACGGGNNGCAGCAG \
-G GGACTACNNGGGTNTCTAAT \
--discard-untrimmed \
--cores 60 \
-o /home/jovyan/work/pvc/Brazil/out_2Feb23/seq_data/no_prim_fastq_gz/br${name}_${name}_L001_R1_001.fastq.gz \
-p /home/jovyan/work/pvc/Brazil/out_2Feb23/seq_data/no_prim_fastq_gz/br${name}_${name}_L001_R2_001.fastq.gz \
${i} \
${i%%"_R"*}_R2_001.fastq.gz

done


#--------------------------------------------------------------------------------------------
# Check quality of with fastQC of trimmid reads
#--------------------------------------------------------------------------------------------
# Decompress files 
cd ${out_dir}/seq_data/no_prim_fastq

gzip -d *.gz 

cd - 


# Concatinate decompressed filtered files 
cat ${out_dir}/seq_data/no_prim_fastq/*_R1.fastq > ${out_dir}/seq_data/fastq_comb/filt_comb_R1.fastq

cat ${out_dir}/seq_data/no_prim_fastq/*_R2.fastq > ${out_dir}/seq_data/fastq_comb/filt_comb_R2.fastq


# Check with fastQC
source activate ${env_dir}/fastqc

fastqc \
-o ${out_dir}/fastqc_comb \
${out_dir}/seq_data/fastq_comb/filt_*

conda deactivate 