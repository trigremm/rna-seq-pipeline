dir="./ensembl_GRCh37_75/"
 /home/adminrig/anaconda3/envs/rna/bin/STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir ${dir} \
  --genomeFastaFiles ${dir}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
  --sjdbGTFfile ${dir}/Homo_sapiens.GRCh37.75.gtf \
  --sjdbOverhang 100 \
  --limitGenomeGenerateRAM 234000000000
