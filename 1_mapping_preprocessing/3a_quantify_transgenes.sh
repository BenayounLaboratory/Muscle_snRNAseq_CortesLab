#cat Cre_sequence.fasta Human_TFEB_sequence.fasta gencode.vM36.pc_transcripts.fa > 2025-03-11_Transgene_and_mouse_transcripts.fa
#
#kallisto index -i 2025-03-11_Transgene_and_mouse_transcripts.idx 2025-03-11_Transgene_and_mouse_transcripts.fa
#
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Male_6mo_runner_1.kallisto_res     --single -l 150 -s 50 ../FASTQ/M/AA1_CKDL240033871-1A_22TJ5TLT3_S4_L001_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Male_6mo_runner_2.kallisto_res     --single -l 150 -s 50 ../FASTQ/M/AA2_CKDL240033871-1A_22TJ5TLT3_S6_L001_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Male_6mo_sed_1.kallisto_res        --single -l 150 -s 50 ../FASTQ/M/AA3_CKDL240033871-1A_22TJ5TLT3_S5_L001_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Male_6mo_sed_2.kallisto_res        --single -l 150 -s 50 ../FASTQ/M/AA4_CKDL240033871-1A_22TJ5TLT3_S2_L001_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Male_6mo_TFEB_1.kallisto_res       --single -l 150 -s 50 ../FASTQ/M/AA5_CKDL240033871-1A_22TJ5TLT3_S1_L001_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Male_6mo_TFEB_2.kallisto_res       --single -l 150 -s 50 ../FASTQ/M/AA6_CKDL240033871-1A_22TJ5TLT3_S3_L001_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Female_6mo_runner_1.kallisto_res   --single -l 150 -s 50 ../FASTQ/F/AA7_CKDL240033872-1A_22TJ5TLT3_S3_L002_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Female_6mo_runner_2.kallisto_res   --single -l 150 -s 50 ../FASTQ/F/AA8_CKDL240033872-1A_22TJ5TLT3_S6_L002_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Female_6mo_sed_1.kallisto_res      --single -l 150 -s 50 ../FASTQ/F/AA9_CKDL240033872-1A_22TJ5TLT3_S5_L002_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Female_6mo_sed_2.kallisto_res      --single -l 150 -s 50 ../FASTQ/F/AA10_CKDL240033872-1A_22TJ5TLT3_S2_L002_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Female_6mo_TFEB_1.kallisto_res     --single -l 150 -s 50 ../FASTQ/F/AA11_CKDL240033872-1A_22TJ5TLT3_S1_L002_R2_001.fastq.gz
kallisto quant -i 2025-03-11_Transgene_and_mouse_transcripts.idx  --plaintext -o Female_6mo_TFEB_2.kallisto_res     --single -l 150 -s 50 ../FASTQ/F/AA12_CKDL240033872-1A_22TJ5TLT3_S4_L002_R2_001.fastq.gz
