#cellranger count --sample AA1_CKDL240033871-1A_22TJ5TLT3 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
#                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FQ_2 \
#                 --id Male_6mo_runner_1 \
#                 --expect-cells 10000 \
#                 --localmem 32 \
#                 --localcores 6 \
#                 --create-bam false

cellranger count --sample AA2_CKDL240033871-1A_22TJ5TLT3 \
                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FQ_2 \
                 --id Male_6mo_runner_2 \
                 --expect-cells 10000 \
                 --localmem 32 \
                 --localcores 6 \
                 --create-bam false

cellranger count --sample AA3_CKDL240033871-1A_22TJ5TLT3 \
                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FQ_2 \
                 --id Male_6mo_sed_1 \
                 --expect-cells 10000 \
                 --localmem 32 \
                 --localcores 6 \
                 --create-bam false

cellranger count --sample AA4_CKDL240033871-1A_22TJ5TLT3 \
                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FQ_2 \
                 --id Male_6mo_sed_2 \
                 --expect-cells 10000 \
                 --localmem 32 \
                 --localcores 6 \
                 --create-bam false

cellranger count --sample AA5_CKDL240033871-1A_22TJ5TLT3 \
                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FQ_2 \
                 --id Male_6mo_TFEB_1 \
                 --expect-cells 10000 \
                 --localmem 32 \
                 --localcores 6 \
                 --create-bam false

#cellranger count --sample AA6_CKDL240033871-1A_22TJ5TLT3 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
#                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FASTQ \
#                 --id Male_6mo_TFEB_2 \
#                 --expect-cells 10000 \
#                 --localmem 32 \
#                 --localcores 6 \
#                 --create-bam false

#cellranger count --sample AA7_CKDL240033872-1A_22TJ5TLT3 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
#                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FASTQ \
#                 --id Female_6mo_runner_1 \
#                 --expect-cells 10000 \
#                 --localmem 32 \
#                 --localcores 6 \
#                 --create-bam false
#
#cellranger count --sample AA8_CKDL240033872-1A_22TJ5TLT3 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
#                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FASTQ \
#                 --id Female_6mo_runner_2 \
#                 --expect-cells 10000 \
#                 --localmem 32 \
#                 --localcores 6 \
#                 --create-bam false
#
#cellranger count --sample AA9_CKDL240033872-1A_22TJ5TLT3 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
#                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FASTQ \
#                 --id Female_6mo_sed_1 \
#                 --expect-cells 10000 \
#                 --localmem 32 \
#                 --localcores 6 \
#                 --create-bam false
#
#cellranger count --sample AA10_CKDL240033872-1A_22TJ5TLT3 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
#                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FASTQ \
#                 --id Female_6mo_sed_2 \
#                 --expect-cells 10000 \
#                 --localmem 32 \
#                 --localcores 6 \
#                 --create-bam false
#
#cellranger count --sample AA11_CKDL240033872-1A_22TJ5TLT3 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
#                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FASTQ \
#                 --id Female_6mo_TFEB_1 \
#                 --expect-cells 10000 \
#                 --localmem 32 \
#                 --localcores 6 \
#                 --create-bam false
#
#cellranger count --sample AA12_CKDL240033872-1A_22TJ5TLT3 \
#                 --transcriptome /home/benayoun/Softwares/cellranger-9.0.0/Genome/refdata-gex-GRCm39-2024-A \
#                 --fastqs /home/benayoun/Projects/2024-10-07_10x_Cortes_lab/FASTQ \
#                 --id Female_6mo_TFEB_2 \
#                 --expect-cells 10000 \
#                 --localmem 32 \
#                 --localcores 6 \
#                 --create-bam false

