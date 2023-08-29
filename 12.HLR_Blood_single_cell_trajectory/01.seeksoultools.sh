# http://seeksoul.seekgene.com/zh/v1.2.0/2.tutorial/1.rna/2.run.html
date

mkdir -p HLR_downsp/
cd HLR_downsp/
seeksoultools rna run \
        --fq1 HLR_downsp_1.fq.gz \
        --fq2 HLR_downsp_2.fq.gz \
        --samplename HLR_downsp \
        --outdir HLR_downsp \
        --genomeDir S_clava_211028/star \
        --gtf S_clava_211028/genes/genes.gtf \
        --chemistry MM \
        --core 20 \
        --star_path cellranger-6.1.2/lib/bin/STAR
date
mkdir -p stcr-2/
cd stcr-2/
seeksoultools rna run \
        --fq1 YF230731-stcr-2_good_1.fq.gz \
        --fq2 YF230731-stcr-2_good_2.fq.gz \
        --samplename stcr-blood-2 \
        --outdir stcr-2 \
        --genomeDir /home/jianga/Project/01.SeekGeneSingleCellDong/01.Reference/S_clava_211028/star \
        --gtf S_clava_211028/genes/genes.gtf \
        --chemistry MM \
        --core 20 \
	--star_path cellranger-6.1.2/lib/bin/STAR

date
