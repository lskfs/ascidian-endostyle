# 127804648   159755810 11772916622 HLR_2.fq
# 127804648   159755810 11772916622 HLR_2.fq

ln -s Unknown_BN590-ZX01-010007_good_1.fq.gz HLR_1.fq.gz
ln -s Unknown_BN590-ZX01-010007_good_2.fq.gz HLR_2.fq.gz

seqtk sample -s100 HLR_1.fq 0.2 > HLR_downsp_1.fq
gzip HLR_downsp_1.fq
seqtk sample -s100 HLR_2.fq 0.2 > HLR_downsp_2.fq
gzip HLR_downsp_2.fq
