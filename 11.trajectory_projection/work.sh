
python3 classify_trajectory_genes.py
python3 extract_TLC_exp.py
python3 trajectory_mapping.py


sort -u /dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/17.downstream/19.zebrafish_trajectory/monocle2/genSmoothCurves_genes_to_anno.orig.txt > zebrafish_trajectory_genes.txt

awk '{print $6}' exp.trajectory_ordered.txt | sort -u | fishInWinter.pl -bf table -ff table -fc 2 - /dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/17.downstream/13.sonic/styela_to_zebraname.txt | cut -f 1 | sort -u | comm -12 - zebrafish_trajectory_genes.txt > comm_trajectory_genes.txt

