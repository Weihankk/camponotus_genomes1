# use the Juicebox Assembly Tools to manually check and correct any obvious issues with the assembly using the
# assembly and hic files


# run the asm pipeline post review tool to input the manually edited assembly file and create the final
# fasta file and hic file for visualization

~/3d-dna/run-asm-pipeline-post-review.sh \
-r /lustre/scratch/jmanthey/13_camponotos_genome/scaffolding/final/camponotus_reordered.final.review.assembly \
--sort-output -i 25000 \
--build-gapped-map \
~/juicer/references/camponotus_reordered.fasta \
/lustre/scratch/jmanthey/13_camponotos_genome/scaffolding/aligned/merged_nodups.txt


# repeat

# use the Juicebox Assembly Tools to manually check and correct any obvious issues with the assembly using the
# assembly and hic files


# run the asm pipeline post review tool to input the manually edited assembly file and create the final
# fasta file and hic file for visualization

~/3d-dna/run-asm-pipeline-post-review.sh \
-r /lustre/scratch/jmanthey/13_camponotos_genome/scaffolding/final2/camponotus_reordered.final.review.assembly \
--sort-output -i 25000 \
--build-gapped-map \
~/juicer/references/camponotus_reordered.fasta \
/lustre/scratch/jmanthey/13_camponotos_genome/scaffolding/aligned/merged_nodups.txt
