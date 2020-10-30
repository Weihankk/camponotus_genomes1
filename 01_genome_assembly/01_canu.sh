qlogin -q omni -P quanah -pe sm 2

cd /lustre/scratch/jmanthey
mkdir 13_camponotos_genome
cd 13_camponotos_genome

# this one didn't work
/lustre/work/jmanthey/canu-1.7.1/Linux-amd64/bin/canu -p camponotus -d camponotus-pb \
genomeSize=330000000 -pacbio-raw /home/jmanthey/camponotus/camponotus_part?.pb.fasta gnuplotTested=true \
gridEngineMemoryOption="-l h_vmem=MEMORY" gridEngineThreadsOption="-pe sm THREADS" \
gridOptions="-P quanah -q omni" correctedErrorRate=0.065 corMhapSensitivity=normal

# this one worked but ended up using the modifications in canu 1.9 below
/home/jmanthey/canu-1.8/Linux-amd64/bin/canu -p camponotus -d camponotus-pb \
genomeSize=330000000 -pacbio-raw /home/jmanthey/camponotus/camponotus_part?.pb.fasta \
gridEngineMemoryOption="-l h_vmem=MEMORY" gridEngineThreadsOption="-pe sm THREADS" \
gridOptions="-P quanah -q omni" correctedErrorRate=0.065 corMhapSensitivity=normal


# correct for AT bias and increase target coverage
/home/jmanthey/canu-1.9/Linux-amd64/bin/canu -p camp -d camp2 \
genomeSize=330000000 -pacbio-raw /home/jmanthey/camponotus/camponotus_part?.pb.fasta \
minReadLength=900 \
gridEngineResourceOption="-l h_vmem=MEMORY -pe sm THREADS" \
gridOptions="-P quanah -q omni" correctedErrorRate=0.045 corMhapSensitivity=normal \
corMaxEvidenceErate=0.15 corOutCoverage=50
