#/bin/csh
#export REFALIGNER=./tools/RefAligner # on h4
export REFALIGNER=./archive/tools_SSE/RefAligner # on h3
export DIGEST=./Irys-scaffolding/KSU_bioinfo_lab/assemble_XeonPhi/third-party/fa2cmap_multi.pl
export OPTMAP=./vu_162_180K.cmap
export REPO=/home/stelo/cowpea_assemblies/91x
export THREADS=12

export PREFIX=abruijn49x_corrected_k19_cov49_nochi
export REF=$REPO/$PREFIX/$PREFIX.fasta
export SILICOMAP=$REPO/$PREFIX/$PREFIX'_'BspQI.cmap
export OUT=$REPO/$PREFIX/$PREFIX'_'BNG_vs_seq
perl $DIGEST -v -i $REF -e BspQI
$REFALIGNER -f -ref $OPTMAP -i $SILICOMAP -o $OUT -endoutlier 1e-3 -outlier 1e-4 -extend 1 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 2.9 -resSD 0.75 -mres 1.2 -A 5 -biaswt 0 -M 3 -Mfast 0 -maxmem 256 -maxthreads $THREADS -deltaX 9 -deltaY 9 -xmapchim 14 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-15 -stdout -stderr -xmaplen -indel

export PREFIX=canu91x_high_default_nochi
export REF=$REPO/$PREFIX/$PREFIX.fasta
export SILICOMAP=$REPO/$PREFIX/$PREFIX'_'BspQI.cmap
export OUT=$REPO/$PREFIX/$PREFIX'_'BNG_vs_seq
perl $DIGEST -v -i $REF -e BspQI
$REFALIGNER -f -ref $OPTMAP -i $SILICOMAP -o $OUT -endoutlier 1e-3 -outlier 1e-4 -extend 1 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 2.9 -resSD 0.75 -mres 1.2 -A 5 -biaswt 0 -M 3 -Mfast 0 -maxmem 256 -maxthreads $THREADS -deltaX 9 -deltaY 9 -xmapchim 14 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-15 -stdout -stderr -xmaplen -indel

export PREFIX=falcon49x_runA_nochi
export REF=$REPO/$PREFIX/$PREFIX.fasta
export SILICOMAP=$REPO/$PREFIX/$PREFIX'_'BspQI.cmap
export OUT=$REPO/$PREFIX/$PREFIX'_'BNG_vs_seq
perl $DIGEST -v -i $REF -e BspQI
$REFALIGNER -f -ref $OPTMAP -i $SILICOMAP -o $OUT -endoutlier 1e-3 -outlier 1e-4 -extend 1 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 2.9 -resSD 0.75 -mres 1.2 -A 5 -biaswt 0 -M 3 -Mfast 0 -maxmem 256 -maxthreads $THREADS -deltaX 9 -deltaY 9 -xmapchim 14 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-15 -stdout -stderr -xmaplen -indel

export PREFIX=quiver91x_low_default_nochi
export REF=$REPO/$PREFIX/$PREFIX.fasta
export SILICOMAP=$REPO/$PREFIX/$PREFIX'_'BspQI.cmap
export OUT=$REPO/$PREFIX/$PREFIX'_'BNG_vs_seq
perl $DIGEST -v -i $REF -e BspQI
$REFALIGNER -f -ref $OPTMAP -i $SILICOMAP -o $OUT -endoutlier 1e-3 -outlier 1e-4 -extend 1 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 2.9 -resSD 0.75 -mres 1.2 -A 5 -biaswt 0 -M 3 -Mfast 0 -maxmem 256 -maxthreads $THREADS -deltaX 9 -deltaY 9 -xmapchim 14 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-15 -stdout -stderr -xmaplen -indel

export PREFIX=quiver91x_low_outcov100_nochi
export REF=$REPO/$PREFIX/$PREFIX.fasta
export SILICOMAP=$REPO/$PREFIX/$PREFIX'_'BspQI.cmap
export OUT=$REPO/$PREFIX/$PREFIX'_'BNG_vs_seq
perl $DIGEST -v -i $REF -e BspQI
$REFALIGNER -f -ref $OPTMAP -i $SILICOMAP -o $OUT -endoutlier 1e-3 -outlier 1e-4 -extend 1 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 2.9 -resSD 0.75 -mres 1.2 -A 5 -biaswt 0 -M 3 -Mfast 0 -maxmem 256 -maxthreads $THREADS -deltaX 9 -deltaY 9 -xmapchim 14 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-15 -stdout -stderr -xmaplen -indel

