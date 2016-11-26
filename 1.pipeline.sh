#/bin/csh
export DIGEST=./Irys-scaffolding/KSU_bioinfo_lab/assemble_XeonPhi/third-party/fa2cmap_multi.pl
export STITCH=./Irys-scaffolding/KSU_bioinfo_lab/stitch/stitch.pl
export PREFIX=quiver91x_low_default_tig0001
export OPTMAP=./vu_162_180K.cmap
export REF=./$PREFIX/$PREFIX.fasta
export REFNEW=./$PREFIX/$PREFIX'_'new.fasta
export REFNEWFORMATTED=./$PREFIX/$PREFIX'_'new'_'formatted.fasta
export REFSCAF=./$PREFIX/$PREFIX'_'superscaffold.fasta

export REFALIGNER=./tools/RefAligner
export SILICOMAP=./$PREFIX/$PREFIX'_'BspQI.cmap
export SILICOMAPNEW=./$PREFIX/$PREFIX'_'new'_'BspQI.cmap
export SILICOMAPSCAF=./$PREFIX/$PREFIX'_'superscaffold'_'BspQI.cmap

export OUT=./$PREFIX/$PREFIX'_'BNG_vs_seq
export OUTNEW=./$PREFIX/$PREFIX'_'new'_'BNG_vs_seq
export OUTSCAF=./$PREFIX/$PREFIX'_'superscaffold

export OUTR=./$PREFIX/$PREFIX'_'new'_'BNG_vs_seq_r.cmap
export XMAPMOD=./$PREFIX/$PREFIX'_'new'_'BNG_vs_seq_modified.xmap
export OUTSTITCH=./$PREFIX/$PREFIX'_'new'_'BNG_vs_seq_stitch

# 1. digest in-silico the assembly
perl $DIGEST -v -i $REF -e BspQI

# 2. run RefAligner to align the in-silico map to the optical map
$REFALIGNER -f -ref $OPTMAP -i $SILICOMAP -o $OUT -endoutlier 1e-3 -outlier 1e-4 -extend 1 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 2.9 -resSD 0.75 -mres 1.2 -A 5 -biaswt 0 -M 3 -Mfast 0 -maxmem 256 -maxthreads 20 -deltaX 9 -deltaY 9 -xmapchim 14 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-15 -stdout -stderr -xmaplen -indel

# 3. check visually the $OUT.cmap using IrysView for chimeric contigs, find the location
# create a CSV file $PREFIX_cut_list.csv with CONTIGID,LOCATION 
# split.py will split to contigs, then store the chimeric-free assembly as $REFNEW
# trust an alignment with confidence score > 15, especially when there are no other sequences aligned to the same region

cd ./$PREFIX
python ../split.py $PREFIX
rm *.flat
rm *.gdx
cd ..

# 4. digest in-silico the chimeric-free assembly
perl $DIGEST -v -i $REFNEW -e BspQI

# 5. run RefAligner to align the new in-silico map to the optical map
$REFALIGNER -f -ref $OPTMAP -i $SILICOMAPNEW -o $OUTNEW -endoutlier 1e-3 -outlier 1e-4 -extend 1 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 2.9 -resSD 0.75 -mres 1.2 -A 5 -biaswt 0 -M 3 -Mfast 0 -maxmem 256 -maxthreads 20 -deltaX 9 -deltaY 9 -xmapchim 14 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-15 -stdout -stderr -xmaplen -indel

# 6. check the alignments in $OUTNEW in IrysView and remove low alignments from the xmap -> store as $XMAPMOD

# 7. use the cleaned xmap to stich
fold -w 80 $REFNEW > $REFNEWFORMATTED
perl $STITCH -r $OUTR -x $XMAPMOD -f $REFNEW -o $OUTSTITCH --f_con 0 --f_algn 0 --s_con 0 --s_algn 0 -n 10000
mv ./$PREFIX/$PREFIX'_'new_BNG_vs_seq_stitch_superscaffold.fasta $REFSCAF

# 8. digest in-silico the superscaffold assembly
perl $DIGEST -v -i $REFSCAF -e BspQI

# 9. run Refaligner again, possibly removing more bad alignments
$REFALIGNER -f -ref $OPTMAP -i $SILICOMAPSCAF -o $OUTSCAF -endoutlier 1e-3 -outlier 1e-4 -extend 1 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 2.9 -resSD 0.75 -mres 1.2 -A 5 -biaswt 0 -M 3 -Mfast 0 -maxmem 256 -maxthreads 20 -deltaX 9 -deltaY 9 -xmapchim 14 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-15 -stdout -stderr -xmaplen -indel

# iterate 6-9 until done
