#/bin/csh 
export DIGEST=./Irys-scaffolding/KSU_bioinfo_lab/assemble_XeonPhi/third-party/fa2cmap_multi.pl
export SEW=./Irys-scaffolding/KSU_bioinfo_lab/stitch/sewing_machine.pl
export REPORT=./Irys-scaffolding/KSU_bioinfo_lab/assemble_XeonPhi/write_report.pl

export PREFIX=canu91x_high_outcov100
export OPTMAP=./vu_162_180K.cmap
export REF=/home/stelo/cowpea_assemblies/91x/$PREFIX/$PREFIX.fasta
export REFALIGNER=./tools/RefAligner
export SILICOMAP=/home/stelo/cowpea_assemblies/91x/$PREFIX/$PREFIX'_'BspQI.cmap
export OUT=/home/stelo/cowpea_assemblies/91x/$PREFIX_superscaffolds

#perl $DIGEST -v -i $REF -e BspQI
perl $SEW -o $OUT -g $OPTMAP -p $PREFIX -e BspQI -f $REF -r $SILICOMAP -a $REFALIGNER
#perl $REPORT -o $OUT -g $OPTMAP -p $PREFIX -e BspQI -f $REF -r $SILICOMAP -a $REFALIGNER --alignment_parameters default_alignment 
