#ATAC correct Tobias - subsample and then perform tn5 shift
#load python/3.7.1
#load samtools
#Yuyan 11-09-19


START=$SECONDS  

MOTIF=$1
INDIR=$2
GENOME=$3
PEAK=$4
OUTDIR=$5
PREFIX=$6


TOBIAS BINDetect --motifs $MOTIF \
                 --signals $INDIR'/'$PREFIX*.sorted.footprint.bw \
                 --genome $GENOME \
                 --peaks $PEAK \
                 --outdir $OUTDIR --cores 4 \
                 --prefix $PREFIX

sleep 10
# hack to ensure job lasts 10 minutes to ensure no throttling
END=$SECONDS
ELAPSED=$((END-START))
echo $ELAPSED
if [ $ELAPSED -lt 600 ]; then
  TOSLEEP=$((600 - ELAPSED))
  sleep $TOSLEEP
fi

