# input bed file
f=$1
perl -ane 'if($_ !~ /^chr/){print "chr$_";}else{print $_;}' $f > $f.wchr
liftOver $f.wchr hg19ToHg38.over.chain $f.wchr.liftOver.build38 $f.wchr.liftOver.build38.unmapped
