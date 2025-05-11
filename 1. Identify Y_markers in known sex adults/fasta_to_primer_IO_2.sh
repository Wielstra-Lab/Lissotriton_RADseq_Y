
fasta=$1
header=$2
out_IO=$3

cat $header >> $out_IO

cat $fasta | while read seq_ID Info; read sequence;
do seq_num=${seq_ID#*>*};
  short_seq=${sequence%*NNNNNNNNNN*};
  echo SEQUENCE_ID=$seq_num-short >> $out_IO
  echo SEQUENCE_TEMPLATE=$short_seq >> $out_IO
  echo = >> $out_IO
  echo SEQUENCE_ID=$seq_num-long >> $out_IO
  echo SEQUENCE_TEMPLATE=$sequence >> $out_IO
  echo SEQUENCE_TARGET=125,50 >> $out_IO
  echo = >> $out_IO
  done
