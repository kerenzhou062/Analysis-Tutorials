#!/usr/bin/sh

# Created on Nov 4, 2019
# @author: Keren Zhou
# @lab: Jianjun Chen Lab, City of Hope

if [[ $# -ne 3 ]]
then
  echo 'Usage: '$0' <CTK_sub_CIMS.pool.bed> <pool.tag.uniq.mutation.txt> <ouput.CT.bed>'
  exit 2
fi

function c2tMutation {
  awk 'BEGIN{OFS="\t";FS="\t";}
  ARGIND==1{
    key = $1"\t"$2"\t"$3"\t"$6;
    arrA[key] = 1;
  }
  ARGIND==2{
    key = $1"\t"$2"\t"$3"\t"$6;
    if (key in arrA) {
      if ($6 == "+") {
        if ($8=="C" && $10=="T") {
          arrB[key] += 1;
        }
      }else{
        if ($8=="G" && $10=="A") {
          arrB[key] += 1;
        }
      }
    }
  }
  ARGIND==3{
    key = $1"\t"$2"\t"$3"\t"$6;
    if (key in arrB) {
      $5 = arrB[key];
      print $0;
    }
  }' $1 $2 $1 | sort -t $'\t' -k1,1 -k2,2n > $3
}

c2tMutation $1 $2 $3
