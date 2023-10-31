#!/bin/bash




## This script takes in a bam file and returns a .hic file.

usage() {
  echo "Usage: createHiCfile.sh [option]"
  echo "          [-b <bam>]"
  echo "          [-n <name>]"
  echo "          [-g <genome>]"
  echo "          [-s <sorted>]"
  echo "          [-c <chr retention>]"
  echo "          [-j Juicer Path]"
  echo "          [-h help]"
  exit 1
 }


bam=""
id=""
genome=""
sorted=false
chr=false

while getopts "b:n:g:j:sch" o; do
    case "${o}" in
        b)
            bam=${OPTARG}
            ;;
        n)
            id=${OPTARG}
            ;;
        g)
            genome=${OPTARG}
            ;;
        j)
            juicer_path=${OPTARG}
            ;;
        s)
            sorted=true
            ;;
        c)  chr=true
            ;;
        h)
            usage
      ;;
        *)
            usage
            ;;
    esac
done


if [[  ! -s $bam ]]
then
  echo "ERROR: $bam doesn't exist or is empty." >&2
  exit 1;
fi

if [[ $id == "" ]]
then
  id=${bam:0:(${#bam}-4)}
fi

if [[ $genome == "" ]]
then
  genome="hg38"
fi

newBam=${bam}

if [[ $sorted != "true" ]]; then
## sort bam by name:
samtools sort -n $bam > $bam.sorted
newBam=${bam}.sorted
fi


##extract read info to make .hic file

if [[ $chr != "true" ]]; then
## removes "chr" substring
samtools view ${newBam} | gawk 'BEGIN {FS="\t"; OFS="\t"} {name1=$1; str1=and($2,16); chr1=substr($3, 4); pos1=$4; mapq1=$5; getline; name2=$1; str2=and($2,16); chr2=substr($3, 4); pos2=$4; mapq2=$5; sub(/\/1\>/, "", name1); sub(/\/2\>/, "", name2); if(name1==name2) { if (chr1>chr2){print name1, str2, chr2, pos2,1, str1, chr1, pos1, 0, mapq2, mapq1} else {print name1, str1, chr1, pos1, 0, str2, chr2, pos2 ,1, mapq1, mapq2}} else {getline}}' > ${id}.Arrowhead.input

else

samtools view ${newBam} | gawk 'BEGIN {FS="\t"; OFS="\t"} {name1=$1; str1=and($2,16); chr1=$3; pos1=$4; mapq1=$5; getline; name2=$1; str2=and($2,16); chr2=$3; pos2=$4; mapq2=$5; sub(/\/1\>/, "", name1); sub(/\/2\>/, "", name2); if(name1==name2) { if (chr1>chr2){print name1, str2, chr2, pos2,1, str1, chr1, pos1, 0, mapq2, mapq1} else {print name1, str1, chr1, pos1, 0, str2, chr2, pos2 ,1, mapq1, mapq2}} else {getline}}' > ${id}.Arrowhead.input
fi

## sort file by columns 3 and 7

sort -k3,3d -k7,7d < ${id}.Arrowhead.input >  ${id}.Arrowhead.input.sorted


## create .hic file from 11 column file
java -jar ${juicer_path} pre -q 10 ${id}.Arrowhead.input.sorted ${id}.hic ${genome}
rm ${id}.Arrowhead.input ${id}.Arrowhead.input.sorted
