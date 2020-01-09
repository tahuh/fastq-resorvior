# fastq-resorvior
Sampling Fastq reads using resorvior sampling

## System tested
Ubuntu 18.04

## Required compiler
gcc >= 4.8

## Required library

zlib

Heng Li's kseq.h (included in this repository)

## Compilation

After clone this respository do the following
```
$ cd fastq-resorvior
$ gcc -O2 -o fastq-resorvior fastq-resorvior.c -lz
```

## Details

```
./fastq-resorvior -m mode [options] in1.fastq[,in2.fastq] out1.fastq[,out2.fastq]

[options]

   -h          This help message
   -f   FLT    Factor of sampling. If -k is also set, then the smaller value using 
               this factor and -k result will be sampled [default=1.0]
   -k   INT    Actual number of reads to sample. If -f is also set, then the smaller 
               values using -f and -k will be used [default=all_reads]
   -m   STR    Sequencing read mode either s or p to denote single end and paried end
               [default=s]

Example command for single end read sampling

./fastq-resorvior in.fastq out.fastq

Example command for paired end read sampling

./fastq-resorvior in1.fastq in2.fastq out1.fastq out2.fastq

```

## Author

Thomas Sunghoon Heo

## License
GPL v3
