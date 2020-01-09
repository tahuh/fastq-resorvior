/*
 * fastq-resorvior.c
 * compilation : gcc -O2 -o fastq-resorvior fastq-resorvior.c -lz
 * Author : Thomas Sunghoon Heo
 * Resoirvior sampling algorith from : https://en.wikipedia.org/wiki/Reservoir_sampling 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <getopt.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)
unsigned long long total_reads = 0;
unsigned long long n_sample = 0;
float factor = 1.0;
kseq_t *read1 = NULL;
kseq_t *read2 = NULL;
gzFile fp1 = NULL;
gzFile fp2 = NULL;
char *fastq1_name = NULL;
char *fastq2_name = NULL;
char *mode = NULL;

typedef struct s_entry{
	char *name;
	char *comment;
	char *seq;
	char *qual;
}t_entry;

t_entry** read1_sampled = NULL;
t_entry** read2_sampled = NULL;

const char *usage = "./fastq-resorvior -m mode [options] in1.fastq[,in2.fastq] out1.fastq[,out2.fastq]\n"
"\n"
"[options]\n"
"\n"
"   -h          This help message\n"
"   -f   FLT    Factor of sampling. If -k is also set, then the smaller value using \n"
"               this factor and -k result will be sampled [default=1.0]\n"
"   -k   INT    Actual number of reads to sample. If -f is also set, then the smaller \n"
"               values using -f and -k will be used [default=all_reads]\n"
"   -m   STR    Sequencing read mode either s or p to denote single end and paried end\n"
"               [default=s]\n"
"\n"
"Example command for single end read sampling\n\n"
"./fastq-resorvior in.fastq out.fastq\n\n"
"Example command for paired end read sampling\n\n"
"./fastq-resorvior in1.fastq in2.fastq out1.fastq out2.fastq\n\n";

// function forward declaration
void flush_result(FILE *ofstream, t_entry** space, unsigned long long n);
void destroy_sample_space(t_entry **space, unsigned long long n);
t_entry** init_sample_space(unsigned long long n);
void assign_or_update_at(char *name, char *comment, char *seq, char *qual, t_entry** bucket, int loc);
unsigned long long count_reads_se(const char *fname);
unsigned long long count_reads_pe(const char *fname1, const char *fname2);
int se_mode(const char *infile, const char *outfile);
int pe_mode(const char *infile1, const char *infile2, const char *outfile1, const char *outfile2);
int init_params(int argc, char **argv);
void clear_params();


int main(int argc, char **argv)
{
	int ind = init_params(argc, argv);
	int retcode = 0;
	if(strcmp(mode , "s") == 0){
		retcode = se_mode(argv[optind], argv[optind+1]);
	} else {
		retcode = pe_mode(argv[optind], argv[optind+1], argv[optind+2], argv[optind+3]);
	}
	clear_params();
	return retcode;
}

// Function body
void flush_result(FILE *ofstream, t_entry** space, unsigned long long n)
{
	unsigned long long i = 0;
	for(i = 0; i < n; i++){
		fprintf(ofstream, "@%s", space[i]->name);
		if(space[i]->comment != NULL){
			fprintf(ofstream, " %s", space[i]->comment);
		}
		fprintf(ofstream, "\n%s\n+\n%s\n", space[i]->seq, space[i]->qual);
	}
}
void destroy_sample_space(t_entry **space, unsigned long long n)
{
	unsigned long long i = 0;
	for(i = 0; i < n; i++){
		assign_or_update_at(NULL,NULL,NULL,NULL, space, i);
		free(space[i]);
	}
	free(space);
}
t_entry ** init_sample_space( unsigned long long n)
{
	t_entry **space = (t_entry **)malloc(sizeof(t_entry*) * n);
	unsigned long long i = 0;
	for(i = 0; i < n; i++){
		space[i] = (t_entry *)malloc(sizeof(t_entry));
		space[i]->name = NULL; space[i]->comment = NULL; space[i]->seq = NULL; space[i]->qual=NULL;
	}
	return space;
}
void assign_or_update_at(char *name, char *comment, char *seq, char *qual, t_entry** bucket, int loc)
{
	t_entry *bkt_ptr = bucket[loc];
	if(bkt_ptr->name != NULL){
		free(bkt_ptr->name); bkt_ptr->name = NULL;
	}
	if(bkt_ptr->comment != NULL){
		free(bkt_ptr->comment); bkt_ptr->comment=NULL;
	}
	if(bkt_ptr->seq != NULL){
		free(bkt_ptr->seq); bkt_ptr->seq = NULL;
	}
	if(bkt_ptr->qual != NULL){
		free(bkt_ptr->qual); bkt_ptr->qual = NULL;
	}
	if(name != NULL){
		bkt_ptr->name = strdup(name);
	}
	if(comment != NULL){
		bkt_ptr->comment =strdup(comment);
	}
	if(seq != NULL){
		bkt_ptr->seq = strdup(seq);
	}
	if(qual != NULL){
		bkt_ptr->qual = strdup(qual);
	}
}

unsigned long long count_reads_se(const char *fname) 
{
    unsigned long long n = 0;
	fp1 = gzopen(fname, "r");
	if(fp1 == NULL){
		fprintf(stderr, "[fastq-resorvior::count_reads_se] Error open file %s\n", fname);
		return 0;
	}
	read1 = kseq_init(fp1);
	if(read1 == NULL){
		fprintf(stderr, "[fastq-resorvior::count_reads_se] Error init kseq\n");
		gzclose(fp1);
		return 0;
	}
	int l;
	while((l=kseq_read(read1)) >= 0){
		n++;
	}
	kseq_destroy(read1);
	gzclose(fp1);
	return n;
}

unsigned long long count_reads_pe(const char *fname1, const char *fname2)
{
	unsigned long long n = 0;
	fp1 = gzopen(fname1, "r");
	if(fp1 == NULL){
		fprintf(stderr, "[fastq-resorvior::count_reads_se] Error open file %s\n", fname1);
		return 0;
	}
	read1 = kseq_init(fp1);
	if(read1 == NULL){
		fprintf(stderr, "[fastq-resorvior::count_reads_se] Error init kseq\n");
		gzclose(fp1);
		return 0;
	}
	fp2 = gzopen(fname2, "r");
	if(fp1 == NULL){
		fprintf(stderr, "[fastq-resorvior::count_reads_se] Error open file %s\n", fname2);
		gzclose(fp1); kseq_destroy(read1);
		return 0;
	}
	read2 = kseq_init(fp2);
	if(read2 == NULL){
		fprintf(stderr, "[fastq-resorvior::count_reads_se] Error init kseq\n");
		gzclose(fp1); kseq_destroy(read1); gzclose(fp2); 
		return 0;
	}
	int l1, l2;
	while( ((l1=kseq_read(read1)) >= 0) && ((l2 =kseq_read(read2))>=0)){
		n++;
	}
	kseq_destroy(read1);
	gzclose(fp1);
	kseq_destroy(read2);
	gzclose(fp2);
	return n;
}

unsigned long long calculate_read_samples_by_factors(unsigned long long n, float f)
{
	unsigned long long cnt = (unsigned long long)(n * f);
	return cnt;
}


int init_params(int argc, char **argv)
{
	if(argc <= 1){
		fprintf(stderr, "%s\n", usage);
		exit(-1);
	}
	int c;
	while((c=getopt(argc, argv, "hf:k:m:")) != -1){
		switch(c){
			case 'f': factor = atof(optarg); break;
			case 'm': mode = strdup(optarg); break;
			case 'k': n_sample = (unsigned long long)atoll(optarg); break;
		}
	}
	if(mode == NULL){
		mode = strdup("s");
	}
	return optind;
}

void clear_params()
{
	if(mode != NULL){
		free(mode);
	}
}
int se_mode(const char *infile, const char *outfile)
{
	fprintf(stderr, "[fastq-resorvior::se_mode] Counting total reads present in FASTQ file\n");
	total_reads = count_reads_se(infile);
	fprintf(stderr, "[fastq-resorvior::se_mode] Total %llu reads present in FASTQ file\n", total_reads);
	if(total_reads == 0){
		fprintf(stderr, "[fastq-resorvior] No reads present. Exit\n");
		return -1;
	}
	if(factor == 1.0){
		if(n_sample == 0){
			FILE *ofp = fopen(outfile, "w");
			fprintf(stderr, "[fastq-resorvior] Sample all reads...\n");
			gzFile tmpfp = gzopen(infile, "r");
			kseq_t *tmp = kseq_init(tmpfp);
			int l = 0;
			while((l=kseq_read(tmp)) >= 0){
				fprintf(ofp, "%s\n%s\n+%s\n", tmp->name.s, tmp->seq.s, tmp->qual.s);
			}
			kseq_destroy(tmp);
			gzclose(tmpfp);
			fclose(ofp);
			return 1;
		}
		else {
			FILE *ofp = fopen(outfile, "w");
			if(n_sample >= total_reads){
				fprintf(stderr, "[fastq-resorvior] Sample all reads...\n");
				gzFile tmpfp = gzopen(infile, "r");
				kseq_t *tmp = kseq_init(tmpfp);
				int l = 0;
				while((l=kseq_read(tmp)) >= 0){
					fprintf(ofp, "%s\n%s\n+%s\n", tmp->name.s, tmp->seq.s, tmp->qual.s);
				}
				kseq_destroy(tmp);
				gzclose(tmpfp);
				fclose(ofp);
				return 1;
			} else {
				read1_sampled = init_sample_space(n_sample);
				fprintf(stderr, "[fastq-resorvior] Sample %llu reads...\n", n_sample);
				gzFile tmpfp = gzopen(infile, "r");
				kseq_t *tmp = kseq_init(tmpfp);
				int l = 0;
				int i = 0;
				while((l=kseq_read(tmp)) >= 0){
					if(i < n_sample){
						assign_or_update_at(tmp->name.s, tmp->comment.s, tmp->seq.s, tmp->qual.s, read1_sampled, i);
					} else {
						int j = rand() % i;
						if( j < n_sample) {
							assign_or_update_at(tmp->name.s, tmp->comment.s, tmp->seq.s, tmp->qual.s, read1_sampled, j);
						}
					}
					
				}
				i++;
				kseq_destroy(tmp);
				gzclose(tmpfp);
				flush_result(ofp, read1_sampled, n_sample);
				fclose(ofp);
				destroy_sample_space(read1_sampled,n_sample);
				return 1;
			}
		}
	} else {
		unsigned long long new_n = (unsigned long long)(total_reads * factor);
		if(n_sample != 0){
			if(new_n <=  n_sample){
				n_sample = new_n;
			}
		} else {
			n_sample = new_n;
		}
		FILE *ofp = fopen(outfile, "w");
		read1_sampled = init_sample_space( n_sample);
		fprintf(stderr, "[fastq-resorvior] Sample %llu reads...\n", n_sample);
		gzFile tmpfp = gzopen(infile, "r");
		kseq_t *tmp = kseq_init(tmpfp);
		int l = 0;
		int i = 0;
		while((l=kseq_read(tmp)) >= 0){
			if(i < n_sample){
				assign_or_update_at(tmp->name.s, tmp->comment.s, tmp->seq.s, tmp->qual.s, read1_sampled, i);
			} else {
				int j = rand() % i;
				if( j < n_sample) {
					assign_or_update_at(tmp->name.s, tmp->comment.s, tmp->seq.s, tmp->qual.s, read1_sampled, j);
				}
			}
			i++;
			
		}
		kseq_destroy(tmp);
		gzclose(tmpfp);
		flush_result(ofp, read1_sampled, n_sample);
		fclose(ofp);
		destroy_sample_space(read1_sampled,n_sample);
		return 1;
	}
    
}


int pe_mode(const char *infile1, const char *infile2, const char *outfile1, const char *outfile2)
{
	fprintf(stderr, "[fastq-resorvior::pe_mode] Counting total reads present in FASTQ file\n");
	total_reads = count_reads_pe(infile1, infile2);
	fprintf(stderr, "[fastq-resorvior::pe_mode] Total %llu reads present in FASTQ file\n", total_reads);
	if(total_reads == 0){
		fprintf(stderr, "[fastq-resorvior] No reads present. Exit\n");
		return -1;
	}
	if(factor == 1.0){
		if(n_sample == 0){
			FILE *ofp = fopen(outfile1, "w");
			FILE *ofp2 = fopen(outfile2, "w");
			fprintf(stderr, "[fastq-resorvior] Sample all reads...\n");
			gzFile tmpfp = gzopen(infile1, "r");
			gzFile tmpfp2 = gzopen(infile2, "r");
			kseq_t *tmp = kseq_init(tmpfp);
			kseq_t *tmp2 = kseq_init(tmpfp2);
			int l = 0;
			int l2 = 0;
			while(((l=kseq_read(tmp)) >= 0)  && ((l2=kseq_read(tmp2)) >= 0) ){
				fprintf(ofp, "%s\n%s\n+%s\n", tmp->name.s, tmp->seq.s, tmp->qual.s);
				fprintf(ofp2, "%s\n%s\n+%s\n", tmp2->name.s, tmp2->seq.s, tmp2->qual.s);
			}
			kseq_destroy(tmp);
			gzclose(tmpfp);
			fclose(ofp);
			kseq_destroy(tmp2);
			gzclose(tmpfp2);
			fclose(ofp2);
			return 1;
		}
		else {
			if(n_sample >= total_reads){
				FILE *ofp = fopen(outfile1, "w");
				FILE *ofp2 = fopen(outfile2, "w");
				fprintf(stderr, "[fastq-resorvior] Sample all reads...\n");
				gzFile tmpfp = gzopen(infile1, "r");
				gzFile tmpfp2 = gzopen(infile2, "r");
				kseq_t *tmp = kseq_init(tmpfp);
				kseq_t *tmp2 = kseq_init(tmpfp2);
				int l = 0;
				int l2 = 0;
				while(((l=kseq_read(tmp)) >= 0)  && ((l2=kseq_read(tmp2)) >= 0) ){
					fprintf(ofp, "%s\n%s\n+%s\n", tmp->name.s, tmp->seq.s, tmp->qual.s);
					fprintf(ofp2, "%s\n%s\n+%s\n", tmp2->name.s, tmp2->seq.s, tmp2->qual.s);
				}
				kseq_destroy(tmp);
				gzclose(tmpfp);
				fclose(ofp);
				kseq_destroy(tmp2);
				gzclose(tmpfp2);
				fclose(ofp2);
				return 1;
			} else {
				read1_sampled = init_sample_space(n_sample);
				read2_sampled = init_sample_space(n_sample);
				fprintf(stderr, "[fastq-resorvior] Sample %llu reads...\n", n_sample);
				gzFile tmpfp = gzopen(infile1, "r");
				gzFile tmpfp2 = gzopen(infile2, "r");
				kseq_t *tmp = kseq_init(tmpfp);
				kseq_t *tmp2 = kseq_init(tmpfp2);
				FILE *ofp = fopen(outfile1, "w");
				FILE *ofp2 = fopen(outfile2, "w");
				int l = 0; int l2 = 0;
				int i = 0;
				while(((l=kseq_read(tmp)) >= 0)  && ((l2=kseq_read(tmp2)) >= 0) ){
					if(i<n_sample){
						assign_or_update_at(tmp->name.s, tmp->comment.s, tmp->seq.s, tmp->qual.s, read1_sampled, i);
						assign_or_update_at(tmp2->name.s, tmp2->comment.s, tmp2->seq.s, tmp2->qual.s, read2_sampled, i);
					} else {
						int j = rand() % i;
						if( j < n_sample) {
							assign_or_update_at(tmp->name.s, tmp->comment.s, tmp->seq.s, tmp->qual.s, read1_sampled, j);
							assign_or_update_at(tmp2->name.s, tmp2->comment.s, tmp2->seq.s, tmp2->qual.s, read2_sampled, j);
						}
					}
					i++;
				}
				
				kseq_destroy(tmp);
				kseq_destroy(tmp2);
				gzclose(tmpfp);
				gzclose(tmpfp2);
				flush_result(ofp, read1_sampled, n_sample);
				flush_result(ofp2, read2_sampled, n_sample);
				fclose(ofp);fclose(ofp2);
				destroy_sample_space(read1_sampled,n_sample);
				destroy_sample_space(read2_sampled,n_sample);
				return 1;
			}
		}
	} else {
		unsigned long long new_n = (unsigned long long)(total_reads * factor);
		if(n_sample != 0){
			if(new_n <=  n_sample){
				n_sample = new_n;
			}
		} else {
			n_sample = new_n;
		}
		read1_sampled=init_sample_space(n_sample);
		read2_sampled=init_sample_space(n_sample);
		fprintf(stderr, "[fastq-resorvior] Sample %llu reads...\n", n_sample);
		gzFile tmpfp = gzopen(infile1, "r");
		gzFile tmpfp2 = gzopen(infile2, "r");
		kseq_t *tmp = kseq_init(tmpfp);
		kseq_t *tmp2 = kseq_init(tmpfp2);
		FILE *ofp = fopen(outfile1, "w");
		FILE *ofp2 = fopen(outfile2, "w");
		int l = 0; int l2 = 0;
		int i = 0;
		while(((l=kseq_read(tmp)) >= 0)  && ((l2=kseq_read(tmp2)) >= 0) ){
			if(i<n_sample){
				assign_or_update_at(tmp->name.s, tmp->comment.s, tmp->seq.s, tmp->qual.s, read1_sampled, i);
				assign_or_update_at(tmp2->name.s, tmp2->comment.s, tmp2->seq.s, tmp2->qual.s, read2_sampled, i);
			} else {
				int j = rand() % i;
				if( j < n_sample) {
					assign_or_update_at(tmp->name.s, tmp->comment.s, tmp->seq.s, tmp->qual.s, read1_sampled, j);
					assign_or_update_at(tmp2->name.s, tmp2->comment.s, tmp2->seq.s, tmp2->qual.s, read2_sampled, j);
				}
			}
			i++;
		}
		
		kseq_destroy(tmp);
		kseq_destroy(tmp2);
		gzclose(tmpfp);
		gzclose(tmpfp2);
		flush_result(ofp, read1_sampled, n_sample);
		flush_result(ofp2, read2_sampled, n_sample);
		fclose(ofp);fclose(ofp2);
		destroy_sample_space(read1_sampled,n_sample);
		destroy_sample_space(read2_sampled,n_sample);
		fprintf(stdout, "[fastq-resorvior::pe_mode] Done\n");
		return 1;
	}
}
