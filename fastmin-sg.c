/*
###############################################################################
# Author: Alex Di Genova 
# Laboratory: ERABLE/INRIA
# Copyright (c)
# year: 2019
###############################################################################
*/

#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include <zlib.h>
#include <pthread.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include "minimap.h"
#include "kseq.h"
#include <time.h>
#include "math.h"

KSEQ_INIT(gzFile, gzread)


//revcomp
//global array for  to compute revcomp
char comp_tab[] = {
        0,   1,       2,       3,       4,   5,       6,       7,       8,   9,  10,  11,      12,  13,  14,  15,
        16,  17,  18,  19,      20,  21,  22,  23,      24,  25,  26,  27,      28,  29,  30,  31,
        32,  33,  34,  35,      36,  37,  38,  39,      40,  41,  42,  43,      44,  45,  46,  47,
        48,  49,  50,  51,      52,  53,  54,  55,      56,  57,  58,  59,      60,  61,  62,  63,
        64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
        'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,      92,  93,  94,  95,
        64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
        'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};


//split_str
char** str_split(char* a_str, const char a_delim)
{
    char** result    = 0;
    size_t count     = 0;
    char* tmp        = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    /* Count how many elements will be extracted. */
    while (*tmp)
    {
        if (a_delim == *tmp)
        {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }

    /* Add space for trailing token. */
    count += last_comma < (a_str + strlen(a_str) - 1);

    /* Add space for terminating null string so caller
       knows where the list of returned strings ends. */
    count++;

    result = malloc(sizeof(char*) * count);

    if (result)
    {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);

        while (token)
        {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }

    return result;
}

//heng li revcomp this is for printing only
void revcomp(char* kmer, int len) {
    int c0, c1;
    int klen=len;
    for (int i = 0; i < klen>>1; ++i) { // reverse complement sequence
        c0 = comp_tab[(int)kmer[i]];
        c1 = comp_tab[(int)kmer[klen - 1 - i]];
        kmer[i] = c1;
        kmer[klen - 1 - i] = c0;
    }
    if (klen & 1) // complement the remaining base
        kmer[klen>>1] = comp_tab[(int)kmer[klen>>1]];
    //return kmer;
}

//funtions for rename long reads
static void stk_printstr(const kstring_t *s, unsigned line_len,FILE* outlongr)
{
    if (line_len != UINT_MAX && line_len != 0) {
        int i, rest = s->l;
        for (i = 0; i < s->l; i += line_len, rest -= line_len) {
            putc('\n',outlongr);
            if (rest > line_len) fwrite(s->s + i, 1, line_len, outlongr);
            else fwrite(s->s + i, 1, rest, outlongr);
        }
        putc('\n',outlongr);
    } else {
        putc('\n',outlongr);
        fputs(s->s,outlongr);
    }
}

static inline void stk_printseq_renamed(const kseq_t *s, int line_len, const char *prefix, int64_t n, FILE* outlongr)
{
    putc('>',outlongr);
    if (n >= 0) {
        if (prefix) fputs(prefix, outlongr);
        fprintf(outlongr,"%lld", (long long)n);
    } else fputs(s->name.s, outlongr);

    stk_printstr(&s->seq, line_len,outlongr);
}

#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))

typedef struct{
    int std;
    int avg;
} pe_stats;


//global variables
int inserts[50];
//total long thread_seqs
int64_t total_long_reads;
//global mutexs
pthread_mutex_t gseq;//for reading sequence
pthread_mutex_t pres;//for printing results
//struct for mapping variables
typedef struct {
    int read_len;
    int minq;
    int kmer_size;
    int minimizer_d;
    int moving_w;
    int qcov;
    int numseq2read;//long reads per thread
    int min_len_long_read;//minimal long read length
    int n_threads;
    int pair_rescue;//flag for pair rescue
} aligment_parameters;

//strcuture for threads
typedef struct {
    //pointer to read the seq from file
    kseq_t* ks;//forward read
    kseq_t* kr;//reverse reads
    //minimap2 index
    mm_idx_t *mi;//pointer to minimap2 database
    mm_mapopt_t mopt;//parameters for db construction in minimap
    mm_tbuf_t *tbuf;    //buffer for minimap2
    //output files
    int libcount;//number of synthetic libs to create
    FILE** outFiles;//synthetic mate pair libraries
    FILE* outlong;//renamed long read file
    //mapping variables
    aligment_parameters *alnp;//pointer to aligment parameters
    pe_stats *pes;
    /*int read_len,minq,kmer_size,moving_w,qcov;
    int numseq2read;//long reads per thread
    int min_len_long_read;//minimal long read length*/
} thread_data_struct;

//structure that save temporary a seq in the local thread
typedef struct{
    char *name;
    char *seq;
    int length;
    int64_t id;
} thread_seqs;

//we need ks->name.s,d,l,fo,mi->seq[fwdh->rid].name,fwdh->rs,fwdh->mapq,kmer_size,fwdh->mlen,fwdh->blen
//dinamical list that saves the hits founds
typedef struct phit{
    //char* lrid;//long read id
    int64_t lrid;//long read unique id
    char* ctgid;//target name
    int32_t plr;
    int32_t orientation;//match orientation
    int32_t pos;//position in target
    int32_t mapq;//mapping quality
    int32_t mlen, blen; // seeded exact match length; seeded alignment block length
    int32_t qcov;
    //kmer seq
    char* k_seq;
    struct phit *next; //next hit result
} phit_t;

//we create a list to handle the hits
void push_hit(phit_t **head, int64_t lrid, char* ctgid,int32_t plr,int32_t o, int32_t pos, int32_t mapq, int32_t mlen, int32_t blen, char* k_seq, int qcov) {
    phit_t * new_node;
    new_node = malloc(sizeof(phit_t));
    //new_node->val = val;
    new_node->plr=plr;
    new_node->orientation=o;
    new_node->pos=pos;
    new_node->mapq=mapq;
    new_node->mlen=mlen;
    new_node->blen=blen;
    new_node->qcov=qcov;
    //we allocate the strings and copy them
    //new_node->lrid=malloc(sizeof(char) * (strlen(lrid)));
    //strcpy(new_node->lrid,lrid);
    new_node->lrid=lrid;
    new_node->ctgid=malloc(sizeof(char*) * (strlen(ctgid)));
    strcpy(new_node->ctgid,ctgid);
    new_node->k_seq=malloc(sizeof(char*) * (strlen(k_seq)));
    strcpy(new_node->k_seq,k_seq);
    //we update the header of the list
    new_node->next = *head;
    *head = new_node;
}

//thread mapping of long reads
void* maplongreads(void *arg)
{
  //we get the given variables
  //thread_data_struct *dargs = arg;
  thread_data_struct *dargs = (thread_data_struct*)arg;
  kseq_t *lks =dargs->ks;
  aligment_parameters *alnp = (aligment_parameters*)dargs->alnp;
  thread_seqs *localseqs=malloc(sizeof(thread_seqs) * (alnp->numseq2read));
  FILE** outFiles = dargs->outFiles;
  FILE* outlongr=dargs->outlong;
  //we set the mapping variables
  int read_len=alnp->read_len;
  int minq=alnp->minq;
  int kmer_size=alnp->kmer_size;
  int moving_w=alnp->moving_w;
  int qcov=alnp->qcov;//query coverage

  mm_idx_t *mi=dargs->mi;
  mm_mapopt_t mopt = dargs->mopt;
  mm_tbuf_t *tbuf=dargs->tbuf;
  int libcount=dargs->libcount;
  //we create a list of hits for each library to store the results
  phit_t** headhits = malloc(sizeof(phit_t*) * (libcount) );
  for (int i = 0; i < libcount; ++i) headhits[i]=NULL;

  //we read the sequence from the reads
  int rseq=0;
  //we create the fwd rev containers and the kmer_seq
  char* fwd=malloc(sizeof(char*) * (read_len));
  char* rev=malloc(sizeof(char*) * (read_len));
  char* kfwdseq=malloc(sizeof(char*) * (kmer_size));
  int fcov=0,rcov=0;
  //I set the values to 0 in the memory
  memset(fwd,0,read_len);
  memset(fwd,0,read_len);
  memset(fwd,0,kmer_size);
  //pointer to regions
  mm_reg1_t *regfwd; //fwd
  mm_reg1_t *regrev; //rev
  int n_regfwd=0; // number of regions to consider fwd
  int n_regrev=0; // number of regions to consider rev
  //pointers to hits
  mm_reg1_t *fwdh;
  mm_reg1_t *revh;
  //local threads variables
  int64_t total_pairs=0;
  int64_t aling_pairs=0;
  int64_t number_long_reads=0;
  int getjob=1;//flag to control the execution of threads jobs

  while(getjob){
      pthread_mutex_lock(&gseq);
      rseq=0;
      //the order is always important
      while(rseq < alnp->numseq2read-1 && kseq_read(lks) >=0 ){
      //we copy the sequence in
      if(lks->seq.l >= alnp->min_len_long_read) {
          stk_printseq_renamed(lks, 60, 0, total_long_reads, outlongr);
          localseqs[rseq].name = malloc(sizeof(char *) * (lks->name.l));
          localseqs[rseq].id = total_long_reads;
          memset(localseqs[rseq].name, 0, lks->name.l);
          memcpy(localseqs[rseq].name, lks->name.s, lks->name.l);
          localseqs[rseq].name[lks->name.l] = '\0';
          localseqs[rseq].length = lks->seq.l;
          localseqs[rseq].seq = malloc(sizeof(char *) * (lks->seq.l));
          memset(localseqs[rseq].seq, 0, lks->seq.l);
          memcpy(localseqs[rseq].seq, lks->seq.s, lks->seq.l);
          localseqs[rseq].seq[lks->seq.l] = '\0';
          //we check that the sequence leenght is equal to the sequence
          assert(strlen(localseqs[rseq].seq) == lks->seq.l);
          rseq++;
      }
          total_long_reads++;
     }
    pthread_mutex_unlock(&gseq);
    if(rseq){
   //we map the reads;
    for(int i=0; i< rseq ; i++){ // each kseq_read() call reads one query sequence
      number_long_reads++;
     for(int ii=0; ii<libcount; ii++){
      int d=inserts[ii];
      int max_d = localseqs[i].length - d;
      if(max_d < 0){
        continue;
      }
      //create the read of the current library
      for(int l=0; l<max_d; l+=moving_w){
        //we clean the statics containers
        memset(fwd,0,read_len);
        memset(rev,0,read_len);
        memset(kfwdseq,0,kmer_size);
          int pa=l+1;
           memcpy(fwd, localseqs[i].seq+pa, read_len);
           fwd[read_len]='\0';
            int pb=l+d-read_len;
            memcpy(rev, localseqs[i].seq+pb, read_len);
            rev[read_len]='\0';
            revcomp(rev,read_len);
            rev[read_len]='\0';
          //we check that the length of fwd, rev is equal to the read_len
          assert(strlen(rev)==strlen(fwd) && strlen(fwd)==read_len);
      n_regfwd=0;
      n_regrev=0;
      total_pairs++;
      //we actually map the reads using minimpa2 API
      regfwd = mm_map(mi, read_len, fwd, &n_regfwd, tbuf, &mopt, 0); // get all hits for the query
      if(n_regfwd > 0){
        regrev = mm_map(mi, read_len, rev, &n_regrev, tbuf, &mopt, 0); // get all hits for the query
        if(n_regrev > 0){
            if(n_regrev == 1 && n_regfwd == 1){
            //we print the aligment results
              fwdh=&regfwd[0];
              revh=&regrev[0];
              //coverage
              fcov=(int)(((float)fwdh->blen/(float)read_len) * 100);
              rcov=(int)(((float)revh->blen/(float)read_len) * 100);
              //printf("%d %d\n",fcov,rcov );
              //we check that the aligment cover the query
                if(fwdh->mapq >= minq && revh->mapq >= minq && fcov >=qcov && rcov>=qcov ){
                  //we add the reverse first to the hits pool
                  strncpy(kfwdseq, rev, kmer_size);
                  kfwdseq[kmer_size]='\0';
                  int fr=0;
                  if(revh->rev){
                    fr=16;
                  }
		//printf("%s %s\n",mi->seq[revh->rid].name,mi->seq[fwdh->rid].name);
                  push_hit(&headhits[ii],localseqs[i].id,mi->seq[revh->rid].name,l,fr,fr > 0 ? revh->re:revh->rs,revh->mapq,revh->mlen,revh->blen,kfwdseq,rcov);
                  //we add the fwd
                    int fo=0;
                    if(fwdh->rev){
                      fo=16;
                    }
                    strncpy(kfwdseq, fwd, kmer_size);
                    kfwdseq[kmer_size]='\0';
                    push_hit(&headhits[ii],localseqs[i].id,mi->seq[fwdh->rid].name,l,fo,fo > 0 ? fwdh->re:fwdh->rs,fwdh->mapq,fwdh->mlen,fwdh->blen,kfwdseq,fcov);
                  aling_pairs++;
                }
              }
          }
        }
      //we clean the hits
      if(n_regfwd > 0)
        free(regfwd);
      if(n_regrev > 0)
        free(regrev);
     }

    }
  }
  //we printnt the results
    pthread_mutex_lock(&pres);
    //we store the hits in the files
    phit_t* nh = NULL;//next hit
    for (int f = 0; f < libcount; ++f)
    {
      //this way print the reverse before the foward
      //the header is not null
      nh=headhits[f];
      while(nh != NULL){
           //fprintf(outFiles[f],"%lld_FG_%d_%d\t%d\t%s\t%d\t%d\t%dM\t*\t%d\t%d\t%s\t%d\n",
             //      (long long)nh->lrid,inserts[f],nh->plr,nh->orientation,nh->ctgid,nh->pos,nh->mapq,kmer_size,nh->mlen,nh->blen,nh->k_seq,nh->qcov);
                     //read_id,orientation,ctg_id,pos,mapq,ksize,qcov     	
                    //fprintf(outshort,"%lld\t%d\t%s\t%d\t%d\t%d\t%d\n",
                     //       (long long)nh->lrid,nh->orientation,nh->ctgid,nh->pos,nh->mapq,kmer_size,nh->qcov);
           fprintf(outFiles[f],"%lld_FG_%d_%d\t%d\t%s\t%d\t%d\t%d\t%d\n",
                   (long long)nh->lrid,inserts[f],nh->plr,nh->orientation,nh->ctgid,nh->pos,nh->mapq,kmer_size,nh->qcov);
          headhits[f]=headhits[f]->next;
          //we free the allocated hits
          //free(nh->lrid);
          free(nh->ctgid);
          free(nh->k_seq);
          free(nh);
          nh=headhits[f];
      }
    }

     //if(number_long_reads%4995==0)
        fprintf(stderr,"Thread_info:Partial LONGREADS=%lld MAPPED_PAIRS=%lld TOTAL_PAIRS=%lld RL=%d MINQ=%d KS=%d MW=%d\n",
                (long long)number_long_reads,(long long)aling_pairs,(long long)total_pairs,read_len,minq,kmer_size,moving_w);

    //we clean the container of reads before allocate more memory in a single threads fashion
    for (size_t i = 0; i < rseq; i++) {
      free(localseqs[i].seq);
      free(localseqs[i].name);
    }
    pthread_mutex_unlock(&pres);

}else{
  getjob=0;
 }

}
//we print some thread information
pthread_mutex_lock(&pres);
fprintf(stderr,"Final_Thread_info:Total LONGREADS=%lld MAPPED_PAIRS=%lld TOTAL_PAIRS=%lld RL=%d MINQ=%d KS=%d MW=%d\n",
        (long long)number_long_reads,(long long)aling_pairs,(long long)total_pairs,read_len,minq,kmer_size,moving_w);
pthread_mutex_unlock(&pres);
//all done in this thread
return NULL;
}


int map_long_read_files(char* linserts, char* queryfile, char* prefix,  aligment_parameters *alnp, mm_idx_t *mi, mm_mapopt_t mopt){
    //we populate the array of libraries
    char* token=strtok(linserts,",");
    int libcount=0;
    while(token != NULL)
    {
        //printf("IS=%d\n",atoi(token));
        inserts[libcount]=atoi(token);
        token=strtok(NULL,",");
        libcount++;
        if(libcount >= 50){
            fprintf(stderr,"Max number of synthetic libraries is 50\n");
            return 1;
        }
    }

    printf("Query file  = %s\n",queryfile);
    // open query file for reading; you may use your favorite FASTA/Q parser
    gzFile f = gzopen(queryfile, "r");
    //printf("%d",f);
    if (! f) {
        fprintf (stderr, "gzopen of '%s' failed: %s.\n", queryfile,
                 strerror (errno));
        exit (EXIT_FAILURE);
    }
    //assert(f);
    kseq_t *ks = kseq_init(f);
    FILE** outFiles = malloc(sizeof(FILE*) * (libcount));
    char filename[100];//max length of file name
    memset(filename, '\0', sizeof(filename));
    sprintf(filename, "longreads.%s.fa",prefix);
    FILE* outlong = fopen(filename,"w");

    for(int i=0; i<libcount; i++){
        //sprintf(filename, "lib.I%d.%s.fastmin-sg.sam", inserts[i],prefix);
        sprintf(filename, "%s.I%d.fm.sam", prefix,inserts[i]);
        printf("Lib : %s\n",filename);
        outFiles[i] = fopen(filename,"w");
    }
    printf("A total of %d syntethic libraries will be created by FastMin-SG\n",libcount);
//we start the mapping of the reads

//we create the threads
    if (pthread_mutex_init(&gseq, NULL) != 0){
        printf("\n mutex for file reading failed\n");
        return 1;
    }
    if (pthread_mutex_init(&pres, NULL) != 0){
        printf("\n mutex for printing  failed\n");
        return 1;
    }
    //we start the mapping
    pthread_t *tid=(pthread_t *)malloc(alnp->n_threads * sizeof(pthread_t ));

    mm_tbuf_t **tbuf;// = mm_tbuf_init(n_threads); // thread buffer; for multi-threading, allocate one tbuf for each thread
    tbuf = (mm_tbuf_t**)calloc(alnp->n_threads, sizeof(mm_tbuf_t*));
    for(int i =0; i < alnp->n_threads; ++i) tbuf[i] = mm_tbuf_init();

    for(int i=0; i < alnp->n_threads; i++){
        thread_data_struct *data = malloc(sizeof *data);
        data->ks=ks;
        data->outFiles=outFiles;
        data->outlong=outlong;
        data->mi=mi;
        data->libcount=libcount;
        data->mopt=mopt;
        data->tbuf=tbuf[i];
        //mapping variables
        data->alnp=alnp;
        int err = pthread_create(&tid[i], NULL, &maplongreads, data);
        if (err != 0){
            printf("\ncan't create thread :[%s]", strerror(err));
        }
    }
    //we wait for the thread to finish
    for(int i=0; i< alnp->n_threads; ++i) pthread_join(tid[i], NULL);
    //we destroy the thhread buffers
    for(int i=0; i<alnp->n_threads; i++) mm_tbuf_destroy(tbuf[i]);

    //end of threading code
    kseq_destroy(ks); // close the query file
    gzclose(f);//flose gzfile
    //close sinthetic library files and fasta file of long reads
    for(int i=0; i<libcount ; i++){
        fclose(outFiles[i]);
    }
    // todo: be a variable file ?
    fclose(outlong);
    return 0;
}



int rescue_pair(mm_reg1_t *r1,mm_reg1_t *ra,int n_reg,pe_stats *ilib){
    //we check return the pair satisfiying the first distance, otherwise we return a negative number indicating a null result
    int d1= r1->rev ? r1->re:r1->rs;
    int min=ilib->avg-(int)2.5*ilib->std;
    int max=ilib->avg+(int)2.5*ilib->std;

    mm_reg1_t *r2;//region that iter for looking for a candidate to rescue
    for(int i=0; i<n_reg; i++){
        r2=&ra[i];
        int d2=r2->rev ? r2->re:r2->rs;
        //same contig
        if(r1->rid == r2->rid){
            int d=abs(d2-d1);
            //expected insert size
            if(d>=min && d<=max){
                //printf("RR:%d %d %d %d %d %d\n",min, d,max,r1->rid,r1->mapq,r2->mapq);
                return i;
            }
        }
    }
    return -1;
}

//thread mapping of short reads
void* mapshortreads(void *arg) {
    //we get the given variables
    //thread_data_struct *dargs = arg;
    thread_data_struct *dargs = (thread_data_struct*)arg;
    kseq_t *lks =dargs->ks;//fwd
    kseq_t *lkr =dargs->kr;//rev
    aligment_parameters *alnp = (aligment_parameters*)dargs->alnp;
    pe_stats *ilib=(pe_stats*)dargs->pes;
    //fprintf(stderr,"Average insert size %d std %d [-10 outlayers,]\n",ilib->avg, ilib->std);
    //fwd pool
    thread_seqs *localseqs=malloc(sizeof(thread_seqs) * (alnp->numseq2read));
    //rev pool
    thread_seqs *localseqr=malloc(sizeof(thread_seqs) * (alnp->numseq2read));
    //FILE** outFiles = dargs->outFiles;
    FILE* outshort=dargs->outlong;
    //we set the mapping variables
    int read_len=0;//todo make it a variable
    int minq=alnp->minq;
    int kmer_size=alnp->kmer_size;
    int moving_w=alnp->moving_w;
    int qcov=alnp->qcov;//query coverage

    mm_idx_t *mi=dargs->mi;
    mm_mapopt_t mopt = dargs->mopt;
    mm_tbuf_t *tbuf=dargs->tbuf;
    //we have oonly one library with short reads
    int libcount=1;

    //we create a list of hits for each library to store the results
    phit_t** headhits = malloc(sizeof(phit_t*) * (libcount) );
    for (int i = 0; i < libcount; ++i) headhits[i]=NULL;

    //we read the sequence from the reads
    int rseq=0;
    //we create the fwd rev containers and the kmer_seq
    //char* fwd=malloc(sizeof(char*) * (read_len));
   // char* rev=malloc(sizeof(char*) * (read_len));
    char* kfwdseq=malloc(sizeof(char*) * (kmer_size));
    int fcov=0,rcov=0;
    //I set the values to 0 in the memory
    //memset(fwd,0,read_len);
    //memset(fwd,0,read_len);
    memset(kfwdseq,0,kmer_size);
    //pointer to regions
    mm_reg1_t *regfwd; //fwd
    mm_reg1_t *regrev; //rev
    int n_regfwd=0; // number of regions to consider fwd
    int n_regrev=0; // number of regions to consider rev
    int64_t rescue=0;
    int64_t prescue=0;
    int pair_rescue=alnp->pair_rescue;//flag to set on off pair rescue;
    //pointers to hits
    mm_reg1_t *fwdh;
    mm_reg1_t *revh;
    //local threads variables
    //int64_t total_pairs=0;
    int64_t aling_pairs=0;
    int64_t number_read_pairs=0;
    //read_len=localseqr[rseq].length;
    int getjob=1;//flag to control the execution of threads jobs

    while(getjob){
        pthread_mutex_lock(&gseq);
        rseq=0;
        //the order is always important
        while(rseq < alnp->numseq2read-1 && kseq_read(lks) >=0 && kseq_read(lkr) >=0 ){
                if(read_len == 0){
                    read_len=lks->seq.l;
                }
                //we copy the sequence in fwd
                localseqs[rseq].name = malloc(sizeof(char *) * (lks->name.l));
                localseqs[rseq].id = total_long_reads;
                memset(localseqs[rseq].name, 0, lks->name.l);
                memcpy(localseqs[rseq].name, lks->name.s, lks->name.l);
                localseqs[rseq].name[lks->name.l] = '\0';
                localseqs[rseq].length = lks->seq.l;
                localseqs[rseq].seq = malloc(sizeof(char *) * (lks->seq.l));
                memset(localseqs[rseq].seq, 0, lks->seq.l);
                memcpy(localseqs[rseq].seq, lks->seq.s, lks->seq.l);
                localseqs[rseq].seq[lks->seq.l] = '\0';
                //we check that the sequence leenght is equal to the sequence
                assert(strlen(localseqs[rseq].seq) == lks->seq.l);
                //we copy the reverse seq
                localseqr[rseq].name = malloc(sizeof(char *) * (lkr->name.l));
                localseqr[rseq].id = total_long_reads;
                memset(localseqr[rseq].name, 0, lkr->name.l);
                memcpy(localseqr[rseq].name, lkr->name.s, lkr->name.l);
                localseqr[rseq].name[lkr->name.l] = '\0';
                localseqr[rseq].length = lkr->seq.l;
                localseqr[rseq].seq = malloc(sizeof(char *) * (lkr->seq.l));
                memset(localseqr[rseq].seq, 0, lkr->seq.l);
                memcpy(localseqr[rseq].seq, lkr->seq.s, lkr->seq.l);
                localseqr[rseq].seq[lkr->seq.l] = '\0';
                //we check that the sequence leenght is equal to the sequence
                assert(strlen(localseqr[rseq].seq) == lkr->seq.l);
                //printf("%s %s\n",localseqs[rseq].name,localseqr[rseq].name);
            rseq++;
            total_long_reads++;
        }
        pthread_mutex_unlock(&gseq);

        if(rseq){
            //we map the reads;
            for(int i=0; i< rseq ; i++){ // each kseq_read() call reads one query sequence
                number_read_pairs++;
                //we actually map the reads using minimpa2 API
                n_regfwd=0,n_regrev=0;
                        regfwd = mm_map(mi, localseqs[i].length, localseqs[i].seq, &n_regfwd, tbuf, &mopt, 0); // get all hits for the query
                        if(n_regfwd > 0){
                            regrev = mm_map(mi, localseqr[i].length, localseqr[i].seq, &n_regrev, tbuf, &mopt, 0); // get all hits for the query
                            if(n_regrev > 0){
                                if(n_regrev == 1 && n_regfwd == 1){
                                    //we print the aligment results
                                    fwdh=&regfwd[0];
                                    revh=&regrev[0];
                                    //coverage
                                    fcov=(int)(((float)fwdh->blen/(float)localseqs[i].length) * 100);
                                    rcov=(int)(((float)revh->blen/(float)localseqr[i].length) * 100);
                                    //printf("%d %d\n",fcov,rcov );
                                    //we check that the aligment cover the query
                                    if(fwdh->mapq >= minq && revh->mapq >= minq && fcov >=qcov && rcov>=qcov && fcov <=110 && rcov<=110){
                                        //we add the reverse first to the hits pool
                                        strncpy(kfwdseq, localseqr[i].seq, kmer_size);
                                        kfwdseq[kmer_size]='\0';
                                        int fr=0;
                                        if(revh->rev){
                                            fr=16;
                                        }
                                        push_hit(&headhits[0],localseqr[i].id,mi->seq[revh->rid].name,2,fr,fr > 0 ? revh->re:revh->rs,revh->mapq,revh->mlen,revh->blen,kfwdseq,rcov);

                                        //we add the fwd
                                        int fo=0;
                                        if(fwdh->rev){
                                            fo=16;
                                        }
                                        strncpy(kfwdseq, localseqs[i].seq, kmer_size);
                                        kfwdseq[kmer_size]='\0';
                                        push_hit(&headhits[0],localseqs[i].id,mi->seq[fwdh->rid].name,1,fo,fo > 0 ? fwdh->re:fwdh->rs,fwdh->mapq,fwdh->mlen,fwdh->blen,kfwdseq,fcov);

                                        //printf("%lld %s %s %d %d %d %d %d %d %d %d %d %d %d %d\n",localseqs[i].id,mi->seq[fwdh->rid].name,mi->seq[revh->rid].name,
                                          //     fwdh->rs,fwdh->re,revh->rs,revh->re,fo,fr,fwdh->qs,fwdh->qe,revh->qs,revh->qe, fo > 0 ? fwdh->re:fwdh->rs, fr > 0 ? revh->re:revh->rs);
                                        aling_pairs++;
                                    }
                                }else{
                                    //we can try rescue the reads
                                    //printf("%d %d\n",n_regfwd,n_regrev);
                                    if((n_regfwd == 1 || n_regrev == 1)  && pair_rescue == 1){
                                        prescue++; //we can rescue 10% of read pairs, focus within a contig rescue only, not between contigs
                                        int p=-1;
                                        //we try to rescue the pair within a contig
                                        //we fix fwd
                                        if(n_regfwd == 1){
                                             p=rescue_pair(&regfwd[0],regrev,n_regrev,ilib);
                                            fwdh=&regfwd[0];
                                            if(p>=0)
                                                revh=&regrev[p];

                                        }else{
                                            //we fix rev
                                            p=rescue_pair(&regrev[0],regfwd,n_regfwd,ilib);
                                            revh=&regrev[0];
                                            if(p>=0)
                                                fwdh=&regfwd[p];
                                        }
                                        //we rescue the pair
                                        if(p>=0){

                                            //printf("Rescue pair %lld %s %s %d %d %d %d %d %d %d %d %d %d %d %d\n",localseqs[i].id,mi->seq[fwdh->rid].name,mi->seq[revh->rid].name,
                                            //   fwdh->rs,fwdh->re,revh->rs,revh->re,fwdh->rev,revh->rev,fwdh->qs,fwdh->qe,revh->qs,revh->qe, fwdh->rev > 0 ? fwdh->re:fwdh->rs, revh->rev > 0 ? revh->re:revh->rs);
                                            //coverage
                                            fcov=(int)(((float)fwdh->blen/(float)localseqs[i].length) * 100);
                                            rcov=(int)(((float)revh->blen/(float)localseqr[i].length) * 100);
                                            //printf("%d %d\n",fcov,rcov );
                                            //we check that the aligment cover the query
                                            if(fcov >=qcov && rcov>=qcov && fcov <=110 && rcov<=110 && fwdh->mapq >= minq && revh->mapq >= minq){
                                                //we add the reverse first to the hits pool
                                                strncpy(kfwdseq, localseqr[i].seq, kmer_size);
                                                kfwdseq[kmer_size]='\0';
                                                int fr=0;
                                                if(revh->rev){
                                                    fr=16;
                                                }
                                                push_hit(&headhits[0],localseqr[i].id,mi->seq[revh->rid].name,2,fr,fr > 0 ? revh->re:revh->rs,revh->mapq,revh->mlen,revh->blen,kfwdseq,rcov);
                                                //we add the fwd
                                                int fo=0;
                                                if(fwdh->rev){
                                                    fo=16;
                                                }
                                                strncpy(kfwdseq, localseqs[i].seq, kmer_size);
                                                kfwdseq[kmer_size]='\0';
                                                push_hit(&headhits[0],localseqs[i].id,mi->seq[fwdh->rid].name,1,fo,fo > 0 ? fwdh->re:fwdh->rs,fwdh->mapq,fwdh->mlen,fwdh->blen,kfwdseq,fcov);

                                                aling_pairs++;
                                                rescue++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        //we clean the hits
                        if(n_regfwd > 0)
                            free(regfwd);
                        if(n_regrev > 0)
                            free(regrev);
            }
            //we printnt the results
            pthread_mutex_lock(&pres);
            //we store the hits in the files
            phit_t* nh = NULL;//next hit
            for (int f = 0; f < libcount; ++f)
            {
                //this way print the reverse before the foward
                //the header is not null
                nh=headhits[f];
                while(nh != NULL){
                    //fprintf(outshort,"%lld/%d\t%d\t%s\t%d\t%d\t%dM\t*\t%d\t%d\t%s\t%d\n",
                    //        nh->lrid,nh->plr,nh->orientation,nh->ctgid,nh->pos,nh->mapq,kmer_size,nh->mlen,nh->blen,nh->k_seq,nh->qcov);
                    //fprintf(outshort,"%lld\t%d\t%s\t%d\t%d\t%dM\t*\t%d\t%d\t%s\t%d\n",
                     //       (long long)nh->lrid,nh->orientation,nh->ctgid,nh->pos,nh->mapq,kmer_size,nh->mlen,nh->blen,nh->k_seq,nh->qcov);
                     //read_id,orientation,ctg_id,pos,mapq,ksize,qcov     	
                    fprintf(outshort,"%lld\t%d\t%s\t%d\t%d\t%d\t%d\n",
                            (long long)nh->lrid,nh->orientation,nh->ctgid,nh->pos,nh->mapq,kmer_size,nh->qcov);
                    headhits[f]=headhits[f]->next;
                    //we free the allocated hits
                    //free(nh->lrid);
                    free(nh->ctgid);
                    free(nh->k_seq);
                    free(nh);
                    nh=headhits[f];
                }
            }

            //if(number_long_reads%4995==0)
            if(pair_rescue)
                fprintf(stderr,"Thread_info:Partial TOTAL_PAIRS=%lld MAPPED_PAIRS=%lld  RL=%d MINQ=%d KS=%d MW=%d PRESCUE=%lld RESCUE=%lld\n",
                    (long long)number_read_pairs,(long long)aling_pairs,read_len,minq,kmer_size,moving_w,(long long)prescue,(long long)rescue);
            else
                fprintf(stderr,"Thread_info:Partial TOTAL_PAIRS=%lld MAPPED_PAIRS=%lld  RL=%d MINQ=%d KS=%d MW=%d\n",
                        (long long)number_read_pairs,(long long)aling_pairs,read_len,minq,kmer_size,moving_w);


            //we clean the container of reads before allocate more memory in a single threads fashion
            for (size_t i = 0; i < rseq; i++) {
                free(localseqs[i].seq);
                free(localseqs[i].name);
                free(localseqr[i].seq);
                free(localseqr[i].name);
            }
            pthread_mutex_unlock(&pres);

        }else{
            getjob=0;
        }
    } //endwhile
//we print some thread information
    pthread_mutex_lock(&pres);
    //print when pair rescue is enabled
    if(pair_rescue)
        fprintf(stderr,"Final_Thread_info:Total TOTAL_PAIRS=%lld MAPPED_PAIRS=%lld  RL=%d MINQ=%d KS=%d MW=%d PRESCUE=%lld RESCUE=%lld\n",
            (long long)number_read_pairs,(long long)aling_pairs,read_len,minq,kmer_size,moving_w,(long long)prescue,(long long)rescue);
    else
        fprintf(stderr,"Final_Thread_info:Total TOTAL_PAIRS=%lld MAPPED_PAIRS=%lld  RL=%d MINQ=%d KS=%d MW=%d\n",
                (long long)number_read_pairs,(long long)aling_pairs,read_len,minq,kmer_size,moving_w);

    pthread_mutex_unlock(&pres);
//all done in this thread
    return NULL;
}


int compare (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}



pe_stats infer_insert_size(char* fwdfile, char* revfile,  aligment_parameters *alnp, mm_idx_t *mi, mm_mapopt_t mopt){
    // open query file for reading; you may use your favorite FASTA/Q parser
    gzFile fwd = gzopen(fwdfile, "r");
    //printf("%d",f);
    if (! fwd) {
        fprintf (stderr, "gzopen of '%s' failed: %s.\n", fwdfile,
                 strerror (errno));
        exit (EXIT_FAILURE);
    }

    gzFile rev = gzopen(revfile, "r");
    //printf("%d",f);
    if (! rev) {
        fprintf (stderr, "gzopen of '%s' failed: %s.\n", revfile,
                 strerror (errno));
        exit (EXIT_FAILURE);
    }

    if (strcmp(fwdfile,revfile) == 0) {
        fprintf (stderr, "fwd and rev file are the same fwd=%s  rev=%s.\n", fwdfile,
                 revfile);
        exit (EXIT_FAILURE);
    }

    //assert(f);
    kseq_t *ks = kseq_init(fwd);
    kseq_t *kr = kseq_init(rev);

    //pointer to regions
    mm_reg1_t *regfwd; //fwd
    mm_reg1_t *regrev; //rev
    //pointers to hits
    mm_reg1_t *fwdh;
    mm_reg1_t *revh;
    int n_regfwd=0; // number of regions to consider fwd
    int n_regrev=0; // number of regions to consider rev
    int rseq=0;
    int mpairs=50000;

    mm_tbuf_t **tbuf;// = mm_tbuf_init(n_threads); // thread buffer; for multi-threading, allocate one tbuf for each thread
    tbuf = (mm_tbuf_t**)calloc(1, sizeof(mm_tbuf_t*));
    tbuf[0] = mm_tbuf_init();
    int insobs[50000] = {0};
    //we read both sequences and know we map them
    while(rseq < mpairs && kseq_read(ks) >=0 && kseq_read(kr) >=0 ){
        regfwd = mm_map(mi, ks->seq.l, ks->seq.s, &n_regfwd, tbuf[0], &mopt, 0); // get all hits for the query
        regrev = mm_map(mi, kr->seq.l, kr->seq.s, &n_regrev, tbuf[0], &mopt, 0); // get all hits for the query

        if(n_regfwd == 1 && n_regrev ==1){
            fwdh=&regfwd[0];
            revh=&regrev[0];
            if(revh->rid == fwdh->rid) {
                //coverage
                int fcov = (int) (((float) fwdh->blen / (float) ks->seq.l) * 100);
                int rcov = (int) (((float) revh->blen / (float) kr->seq.l) * 100);

                if (fcov > 65 && rcov > 65 && fwdh->mapq >= 30 && revh->mapq >= 30) {
                    // us a good pair aligned to the same cluster
                    int d1= fwdh->rev ? fwdh->re:fwdh->rs;
                    int d2= revh->rev ? revh->re:revh->rs;
                    //printf("%d %d %d %d %d\n", fwdh->rid, revh->rid, d1, d2, abs(d1-d2));
                    if(rseq < mpairs)
                        insobs[rseq]=abs(d1-d2);
                    rseq++;// is a good hit
                    //we save the distance among the pairs
                }
            }
        }

    }
    //we close the files and destriy the buffer;
    //we destroy the thhread buffers
    for(int i=0; i<1; i++) mm_tbuf_destroy(tbuf[i]);
    //end of threading code
    kseq_destroy(ks); // close the fwd file
    kseq_destroy(kr); // close the rev file
    gzclose(fwd);//flose gzfile
    gzclose(rev);//flose gzfile

    //we should return the std and avg of the insert size
    qsort (insobs, mpairs, sizeof(int), compare);
    //we remove 10% outlayers
    int beg=(int)(mpairs * 0.1);
    int end=mpairs-beg;
    //we compute the average
    uint64_t  avg=0;
    for (int n=beg; n<end+1; n++) {
        //printf("%d %d\n", n, insobs[n]);
        avg+=insobs[n];
    }

    int average = (int) (avg / (end - beg));
    int64_t std=0;
    //we compute the std of the sample
    for (int n=beg; n<end+1; n++) {
        std+=pow(insobs[n]-average,2);
    }
    int standar_dev=(int)sqrt(std/(end - beg));

    fprintf(stderr,"Average insert size %d std %d [-10 outlayers]\n",average, standar_dev);

    pe_stats a;
    a.avg=average;
    a.std=standar_dev;

    return a;
}



int map_short_read_files(char* prefix, char* fwdfile, char* revfile,  aligment_parameters *alnp, mm_idx_t *mi, mm_mapopt_t mopt){
    printf("Forward file  = %s\n",fwdfile);
    printf("Reverse file  = %s\n",revfile);
    // open query file for reading; you may use your favorite FASTA/Q parser
    gzFile fwd = gzopen(fwdfile, "r");
    //printf("%d",f);
    if (! fwd) {
        fprintf (stderr, "gzopen of '%s' failed: %s.\n", fwdfile,
                 strerror (errno));
        exit (EXIT_FAILURE);
    }

    gzFile rev = gzopen(revfile, "r");
    //printf("%d",f);
    if (! rev) {
        fprintf (stderr, "gzopen of '%s' failed: %s.\n", revfile,
                 strerror (errno));
        exit (EXIT_FAILURE);
    }

    if (strcmp(fwdfile,revfile) == 0) {
        fprintf (stderr, "fwd and rev file are the same fwd=%s  rev=%s.\n", fwdfile,
                 revfile);
        exit (EXIT_FAILURE);
    }


    //we infer the values of the library from 50000 pairs mapped within contigs
    pe_stats ilib=infer_insert_size(fwdfile,revfile,alnp,mi,mopt);
    //fprintf(stderr,"Average insert size %d std %d [-10 outlayers,]\n",ilib.avg, ilib.std);

    //exit(0);
    //assert(f);
    kseq_t *ks = kseq_init(fwd);
    kseq_t *kr = kseq_init(rev);

    char filename[100];//max length of file name
    memset(filename, '\0', sizeof(filename));
    //sprintf(filename, "libshort.%s.fastmin-sg.sam",prefix);
    sprintf(filename, "%s.fm.sam",prefix);
    FILE* outlong = fopen(filename,"w");
//we start the mapping of the reads

//we create the threads
    if (pthread_mutex_init(&gseq, NULL) != 0){
        printf("\n mutex for file reading failed\n");
        return 1;
    }
    if (pthread_mutex_init(&pres, NULL) != 0){
        printf("\n mutex for printing  failed\n");
        return 1;
    }
    //we start the mapping
    pthread_t *tid=(pthread_t *)malloc(alnp->n_threads * sizeof(pthread_t ));

    mm_tbuf_t **tbuf;// = mm_tbuf_init(n_threads); // thread buffer; for multi-threading, allocate one tbuf for each thread
    tbuf = (mm_tbuf_t**)calloc(alnp->n_threads, sizeof(mm_tbuf_t*));
    for(int i =0; i < alnp->n_threads; ++i) tbuf[i] = mm_tbuf_init();

    for(int i=0; i < alnp->n_threads; i++){
        thread_data_struct *data = malloc(sizeof *data);
        data->ks=ks;
        data->kr=kr;
        //data->outFiles=outFiles;
        data->outlong=outlong;
        data->mi=mi;
        //data->libcount=libcount;
        data->mopt=mopt;
        data->tbuf=tbuf[i];
        //mapping variables
        data->alnp=alnp;
        data->pes=&ilib;
        int err = pthread_create(&tid[i], NULL, &mapshortreads, data);
        if (err != 0){
            printf("\ncan't create thread :[%s]", strerror(err));
        }
    }
    //we wait for the thread to finish
    for(int i=0; i< alnp->n_threads; ++i) pthread_join(tid[i], NULL);
    //we destroy the thhread buffers
    for(int i=0; i<alnp->n_threads; i++) mm_tbuf_destroy(tbuf[i]);
    //end of threading code
    kseq_destroy(ks); // close the fwd file
    kseq_destroy(kr); // close the rev file
    gzclose(fwd);//flose gzfile
    gzclose(rev);//flose gzfile
    // todo: be a variable file ?
    fclose(outlong);//close the output samfile
    return 0;
}



/* main function */
static int usage() {
    fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   FastMin-SG <preset> [options] <reference> <query>\n");
	fprintf(stderr, "Version: 0.1\n\n");
	fprintf(stderr, "Preset:  ontraw    build synthetic libraries from Oxford Nanopore reads in FASTA/Q\n");
    fprintf(stderr, "         ontlon    build synthetic libraries from Ultra Long Nanopore reads in FASTA/Q\n");
	fprintf(stderr, "         pacraw    build synthetic libraries from PacBio reads in FASTA/Q\n");
	fprintf(stderr, "         pacccs    build synthetic libraries from PacBio-CCS reads in FASTA/Q\n");
    fprintf(stderr, "         shortr    Pair-end short reads mapping\n\n");
    fprintf(stderr, "Options: (Override the presets)\n");
    fprintf(stderr, "  Indexing:\n");
    fprintf(stderr, "    -H           use homopolymer-compressed k-mer (preferrable for PacBio)\n");
    fprintf(stderr, "    -k INT       k-mer size (no larger than 28) [15]\n");
    fprintf(stderr, "    -w INT       minizer window size [10]\n");
    //fprintf(stderr, "    -I NUM       split index for every ~NUM input bases [4G]\n");
    fprintf(stderr, "  Mapping synthetic read-pairs:\n");
    fprintf(stderr, "    -I list      Insert sizes for synthetic libraries [i.e 500,1000,2000,3000,4000,5000, .. ,20000]\n");
    fprintf(stderr, "    -L INT       Minimum length of long reads [>2000]\n");
    fprintf(stderr, "    -l INT       Length of synthetic reads[200]\n");
    fprintf(stderr, "    -q INT       Minimum quality score (no larger than 60) [40]\n");
    fprintf(stderr, "    -m INT       Moving windown [100]\n");
    fprintf(stderr, "    -c INT       Minimum aligment coverage[65]\n");
    fprintf(stderr, "    -s INT       Pair rescue in short-read alignment [1=true,0=false]\n");
    fprintf(stderr, "    -r INT       Number of read loaded to memory for mapping [long=1000,short=100000]\n");
    fprintf(stderr, "    -n INT       Start long-read counter [def=0]\n");
    fprintf(stderr, "  Other options:\n");
    fprintf(stderr, "    -t INT       Number of thread[1]\n");
    fprintf(stderr, "    -p STR       Ouput prefix [FM]\n");
	fprintf(stderr, "\n");
	return 1;
}


int main(int argc, char *argv[])
{
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;

    aligment_parameters *alnp = malloc(sizeof(aligment_parameters));
    alnp->n_threads=1;
    alnp->read_len=250;
    alnp->minq=40;
    alnp->kmer_size=15;
    alnp->minimizer_d=10;
    alnp->moving_w=150;
    alnp->qcov=65;
    alnp->numseq2read=500;
    alnp->min_len_long_read=2000;
    alnp->pair_rescue=1;


    char linserts[1000]; //= "500,1000,2000,3000,4000,5000,6000,7000,8000,10000,15000,20000,30000,40000,50000"; //libraries
    char linsertscpy[1000];
    char prefix[100];// = "FM"; //output prefix
    memset(linserts, '\0', sizeof(linserts));
    memset(prefix, '\0', sizeof(prefix));
    //we copy the prefix
    strcpy(prefix, "FM");
    strcpy(linserts, "500,1000,2000,3000,4000,5000,6000,7000,8000,10000,15000,20000,30000,40000,50000");
    //set the global variable to rename the long reads
    total_long_reads=0;


    mm_verbose = 3; // disable message output to stderr
	mm_set_opt(0, &iopt, &mopt);
    //we set the size of the batch 4Gb
    iopt.mini_batch_size = 4000000000ULL;

  if (argc < 4) return usage();
  if (strcmp(argv[1], "ontraw") == 0){
        iopt.k = 15, iopt.w = 5;
        alnp->kmer_size=15;
        alnp->minimizer_d=5;
        alnp->numseq2read=1000; //more longer reads
  }else if (strcmp(argv[1], "ontlon") == 0){
        iopt.k = 15, iopt.w = 5;
        alnp->kmer_size=15, alnp->minimizer_d=5;
        alnp->numseq2read=500; //more longer reads
      memset(linserts, '\0', sizeof(linserts));
      strcpy(linserts, "500,1000,2000,3000,4000,5000,6000,7000,8000,10000,15000,20000,30000,40000,50000,60000,70000,80000,90000,100000,120000,150000,180000,200000");
  }
  else if (strcmp(argv[1], "pacraw") == 0){
      iopt.k = 19, iopt.w = 5;
      alnp->kmer_size=19, alnp->minimizer_d=5;
      //compressed homopolymers enabled for raw pacbio reads
      iopt.flag |= MM_I_HPC;
      alnp->numseq2read=1000;
      memset(linserts, '\0', sizeof(linserts));
      strcpy(linserts, "500,1000,2000,3000,4000,5000,6000,7000,8000,10000,15000,20000");
  }
  else if (strcmp(argv[1], "pacccs") == 0){
    //should be similar to short reads mapping
    iopt.k = 21, iopt.w = 10;
    alnp->kmer_size=21, alnp->minimizer_d=10;
    alnp->numseq2read=1000;
    memset(linserts, '\0', sizeof(linserts));
    strcpy(linserts, "500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000");
  }
  else if (strcmp(argv[1], "shortr") == 0){
    //fprintf(stdout,"Preset=%s\n",argv[1]);
    iopt.k = 21, iopt.w = 10; alnp->numseq2read=50000;
    alnp->kmer_size=21, alnp->minimizer_d=10;
    alnp->minq=20;
    alnp->qcov=50; //min qcov for short reads
    mopt.mid_occ = 1000; //values from minimap2 for short reads
    mopt.max_occ = 5000;
    memset(linserts, '\0', sizeof(linserts));
    strcpy(linserts, "");
      //iopt.flag |= MM_F_SR | MM_F_FRAG_MODE | MM_F_NO_PRINT_2ND | MM_F_2_IO_THREADS | MM_F_HEAP_SORT;
  }
  else {
    fprintf(stderr, "[main] unrecognized mode '%s'. Abort!\n", argv[1]);
    return 1;
  }

  //we should check for advanced options to override the presets and defaults variables
  int c=0;
  extern char *optarg;
  extern int optind, optopt;
  while ((c = getopt (argc-1, argv+1, "Hk:w:I:l:q:m:c:r:t:p:s:n:")) != -1) {
      //fprintf(stderr,"%c=%s\n",c,optarg);
  		switch (c) {
  			case 'H': iopt.flag |= MM_I_HPC; break;
  			case 'k': alnp->kmer_size = atoi(optarg); iopt.k = alnp->kmer_size; break;
  			case 'w': alnp->minimizer_d = atoi(optarg); iopt.w = alnp->minimizer_d; break;
  			case 'I': memset(linserts, '\0', sizeof(linserts)); strcpy(linserts, optarg); break;
  			case 'l': alnp->read_len = atoi(optarg); break;
            case 'L': alnp->min_len_long_read= atoi(optarg); break;
            case 'q': alnp->minq = atoi(optarg); break;
            case 'm': alnp->moving_w = atoi(optarg); break;
            case 'c': alnp->qcov = atoi(optarg); break;
            case 's': alnp->pair_rescue = atoi(optarg); break;
            case 'r': alnp->numseq2read = atoi(optarg); break;
            case 't': alnp->n_threads = atoi(optarg); break;
            case 'n': total_long_reads = (int64_t)(atoi(optarg)); break;
            case 'p': memset(prefix, '\0', sizeof(prefix)); strcpy(prefix, optarg); break;
            case '?':
                    fprintf(stderr, "[%s] unrecognized variable '%c'. Abort!\n",argv[1],optopt);
            default:
                    return 1;
  		}
  }

  //check parameters
  if(argc-1-optind!= 2){
      usage();
      return 1;
  }else{
      fprintf(stderr,"LOG: Mapping mode =%s H=%d k=%d w=%d L=%d l=%d q=%d m=%d c=%d r=%d t=%d o=%s I=%s s=%d\n",
              argv[1],iopt.flag,alnp->kmer_size,alnp->minimizer_d,alnp->min_len_long_read,alnp->read_len,alnp->minq,
              alnp->moving_w,alnp->qcov,alnp->numseq2read,alnp->n_threads,prefix,linserts, alnp->pair_rescue);
      printf("reference=%s queryfile=%s\n",argv[optind+1],argv[optind+2]);
  }

  fprintf(stderr,"Building contig index\n");
    clock_t begin = clock();
    clock_t tbegin = clock();
    //we check the reference and the queryfiles
  // open index reader
	mm_idx_reader_t *r = mm_idx_reader_open(argv[optind+1], &iopt, 0);
	mm_idx_t *mi;
  //we create the database index
  mi = mm_idx_reader_read(r, alnp->n_threads);
  mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!

  if(r->n_parts != 1){
      fprintf(stderr,"Index of database is different from 1 = %d\n",r->n_parts);
      return 1;
  }
   //we print index stat to STDERR
    mm_idx_stat(mi);	
    clock_t end = clock();
    //fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
     //       __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
    fprintf(stderr,"Index construction time: %f seconds for %d target sequence(s)\n",(double)(end - begin) / CLOCKS_PER_SEC,mi->n_seq);

    //fprintf(stderr,"Contig index done!\n");

    //we will read the
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(argv[optind+2], "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    while ((read = getline(&line, &len, fp)) != -1) {
        //we remove the \n and split the line of the
        line[strcspn(line, "\n")] = 0;
        char** tokens = str_split(line, ' ');
        int i;
        for (i = 0; *(tokens + i); i++);

        if(strcmp(argv[1], "shortr") != 0){
            //we compute the size
            if(i == 2) {
                strcpy(linsertscpy, linserts);
                begin = clock();
                map_long_read_files(linsertscpy,tokens[1],tokens[0],alnp, mi,mopt);
                end = clock();
                fprintf(stderr,"Mapping (CPU) time: %f seconds using %d thread(s) for query %s file \n",(double)(end - begin) / CLOCKS_PER_SEC,alnp->n_threads,tokens[1]);
            }else{
                fprintf(stderr,"Error format of long-read files!\n");
                fprintf(stderr,"<LID> <LONGREAD>.fq.gz\n");
                return 1;
            }
            //free(tokens);
        }else{
            //printf("Implementing shortr....\n");
            //we compute the size
            if(i == 3) {
                begin = clock();
                map_short_read_files(tokens[0],tokens[1],tokens[2],alnp,mi,mopt);
                end = clock();
                fprintf(stderr,"Mapping (CPU) time: %f seconds using %d thread(s) for query %s %s file \n",(double)(end - begin) / CLOCKS_PER_SEC,alnp->n_threads,tokens[1],tokens[2]);
            }else{
                fprintf(stderr,"Error format of short-read files!\n");
                fprintf(stderr,"<ID> <FOWARD>.fq.gz <REVERSE>.fq.gz\n");
                return 1;
            }
            //return 1;
        }

        free(tokens);
    }

    fclose(fp);

    //end of fastmin-sg
    mm_idx_destroy(mi);
    mm_idx_reader_close(r); // close the index reader
    clock_t tend = clock();
    fprintf(stderr,"Total time: %f seconds using %d thread(s)\n",(double)(tend - tbegin) / CLOCKS_PER_SEC,alnp->n_threads);
    return 0;

}

/* API DOCS
typedef struct {
int32_t id;             // ID for internal uses (see also parent below)
int32_t cnt;            // number of minimizers; if on the reverse strand
int32_t rid;            // reference index; if this is an alignment from inversion rescue
int32_t score;          // DP alignment score
int32_t qs, qe, rs, re; // query start and end; reference start and end
int32_t parent, subsc;  // parent==id if primary; best alternate mapping score
int32_t as;             // offset in the a[] array (for internal uses only)
int32_t mlen, blen;     // seeded exact match length; seeded alignment block length
int32_t n_sub;          // number of suboptimal mappings
int32_t score0;         // initial chaining score (before chain merging/spliting)
uint32_t mapq:8, split:2, rev:1, inv:1, sam_pri:1, proper_frag:1, pe_thru:1, seg_split:1, seg_id:8, split_inv:1, dummy:7;
uint32_t hash;
float div;
mm_extra_t *p;
} mm_reg1_t;

//alingment presents
 if (strcmp(preset, "ava-ont") == 0) {
             io->flag = 0, io->k = 15, io->w = 5;
             mo->flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN;
             mo->min_chain_score = 100, mo->pri_ratio = 0.0f, mo->max_gap = 10000, mo->max_chain_skip = 25;
             mo->bw = 2000;
     } else if (strcmp(preset, "ava-pb") == 0) {
             io->flag |= MM_I_HPC, io->k = 19, io->w = 5;
             mo->flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN;
             mo->min_chain_score = 100, mo->pri_ratio = 0.0f, mo->max_gap = 10000, mo->max_chain_skip = 25;
    }
*/
