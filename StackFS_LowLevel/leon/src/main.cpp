#define FUSE_USE_VERSION 30
#define _XOPEN_SOURCE 500
//#define _GNU_SOURCE
#include <stdarg.h>
#include <unistd.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <errno.h>
#include <fuse3/fuse.h>
#include <assert.h>
#include <fuse3/fuse_lowlevel.h>
#include <stddef.h>
#include <fcntl.h> /* Definition of AT_* constants */
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <pthread.h>
#include <sys/xattr.h>
#include <sys/syscall.h>
#include <glib.h>
#include <Leon.hpp>
#include <unordered_set>

using namespace std;
#define EXT_LEN 5
#define THRESHOLD 21529800
#define CACHE_SIZE 20
#define EXT ".leon"
#define CONFIG "/.biofs.dirinfo"
#define TMP_FILE "/.file_compress"
#define CONFIG_LEN 15
#define TXT_EXT_LEN 4
#define TXT_EXT ".txt"
FILE *logfile;
#define ATTR "user.format"
#define ORIG_ATTR "user.orig_format"
#define FASTA_SIZE_ATTR "user.fasta_size"
#define FASTQ_SIZE_ATTR "user.fastq_size"
#define SIZE_LEN 15
#define ATTR_SIZE 6
#define FASTA "fasta"
#define FASTQ "fastq"
#define LEON "leon"
#define FASTA_F 2
#define FASTQ_F 1
#define LEON_F 3
#define NONE_F 0
#define USE_SPLICE 0
#define PG_SIZE 4096
#define TRACE_FILE "/trace_stackfs.log"
#define TRACE_FILE_LEN 18

pthread_spinlock_t spinlock; /* Protecting the above spin lock */
char banner[4096];
GHashTable* file_table;
GHashTable* file_cache;
unordered_map<string , unordered_set<string>> directory_list_map;
unordered_map<string, int> created_files_map; // new files and format
struct file_node{
	char file_name[PATH_MAX];
	Leon *leon;
	int total_size;
	bool isFastq;
	int format;
	struct file_node* next;
	struct file_node* prev;
	struct file_page* head;
};
struct file_node* file_node_head, *file_node_end;
struct file_pages{
	char* file_page;
	Leon *leon;
	int s_off;
	int e_off;
	int block_id;
	bool isFastq;
	int format;
	long size;
	int dirty;
	bool sequenceAdded;
        int seq_size_A;
        int last_file_size_A;
 	int seq_size_Q;
        int last_file_size_Q;
	struct file_pages* prev_page;
	struct file_pages* next_page;
};
static void splitFiles(int argc, char* argv[]);

void print_usage(void)
{
	printf("USAGE	: ./StackFS_ll -r <rootDir>|-rootdir=<rootDir> ");
	printf("[--attrval=<time(secs)>] [--statsdir=<statsDirPath>] ");
	printf("<mountDir> [FUSE options]\n"); /* For checkPatch.pl */
	printf("<rootDir>  : Root Directory containg the Low Level F/S\n");
	printf("<attrval>  : Time in secs to let kernel know how muh time ");
	printf("the attributes are valid\n"); /* For checkPatch.pl */
	printf("<statsDirPath> : Path for copying any statistics details\n");
	printf("<mountDir> : Mount Directory on to which the F/S should be ");
	printf("mounted\n"); /* For checkPatch.pl */
	printf("Example    : ./StackFS_ll -r rootDir/ mountDir/\n");
}

int log_open(char *statsDir)
{
	char *trace_path = NULL;

	if (statsDir) {
		trace_path = (char *)malloc(strlen(statsDir) +
				TRACE_FILE_LEN + 1);
		memset(trace_path, 0, strlen(statsDir) + TRACE_FILE_LEN + 1);
		strncpy(trace_path, statsDir, strlen(statsDir));
		strncat(trace_path, TRACE_FILE, TRACE_FILE_LEN);
	} else {
		trace_path = (char *)malloc(TRACE_FILE_LEN);
		memset(trace_path, 0, TRACE_FILE_LEN);
		strncpy(trace_path, TRACE_FILE + 1, TRACE_FILE_LEN-1);
	}
	printf("Trace file location : %s\n", trace_path);
	logfile = fopen(trace_path, "w");
	if (logfile == NULL) {
		perror("logfile");
		free(trace_path);
		return -1;
	}
	free(trace_path);
	setvbuf(logfile, NULL, _IOLBF, 0);
	return 0;
}

string formatSize(int size){
	string ret = to_string(size);
	for(int i=ret.size();i<SIZE_LEN;i++){
		ret = "0"+ret;
	}
	printf("formatted file size : %s\n", ret.c_str());
	return ret;
}


//Not Used currently
void splitCopy(int argc, char* argv[],const char* tmp_file, int block_id){
	Leon *leon = new Leon();
	leon->run (argc, argv);
	string dir = System::file().getDirectory(leon->_inputFilename);
	int block_count = block_id;
	bool isFasta = true, noHeader = false;
	IBank* whole_bank = Bank::open(tmp_file);
	int64_t seqCount = whole_bank->estimateNbItems();
	string temp_file(".tmp_file_compress");
	if(leon->_inputFilename.find(".fq") !=  string::npos || leon->_inputFilename.find(".fastq") !=  string::npos){
		if(! leon->getParser()->saw (Leon::STR_DNA_ONLY) && ! leon->getParser()->saw (Leon::STR_NOQUAL)){
			isFasta = false;
		}
	}
	if(!isFasta)
		temp_file += ".fastq";
	else
		temp_file += ".fasta";
	if(leon->getParser()->saw (Leon::STR_NOHEADER))
		noHeader = true;
	temp_file = dir+"/"+temp_file;
	Iterator<Sequence>* itSeq = leon->createIterator<Sequence> (whole_bank->iterator(),seqCount,
			"Creating Temp Files"
			);
	string line;
	//leon->orig_block_size->clear();
	//leon->seq_per_block->clear();
	bool more_lines = true;
	bool first_time = true;
	itSeq->first();
	while(! itSeq->isDone()){
		Leon *leon1  = new Leon();
		leon1->run (argc, argv);
		ofstream fout;
		int j = 0, size = 0, readid = 0;
		bool  reading = true;
		fout.open(temp_file.c_str());   //create a new file for each block
		if (!fout.good()){
			cerr << "I/O error while reading file." << endl;
		}
		while(! itSeq->isDone() && j < leon->READ_PER_BLOCK)
		{
			stringstream sint;
			sint << readid;
			if(!noHeader)
			{
				string line = itSeq->item().getComment();
				if(isFasta)
					fout<<">";
				else
					fout<<"@";
				fout<<line<<'\n';
				size += 1+line.size()+1;
			}
			else
			{
				if(isFasta)
					fout<< "> "<<sint.str()<<'\n';
				else
					fout<< "@ "<<sint.str()<<'\n';
				readid++;
				size+= 2+sint.str().size()+1;
			}
			int dna_len = itSeq->item().getDataSize();
			char dna[dna_len+1];
			strncpy(dna, itSeq->item().getDataBuffer(), dna_len);
			dna[dna_len] = '\0';
			fout<<dna<<'\n';
			size+= dna_len+1;
			if( !isFasta)
			{
				string line = itSeq->item().getQuality();
				fout<<"+\n";
				fout<<line<<'\n';
				size+= 2+line.size()+1;
			}
			j++;
			itSeq->next();
		}
		//leon->orig_block_size->push_back(size);
		//leon->seq_per_block->push_back(j);
		fout.close();
		leon1->executeCompression(block_count, temp_file.c_str());
		block_count++;
		delete leon1;
	}
	//leon->saveConfig();
	delete leon;
}

void log_close(void)
{

	if (logfile)
		fclose(logfile);
}

static char** prep_args(char * name, bool compress, bool isFasta, bool isBlock, int &count){
	count = 5;
	if(isFasta)
		count++;
	char** argv = new char*[count];
	argv[0] = new char[7];
	strcpy(argv[0],"./leon");
	argv[1] = new char[6];
	strcpy(argv[1],"-file");
	argv[2] = new char[strlen(name)+1];
	strcpy(argv[2],name);
	argv[3] =  new char[3];
	if(compress){
		strcpy(argv[3], "-c");
	}else{
		strcpy(argv[3], "-d");
	}
	argv[4] = new char[3];
	if(isBlock)
		strcpy(argv[4], "-b");
	else
		strcpy(argv[4], "-e");
	if(isFasta){
		argv[5] = new char[8];
		strcpy(argv[5], "-noqual");
	}	
	return argv;
}

static char* createFile( char * name, char* ext){
	char* point = NULL;
	char* file_name = NULL;
	if((point = strrchr(name,'/')) != NULL) {
		struct stat buffer;
		int file_name_len = strlen(name)-strlen(point);
		file_name = (char *) malloc(file_name_len+strlen(ext)+1);
		strncpy(file_name, name, file_name_len);
		strncpy(file_name+file_name_len, ext, strlen(ext));
		file_name[file_name_len+strlen(ext)] = '\0';
	}else{
		file_name = (char *) malloc(strlen(ext)+1);
		strcpy(file_name, ext);
	}
	return file_name;
}

int64_t print_timer(void)
{
	struct timespec tms;

	if (clock_gettime(CLOCK_REALTIME, &tms)) {
		printf("ERROR\n");
		return 0;
	}
	int64_t micros = tms.tv_sec * 1000000;

	micros += tms.tv_nsec/1000;
	if (tms.tv_nsec % 1000 >= 500)
		++micros;
	return micros;
}

/* called with file lock */
int print_banner(void)
{
	int len;
	int64_t time;
	int pid;
	unsigned long tid;

	banner[0] = '\0';
	time = print_timer()/1000;
	pid = getpid();
	tid = syscall(SYS_gettid);
	if (time == 0)
		return -1;
	len = sprintf(banner, "Time : %"PRId64" Pid : %d Tid : %lu ",
			time, pid, tid);
	(void) len;
	fputs(banner, logfile);
	return 0;
}

void StackFS_trace(const char *format, ...)
{
	va_list ap;
	int ret = 0;

	/*lock*/
	pthread_spin_lock(&spinlock);
	if (logfile) {
		/*Banner : time + pid + tid*/
		ret = print_banner();
		if (ret)
			goto trace_out;
		/*Done with banner*/
		va_start(ap, format);
		vfprintf(logfile, format, ap);
		/*Done with trace*/
		fprintf(logfile, "\n");
	}
trace_out:
	/*unlock*/
	pthread_spin_unlock(&spinlock);
}

/*=============Hash Table implementation==========================*/

/* The node structure that we maintain as our local cache which maps
 * the ino numbers to their full path, this address is stored as part
 * of the value of the hash table */
struct lo_inode {
	struct lo_inode *next;
	struct lo_inode *prev;
	/* Full path of the underlying ext4 path
	 * correspoding to its ino (easy way to extract back) */
	char *name;
	/* Inode numbers and dev no's of
	 * underlying EXT4 F/s for the above path */
	ino_t ino;
	dev_t dev;
	/* inode number sent to lower F/S */
	ino_t lo_ino;
	/* Lookup count of this node */
	uint64_t nlookup;
};

#define HASH_TABLE_MIN_SIZE 8192

/* The structure is used for maintaining the hash table
 * 1. array	--> Buckets to store the key and values
 * 2. use	--> Current size of the hash table
 * 3. size	--> Max size of the hash table
 * (we start with NODE_TABLE_MIN_SIZE)
 * 4. split	--> used to resize the table
 * (this is how fuse-lib does) */
struct node_table {
	struct lo_inode **array;
	size_t use;
	size_t size;
	size_t split;
};

static int hash_table_init(struct node_table *t)
{
	t->size = HASH_TABLE_MIN_SIZE;
	t->array = (struct lo_inode **) calloc(1,
			sizeof(struct lo_inode *) * t->size);
	if (t->array == NULL) {
		fprintf(stderr, "fuse: memory allocation failed\n");
		return -1;
	}
	t->use = 0;
	t->split = 0;

	return 0;
}

void hash_table_destroy(struct node_table *t)
{
	free(t->array);
}

static int hash_table_resize(struct node_table *t)
{
	size_t newsize = t->size * 2;
	void *newarray = NULL;

	newarray = realloc(t->array, sizeof(struct lo_inode *) * newsize);
	if (newarray == NULL) {
		fprintf(stderr, "fuse: memory allocation failed\n");
		return -1;
	}

	t->array = (struct lo_inode**)newarray;
	/* zero the newly allocated space */
	memset(t->array + t->size, 0, t->size * sizeof(struct lo_inode *));
	t->size = newsize;
	t->split = 0;

	return 0;
}

/* The structure which is used to store the hash table
 * and it is always comes as part of the req structure */
struct lo_data {
	/* hash table mapping key (inode no + complete path) -->
	 *  value (linked list of node's - open chaining) */
	struct node_table hash_table;
	/* protecting the above hash table */
	pthread_spinlock_t spinlock;
	/* put the root Inode '/' here itself for faster
	 * access and some other useful raesons */
	struct lo_inode root;
	/* do we still need this ? let's see*/
	double attr_valid;
};

struct lo_dirptr {
	DIR *dp;
	struct dirent *entry;
	off_t offset;
};

static struct lo_dirptr *lo_dirptr(struct fuse_file_info *fi)
{
	return ((struct lo_dirptr *) ((uintptr_t) fi->fh));
}


static struct lo_data *get_lo_data(fuse_req_t req)
{
	return (struct lo_data *) fuse_req_userdata(req);
}

static struct lo_inode *lo_inode(fuse_req_t req, fuse_ino_t ino)
{
	if (ino == FUSE_ROOT_ID)
		return &get_lo_data(req)->root;
	else
		return (struct lo_inode *) (uintptr_t) ino;
}

static char *lo_name(fuse_req_t req, fuse_ino_t ino)
{
	return lo_inode(req, ino)->name;
}

/* This is what given to the kernel FUSE F/S */
static ino_t get_lower_fuse_inode_no(fuse_req_t req, fuse_ino_t ino) {
	return lo_inode(req, ino)->lo_ino;
}

/* This is what given to the user FUSE F/S */
static ino_t get_higher_fuse_inode_no(fuse_req_t req, fuse_ino_t ino) {
	return lo_inode(req, ino)->ino;
}


static double lo_attr_valid_time(fuse_req_t req)
{
	return ((struct lo_data *) fuse_req_userdata(req))->attr_valid;
}

static void construct_full_path(fuse_req_t req, fuse_ino_t ino,
		char *fpath, const char *path)
{
	strcpy(fpath, lo_name(req, ino));
	strncat(fpath, "/", 1);
	strncat(fpath, path, PATH_MAX);
}

/*======================End=======================================*/

/* Function which generates the hash depending on the ino number
 * and full path */
static size_t name_hash(struct lo_data *lo_data, fuse_ino_t ino,
		const char *fullpath)
{
	uint64_t hash = ino;
	uint64_t oldhash;
	const char *name;

	name = fullpath;

	for (; *name; name++)
		hash = hash * 31 + (unsigned char) *name;

	hash %= lo_data->hash_table.size;
	oldhash = hash % (lo_data->hash_table.size / 2);
	if (oldhash >= lo_data->hash_table.split)
		return oldhash;
	else
		return hash;
}

static void remap_hash_table(struct lo_data *lo_data)
{
	struct node_table *t = &lo_data->hash_table;
	struct lo_inode **nodep;
	struct lo_inode **next;
	struct lo_inode *prev;
	size_t hash;

	if (t->split == t->size / 2)
		return;

	// split this bucket by recalculating the hash 
	hash = t->split;
	t->split++;

	for (nodep = &t->array[hash]; *nodep != NULL; nodep = next) {
		struct lo_inode *node = *nodep;
		size_t newhash = name_hash(lo_data, node->ino, node->name);

		if (newhash != hash) {
			prev = node->prev;
			*nodep = node->next;
			if (*nodep)
				(*nodep)->prev = prev;

			node->prev = NULL;
			node->next = t->array[newhash];
			if (t->array[newhash])
				(t->array[newhash])->prev = node;
			t->array[newhash] = node;
			next = nodep;
		} else {
			next = &node->next;
		}
	}

	//If we have reached the splitting to half of the size
	// then double the size of hash table 
	if (t->split == t->size / 2)
		hash_table_resize(t);
}

static int insert_to_hash_table(struct lo_data *lo_data,
		struct lo_inode *lo_inode)
{
	size_t hash = name_hash(lo_data, lo_inode->ino, lo_inode->name);
	StackFS_trace("The insert file name: %s", lo_inode->name);
	lo_inode->next = lo_data->hash_table.array[hash];
	if (lo_data->hash_table.array[hash])
		(lo_data->hash_table.array[hash])->prev = lo_inode;
	lo_data->hash_table.array[hash] = lo_inode;
	lo_data->hash_table.use++;

	if (lo_data->hash_table.use >= lo_data->hash_table.size / 2)
		remap_hash_table(lo_data);

	return 0;
}
static void hash_table_reduce(struct node_table *t)
{
	size_t newsize = t->size / 2;
	void* newarray;

	if (newsize < HASH_TABLE_MIN_SIZE)
		return;

	newarray = realloc(t->array, sizeof(struct node *) * newsize);
	if (newarray != NULL)
		t->array = (struct lo_inode **)newarray;

	t->size = newsize;
	t->split = t->size / 2;
}

static void remerge_hash_table(struct lo_data *lo_data)
{
	struct node_table *t = &lo_data->hash_table;
	int iter;

	/*This means all the hashes would be under the half size
	 * of table (so simply make it half) */
	if (t->split == 0)
		hash_table_reduce(t);

	for (iter = 8; t->split > 0 && iter; iter--) {
		struct lo_inode **upper;

		t->split--;
		upper = &t->array[t->split + t->size / 2];
		if (*upper) {
			struct lo_inode **nodep;
			struct lo_inode *prev = NULL;

			for (nodep = &t->array[t->split];
					*nodep; nodep = &(*nodep)->next)
				prev = *nodep;

			*nodep = *upper;
			(*upper)->prev = prev;
			*upper = NULL;
			break;
		}
	}
}

static int delete_from_hash_table(struct lo_data *lo_data,
		struct lo_inode *lo_inode)
{
	struct lo_inode *prev, *next;

	prev = next = NULL;
	size_t hash = 0;

	pthread_spin_lock(&lo_data->spinlock);

	prev = lo_inode->prev;
	next = lo_inode->next;

	if (prev) {
		prev->next = next;
		if (next)
			next->prev = prev;
		goto del_out;
	} else {
		hash = name_hash(lo_data, lo_inode->ino, lo_inode->name);

		if (next)
			next->prev = NULL;

		lo_data->hash_table.array[hash] = next;
	}

del_out:
	// free the lo_inode  
	lo_inode->prev = lo_inode->next = NULL;
	free(lo_inode->name);
	free(lo_inode);

	lo_data->hash_table.use--;
	if (lo_data->hash_table.use < lo_data->hash_table.size / 4)
		remerge_hash_table(lo_data);

	pthread_spin_unlock(&lo_data->spinlock);
	return 0;
}

/* Function which checks the inode in the hash table
 * by calculating the hash from ino and full path */
static struct lo_inode *lookup_lo_inode(struct lo_data *lo_data,
		struct stat *st, const char *fullpath)
{
	size_t hash = name_hash(lo_data, st->st_ino, fullpath);
	struct lo_inode *node;

	for (node = lo_data->hash_table.array[hash]; node != NULL;
			node = node->next) {
		if ((node->ino == st->st_ino) && (node->dev == st->st_dev) &&
				(strcmp(node->name, fullpath) == 0))
			return node;
	}

	return NULL;
}

void free_hash_table(struct lo_data *lo_data)
{
	struct lo_inode *node, *next;
	int i;
	int size = lo_data->hash_table.size;
	for (i = 0; i < size; i++) {
		node = lo_data->hash_table.array[i];
		while (node) {
			next = node->next;
			/* free up the node */
			free(node->name);
			free(node);
			node = next;
		}
	}
}

/* A function which checks the hash table and returns the lo_inode
 * otherwise a new lo_inode is created and inserted into the hashtable
 * req		--> for the hash_table reference
 * st		--> to check against the ino and dev_id
 *			when navigating the bucket chain
 * fullpath	--> full path is used to construct the key */
struct lo_inode *find_lo_inode(fuse_req_t req, struct stat *st, char *fullpath)
{
	struct lo_data *lo_data;
	struct lo_inode *lo_inode;
	int res = 0;

	lo_data = get_lo_data(req);

	pthread_spin_lock(&lo_data->spinlock);

	lo_inode = lookup_lo_inode(lo_data, st, fullpath);

	if (lo_inode == NULL) {
		/* create the node and insert into hash_table */
		lo_inode =  (struct lo_inode *)calloc(1, sizeof(struct lo_inode));
		if (!lo_inode)
			goto find_out;
		lo_inode->ino = st->st_ino;
		lo_inode->dev = st->st_dev;
		lo_inode->name = strdup(fullpath);
		/* store this for mapping (debugging) */
		lo_inode->lo_ino = (uintptr_t) lo_inode;
		lo_inode->next = lo_inode->prev = NULL;

		/* insert into hash table */
		res = insert_to_hash_table(lo_data, lo_inode);
		if (res == -1) {
			free(lo_inode->name);
			free(lo_inode);
			lo_inode = NULL;
			goto find_out;
		}
	}
	lo_inode->nlookup++;
find_out:
	pthread_spin_unlock(&lo_data->spinlock);
	return lo_inode;
}

static int getSeqSize(string s,string file_format){
	string line="";
        stringstream ss(s);
	int lineCount = 0;
	int seq_size = 0;
	if(file_format.find(FASTA)!=string::npos){
		lineCount = 2;
	}else{
		lineCount = 4;
	}
	while(lineCount>0){
		if(getline(ss, line)){
			seq_size += line.size();
		}
		lineCount--;
	}
        return seq_size;
}

static void getFileCount(char* filebase, vector<int>& block_sizes, bool isFasta){
	string base(filebase);
        string filename;
        struct stat buffer;
	int count = 0;
	block_sizes.clear();
	string attr_name;
	if(isFasta)
		attr_name = FASTA_SIZE_ATTR;
	else
		attr_name = FASTQ_SIZE_ATTR;
        filename = base+"_0"+EXT;
        while(stat(filename.c_str(), &buffer)==0){
		char value[SIZE_LEN];
        	int ret = lgetxattr(filename.c_str(), attr_name.c_str(), value, SIZE_LEN);
		if(ret >0){
			block_sizes.push_back(stoi(value));
		}
		//StackFS_trace("Final file name : %s, count : %d",filename.c_str(), block_sizes[count]);
		count++;
		filename = base + "_" + to_string(count)+EXT;
	}
}

static long get_file_size(char* fullPath, string val){
	long ret =0l;
	vector<int> block_sizes;
	if(val.find(FASTA)!=string::npos){
		getFileCount(fullPath, block_sizes, true); //is Fasta true
	}else{
		getFileCount(fullPath, block_sizes, false); // is Fasta false
	}
	for(int i=0;i<block_sizes.size();i++){
		ret += (long)block_sizes[i];
	}
        //StackFS_trace("Final file name : %s, size : %d",fullPath, ret);
	return ret;
}

static void stackfs_ll_lookup(fuse_req_t req, fuse_ino_t parent,
		const char *name)
{
	struct fuse_entry_param e;
	int res;
	char *fullPath = NULL;
	char *file_name = NULL;
	double attr_val;

	StackFS_trace("Lookup called on name : %s, parent ino : %llu",
			name, parent);
	fullPath = (char *)malloc(PATH_MAX);
	construct_full_path(req, parent, fullPath, name);

	attr_val = lo_attr_valid_time(req);
	memset(&e, 0, sizeof(e));

	e.attr_timeout = attr_val;
	e.entry_timeout = 1.0; /* dentry timeout */

	generate_start_time(req);
	char * point = NULL;
	string filename(fullPath);
	char* fname = (char*) malloc(PATH_MAX);
	strcpy(fname, fullPath);
	strcat(fname, "_0.leon");
	char value[ATTR_SIZE];
	int ret = lgetxattr(fname, ATTR, value, ATTR_SIZE);
	string val(value);
	free(fname);
	fname = NULL;
	StackFS_trace("The file name: %s is looked up %s", fullPath, value);
	if(ret > 0 && (val.find(FASTA)!=string::npos || val.find(FASTQ)!=string::npos )) {
		struct stat buffer;
		bool isFastq = true;
		int len = strlen(fullPath)+EXT_LEN + 3;
		file_name = (char*) malloc(PATH_MAX);
		strcpy(file_name, fullPath);
		strcat(file_name, "_0");
		strcat(file_name, EXT);
		file_name[len-1] = '\0';
		if(stat(file_name, &buffer) == 0) {
			res = stat(file_name, &e.attr);
			e.attr.st_size = get_file_size(fullPath, val);
		}else{
			res = stat(fullPath, &e.attr);
		}
	}else if(ret > 0 && val.find(LEON)!=string::npos ){
                struct stat buffer;
                bool isFastq = false;
                file_name = (char*) malloc(PATH_MAX);
                /*string base_name(filename.substr(0, filename.size()-5));
                string fastq_name(base_name);
                fastq_name += ".fastq";
		time_t last;*/
                string real_name = filename + "_0.leon";
		string qual_name = filename + "_0.qual";
                if(stat(real_name.c_str(),&buffer) == 0){
                        int len = strlen(name) + 8;
                        isFastq = true;
                        strcpy(file_name, filename.c_str());
                        file_name[len-1] = '\0';
                        res = stat(real_name.c_str(),&e.attr);
			if(stat(qual_name.c_str(),&buffer) == 0){
				e.attr.st_size += buffer.st_size;
			}
                        real_name = filename + "_1.leon";
			qual_name = filename + "_1.qual";
                        int block_no = 1;
                        while(stat(real_name.c_str(),&buffer) == 0){
                               	e.attr.st_size += buffer.st_size;
                                block_no++;
                                real_name = filename + "_" +to_string(block_no)+ "./leon";
				qual_name = filename + "_" +to_string(block_no)+ "./qual";
				if(stat(qual_name.c_str(),&buffer) == 0){
                                	e.attr.st_size += buffer.st_size;
				}
                        }
                }
	}else{
		res = stat(fullPath, &e.attr);
	}
	generate_end_time(req);
	populate_time(req);
	//StackFS_trace(" The result of look up: %d %s", res, file_name);
	if (res == 0) {
		struct lo_inode *inode;
		inode = find_lo_inode(req, &e.attr, fullPath);

		if (fullPath)
			free(fullPath);

		if(file_name!=NULL){
			free(file_name);
			file_name = NULL;
		}

		if (!inode)
			fuse_reply_err(req, ENOMEM);
		else {
			/* store this address for faster path conversations */
			StackFS_trace("Successfully looked up");
			e.ino = inode->lo_ino;
			fuse_reply_entry(req, &e);
		}
	} else {
		if (fullPath)
			free(fullPath);
		fuse_reply_err(req, ENOENT);
	}
	if(file_name!=NULL){
		free(file_name);
		file_name = NULL;
	}
}

static void stackfs_ll_getattr(fuse_req_t req, fuse_ino_t ino,
		struct fuse_file_info *fi)
{
	int res;
	struct stat buf;
	(void) fi;
	double attr_val;
	char * name = lo_name(req, ino);
	char* file_name = NULL;
	StackFS_trace("Getattr called on name : %s and inode : %llu",
			lo_name(req, ino), lo_inode(req, ino)->ino);
	attr_val = lo_attr_valid_time(req);
	generate_start_time(req);
	string filename (name);
	char* fname = (char*) malloc(PATH_MAX);
        strcpy(fname, name);
        strcat(fname, "_0.leon");
	char value[ATTR_SIZE];
        int ret = lgetxattr(fname, ATTR, value, ATTR_SIZE);
	string val(value);
	free(fname);
	fname = NULL;
	if(ret > 0 && (val.find(FASTA)!=string::npos || val.find(FASTQ)!=string::npos )) {
		struct stat buffer;
		int len = strlen(name)+EXT_LEN + 3;
		file_name = (char*) malloc(PATH_MAX);
		strcpy(file_name, name);
		strcat(file_name, "_0");
		strcat(file_name, EXT);
		file_name[len-1] = '\0';
		if(stat(file_name, &buffer) == 0) {
			res = stat(file_name, &buf);
                        buf.st_size = get_file_size(name, val);
		}else{
			res = stat(name, &buf);
		}
	}else if(ret > 0 && val.find(LEON)!=string::npos ){
                struct stat buffer;
                bool isFastq = false;
                file_name = (char*) malloc(PATH_MAX);
                /*string base_name(filename.substr(0, filename.size()-5));
                string fastq_name(base_name);
                fastq_name += ".fastq";
		time_t last;*/
                string real_name = filename + "_0.leon";
                string qual_name = filename + "_0.qual";
                if(stat(real_name.c_str(),&buffer) == 0){
                        int len = strlen(name) + 8;
                        isFastq = true;
                        strcpy(file_name, filename.c_str());
                        file_name[len-1] = '\0';
                        res = stat(real_name.c_str(),&buf);
                        if(stat(qual_name.c_str(),&buffer) == 0){
                                buf.st_size += buffer.st_size;
			}
                        real_name = filename + "_1.leon";
                        qual_name = filename + "_1.qual";
                        int block_no = 1;
                        while(stat(real_name.c_str(),&buffer) == 0){
                                buf.st_size += buffer.st_size;
                                block_no++;
                                real_name = filename + "_" +to_string(block_no)+ "./leon";
                                qual_name = filename + "_" +to_string(block_no)+ "./qual";
                                if(stat(qual_name.c_str(),&buffer) == 0){
                                        buf.st_size += buffer.st_size;
				}
                        }
                }
	}else{
		res = stat(name, &buf);
	}
	if(file_name!=NULL){
		free(file_name);
		file_name = NULL;
	}
	generate_end_time(req);
	populate_time(req);
	if (res == -1)
		return (void) fuse_reply_err(req, errno);

	fuse_reply_attr(req, &buf, attr_val);
}

static void stackfs_ll_setattr(fuse_req_t req, fuse_ino_t ino,
		struct stat *attr, int to_set, struct fuse_file_info *fi)
{
	int res;
	(void) fi;
	struct stat buf;
	double attr_val;

	StackFS_trace("Setattr called on name : %s and inode : %llu",
			lo_name(req, ino), lo_inode(req, ino)->ino);
	attr_val = lo_attr_valid_time(req);
	generate_start_time(req);
	if (to_set & FUSE_SET_ATTR_SIZE) {
		/*Truncate*/
		res = truncate(lo_name(req, ino), attr->st_size);
		if (res != 0) {
			generate_end_time(req);
			populate_time(req);
			return (void) fuse_reply_err(req, errno);
		}
	}

	if (to_set & (FUSE_SET_ATTR_ATIME | FUSE_SET_ATTR_MTIME)) {
		/* Update Time */
		struct utimbuf tv;

		tv.actime = attr->st_atime;
		tv.modtime = attr->st_mtime;
		res = utime(lo_name(req, ino), &tv);
		if (res != 0) {
			generate_end_time(req);
			populate_time(req);
			return (void) fuse_reply_err(req, errno);
		}
	}

	memset(&buf, 0, sizeof(buf));
	res = stat(lo_name(req, ino), &buf);
	generate_end_time(req);
	populate_time(req);
	if (res != 0)
		return (void) fuse_reply_err(req, errno);

	fuse_reply_attr(req, &buf, attr_val);
}

static void stackfs_ll_create(fuse_req_t req, fuse_ino_t parent,
		const char *name, mode_t mode, struct fuse_file_info *fi)
{
	int fd, res;
	struct fuse_entry_param e;
	char *fullPath = NULL;
	double attr_val;

	StackFS_trace("Create called on %s and parent ino : %llu",
			name, lo_inode(req, parent)->ino);

	fullPath = (char *)malloc(PATH_MAX);
	construct_full_path(req, parent, fullPath, name);
	attr_val = lo_attr_valid_time(req);

	generate_start_time(req);

	fd = creat(fullPath, mode);

	if (fd == -1) {
		if (fullPath)
			free(fullPath);
		generate_end_time(req);
		populate_time(req);
		return (void)fuse_reply_err(req, errno);
	}

	memset(&e, 0, sizeof(e));

	e.attr_timeout = attr_val;
	e.entry_timeout = 1.0;

	res = stat(fullPath, &e.attr);
	generate_end_time(req);
	populate_time(req);

	if (res == 0) {
		/* insert lo_inode into the hash table */
		struct lo_data *lo_data;
		struct lo_inode *lo_inode;

		lo_inode = (struct lo_inode*)calloc(1, sizeof(struct lo_inode));
		if (!lo_inode) {
			if (fullPath)
				free(fullPath);

			return (void) fuse_reply_err(req, errno);
		}

		lo_inode->ino = e.attr.st_ino;
		lo_inode->dev = e.attr.st_dev;
		lo_inode->name = strdup(fullPath);
		/* store this for mapping (debugging) */
		lo_inode->lo_ino = (uintptr_t) lo_inode;
		lo_inode->next = lo_inode->prev = NULL;
		free(fullPath);

		lo_data = get_lo_data(req);
		pthread_spin_lock(&lo_data->spinlock);

		res = insert_to_hash_table(lo_data, lo_inode);

		pthread_spin_unlock(&lo_data->spinlock);

		if (res == -1) {
			free(lo_inode->name);
			free(lo_inode);
			fuse_reply_err(req, EBUSY);
		} else {
			lo_inode->nlookup++;
			e.ino = lo_inode->lo_ino;
			StackFS_trace("Create called, e.ino : %llu", e.ino);
			fi->fh = fd;
			fuse_reply_create(req, &e, fi);
			created_files_map[lo_inode->name] = NONE_F;
		}

	} else {
		if (fullPath)
			free(fullPath);
		fuse_reply_err(req, errno);
	}
}

static void stackfs_ll_mkdir(fuse_req_t req, fuse_ino_t parent,
		const char *name, mode_t mode)
{
	int res;
	struct fuse_entry_param e;
	char *fullPath = NULL;
	double attr_val;

	StackFS_trace("Mkdir called with name : %s, parent ino : %llu",
			name, lo_inode(req, parent)->ino);

	fullPath = (char *)malloc(PATH_MAX);
	construct_full_path(req, parent, fullPath, name);
	attr_val = lo_attr_valid_time(req);

	generate_start_time(req);
	res = mkdir(fullPath, mode);

	if (res == -1) {
		/* Error occurred while creating the directory */
		if (fullPath)
			free(fullPath);

		generate_end_time(req);
		populate_time(req);

		return (void)fuse_reply_err(req, errno);
	}

	/* Assign the stats of the newly created directory */
	memset(&e, 0, sizeof(e));
	e.attr_timeout = attr_val;
	e.entry_timeout = 1.0; /* may be attr_val */
	res = stat(fullPath, &e.attr);
	generate_end_time(req);
	populate_time(req);

	if (res == 0) {
		/* insert lo_inode into the hash table */
		struct lo_data *lo_data;
		struct lo_inode *lo_inode;

		lo_inode = (struct lo_inode*)calloc(1, sizeof(struct lo_inode));
		if (!lo_inode) {
			if (fullPath)
				free(fullPath);

			return (void) fuse_reply_err(req, errno);
		}

		lo_inode->ino = e.attr.st_ino;
		lo_inode->dev = e.attr.st_dev;
		lo_inode->name = strdup(fullPath);
		/* store this for mapping (debugging) */
		lo_inode->lo_ino = (uintptr_t) lo_inode;
		lo_inode->next = lo_inode->prev = NULL;
		free(fullPath);

		lo_data = get_lo_data(req);

		pthread_spin_lock(&lo_data->spinlock);

		res = insert_to_hash_table(lo_data, lo_inode);

		pthread_spin_unlock(&lo_data->spinlock);

		if (res == -1) {
			free(lo_inode->name);
			free(lo_inode);
			fuse_reply_err(req, EBUSY);
		} else {
			lo_inode->nlookup++;
			e.ino = lo_inode->lo_ino;
			fuse_reply_entry(req, &e);
		}
	} else {
		if (fullPath)
			free(fullPath);
		fuse_reply_err(req, errno);
	}
}

static char* change_file_ext(char * name, char* ext){
	char * point = NULL;
	char * file_name = NULL;
	int ext_len = sizeof(ext);
	if((point = strrchr(name,'.')) != NULL) {
		int file_name_len = strlen(name)-strlen(point);
		file_name = (char *) malloc(file_name_len+ext_len+1);
		strncpy(file_name, name, file_name_len);
		strncpy(file_name+file_name_len, ext, ext_len);
		file_name[file_name_len+ext_len] = '\0';
	}
	return file_name;
}
// name will contain the existing file name portion, an existing file XXXX.fastq_1.leon -> name XXXX.fastq
// we only require the orginal file portion and isFasta or Fastq information to operate
static void decompress_store(char * name, int fd, int format){
	int bytes_read = 0;
	char* key = (char *)malloc(strlen(name)+1);
	int offSet = 0;
	strcpy(key, name);
	if(!g_hash_table_contains(file_table,key)){
		struct file_pages * fp = (struct file_pages*) malloc(sizeof(struct file_pages));
		struct file_pages * current = fp;
		fp->prev_page = NULL;
		fp->next_page = NULL;
		Leon* leon = new Leon();
                string file_format;
		int count = 0;
		char** args;
		if(format ==FASTA_F){
			args = prep_args(name,false, true, true,count);
                        file_format = "FASTA";
		}else{
			args = prep_args(name,false, false, true,count);
                        file_format = "FASTQ";
		}
		leon->run(count, args);
		vector<string>* out =leon->executeDecompression(0);
		current->file_page = (char *) malloc ((*out)[0].size()+1);
		strcpy(current->file_page,(*out)[0].c_str());
		offSet = strlen(current->file_page)+1;
		current->s_off = 0;
		current->seq_size_A = -1;
		current->seq_size_Q = -1;
		current->last_file_size_A = -1;
		current->last_file_size_Q = -1;
		if(format ==FASTA_F){
                	current->seq_size_A = getSeqSize((*out)[0], file_format);
		}else{
			current->seq_size_Q = getSeqSize((*out)[0], file_format);
		}
		current->e_off = (*out)[0].size();
		current->block_id = 0;
		current->leon = leon;
		current->format = format;
		current->dirty = 0;
		current->size = (*out)[0].size();
		current->sequenceAdded = false;
		pthread_spin_lock(&spinlock);
		delete out;
		g_hash_table_insert(file_table, key, fp);
		pthread_spin_unlock(&spinlock);
	}else{
		struct file_node * fn = (struct file_node*) g_hash_table_lookup(file_cache, key);
		if(fn!=NULL){
			StackFS_trace("Removing from Linked list: %s", key);
			pthread_spin_lock(&spinlock);
			fn->prev->next = fn->next;
			fn->next->prev = fn->prev;
			g_hash_table_remove(file_cache, key);
			free(fn);
			fn = NULL;
			pthread_spin_unlock(&spinlock);
		}
		if(key!=NULL){
			free(key);
			key = NULL;
		}
	}
}

static void splitFiles(int argc, char* argv[], int fast){
	Leon *leon = new Leon();
	leon->run (argc, argv);
	string dir = System::file().getDirectory(leon->_inputFilename);
	int block_count = 0;
	bool isFasta = true, noHeader = false;
	string temp_file(".tmp_file_compress");
	string fast_name = leon->_inputFilename;
        if(fast == FASTQ_F)
                isFasta = false;
        if(isFasta){
                temp_file += ".fasta";
		fast_name += ".fasta";
        }else{
                temp_file += ".fastq";
		fast_name += ".fastq";
	}
	rename(leon->_inputFilename.c_str(), fast_name.c_str());
	IBank* whole_bank = Bank::open(fast_name);
	int64_t seqCount = whole_bank->estimateNbItems();
	if(leon->getParser()->saw (Leon::STR_NOHEADER))
		noHeader = true;
	temp_file = dir+"/"+temp_file;
	Iterator<Sequence>* itSeq = leon->createIterator<Sequence> (whole_bank->iterator(),seqCount,
			"Creating Temp Files"
			);
	string line;
	//leon->orig_block_size->clear();
	//leon->fastq_block_size->clear();
	//leon->seq_per_block->clear();
	bool more_lines = true;
	bool first_time = true;
	itSeq->first();
	while(! itSeq->isDone()){
		Leon *leon1  = new Leon();
		leon1->run (argc, argv);
		ofstream fout;
		int j = 0, size = 0, readid = 0, fasta_size = 0;
		bool  reading = true;
		fout.open(temp_file.c_str());   //create a new file for each block
		if (!fout.good()){
			cerr << "I/O error while reading file." << endl;
		}
		while(! itSeq->isDone() && j < leon->READ_PER_BLOCK)
		{
			stringstream sint;
			sint << readid;
			if(!noHeader)
			{
				string line = itSeq->item().getComment();
				if(isFasta)
					fout<<">";
				else
					fout<<"@";
				fout<<line<<'\n';
				size += 1+line.size()+1;
				fasta_size += 1+line.size()+1;
			}
			else
			{
				if(isFasta)
					fout<< "> "<<sint.str()<<'\n';
				else
					fout<< "@ "<<sint.str()<<'\n';
				readid++;
				size+= 2+sint.str().size()+1;
				fasta_size += 2+sint.str().size()+1;
			}
			int dna_len = itSeq->item().getDataSize();
			char dna[dna_len+1];
			strncpy(dna, itSeq->item().getDataBuffer(), dna_len);
			dna[dna_len] = '\0';
			fout<<dna<<'\n';
			size+= dna_len+1;
			fasta_size += dna_len + 1;
			if( !isFasta)
			{
				string line = itSeq->item().getQuality();
				fout<<"+\n";
				fout<<line<<'\n';
				size+= 2+line.size()+1;
			}
			j++;
			itSeq->next();
		}
		//leon->fastq_block_size->push_back(size);
		//leon->orig_block_size->push_back(fasta_size);
		//leon->seq_per_block->push_back(j);
		fout.close();
		leon1->executeCompression(block_count, temp_file.c_str());
		string fasta_value, fastq_value;
		fasta_value = formatSize(fasta_size);
		fastq_value = formatSize(size);
                string filename = leon->_inputFilename+"_"+to_string(block_count)+".leon";
                if(lsetxattr(filename.c_str(), FASTA_SIZE_ATTR, fasta_value.c_str(), SIZE_LEN, 0)==-1)
                        StackFS_trace("Error in setxattr");
		if(lsetxattr(filename.c_str(), FASTQ_SIZE_ATTR, fastq_value.c_str(), SIZE_LEN, 0)==-1)
                        StackFS_trace("Error in setxattr");
		block_count++;
	}
	//leon->saveConfig();
}	

//Not Used Currently
static void compress_save(char * name){
	char * point = NULL;
	StackFS_trace("The file name to be saved: %s", name);
	if((point = strrchr(name,'.')) != NULL) {
		struct stat buffer;
		int compression_type = 1;
		string file_name (name);
		file_name = file_name+"_0" + EXT;
		if(stat(file_name.c_str(), &buffer) == 0) {
			struct file_pages * fp = (struct file_pages*) g_hash_table_lookup(file_table, name);
			struct file_pages * current = fp;
			if(compression_type == 1){
				//pthread_spin_lock(&spinlock);
				char* tmp_name = createFile(name, TMP_FILE);
				while(current!=NULL){
					if(current->dirty != 0){
						string temp (tmp_name);
						string name_check(name);
						if(name_check.find(".fq")!=string::npos 
								|| name_check.find(".fastq")!=string::npos){
							temp = temp+".fq";
						}else{
							temp = temp+".fa";
						}
						ofstream fout;
						fout.open(temp.c_str());
						fout<<current->file_page;
						fout.close();
						int count = 0;
						char** args = prep_args(name, true, false, true,count);
						splitCopy(count, args, temp.c_str(), current->block_id);
						current->dirty = 0;
					}
					current = current->next_page;
				}
				// pthread_spin_unlock(&spinlock);
			}else{
				int fd = open(file_name.c_str(), O_RDWR);
				int offset = 0;
				while(current!=NULL && current->next_page!=NULL){
					int off = 0;
					while(off<PG_SIZE){
						current->file_page[off] -= 1;
						off++;
					}
					offset+=pwrite(fd,current->file_page, PG_SIZE, offset);
					StackFS_trace("The flushed page: %s", current->file_page);
					current = current->next_page;
				}
				if(current!=NULL){
					int off = 0;
					int len = strlen(current->file_page);
					while(off<len){
						current->file_page[off] -= 1;
						off++;
					}
					offset+=pwrite(fd, current->file_page, len, offset);
					StackFS_trace("The flushed page: %s", current->file_page);
					current = current->next_page;
				} 
			}
		}
	}	
}

static void safe_remove(char* name){
	char* key = (char *)malloc( strlen(name)+1 );
	strcpy(key, name);
	StackFS_trace("File to be removed from cache: %s", key);
	if(g_hash_table_contains(file_table,key)){
		struct file_pages * fp = (struct file_pages*) g_hash_table_lookup(file_table, key);
		StackFS_trace("File found : %s removing %p", key, fp);
		struct file_pages * current = fp;
		while(current!=NULL && current->next_page!=NULL){
			current = current->next_page;
		}
		while(current!=NULL && current!=fp){
			current = current->prev_page;
			if(current->next_page!=NULL){
				free(current->next_page->file_page);
				free(current->next_page);
				current->next_page = NULL;
			}
		}
		delete fp->leon;
		if(fp!=NULL){
			free(fp->file_page);
			free(fp);
			fp = NULL;
		}
		if(g_hash_table_contains(file_cache,key)){
			struct file_node * fn = (struct file_node*) g_hash_table_lookup(file_cache, key);
                	fn->prev->next = fn->next;
			fn->next->prev = fn->prev;
                        free(fn);
			pthread_spin_lock(&spinlock);
			g_hash_table_remove(file_cache, key);
			pthread_spin_unlock(&spinlock);
                }
		pthread_spin_lock(&spinlock);
		g_hash_table_remove(file_table,key);
		pthread_spin_unlock(&spinlock);
	}
	if(key!=NULL){
		free(key);
		key = NULL;
	}
}

static void add_to_cache(char* name, int fd){
	char * key = (char*) malloc(strlen(name)+1);
	strcpy(key, name);
	bool cancelAdd = false;
	if(g_hash_table_contains(file_table,key)){
		struct file_pages * fp = (struct file_pages*) g_hash_table_lookup(file_table, key);
		struct file_pages * current = fp;
		long size = 0;
		while(current!=NULL){
			size += current->size;
			if(size>THRESHOLD){
				cancelAdd = true;
				safe_remove(name);
				break;
			}
                        current = current->next_page;
                }
	}
	if(!cancelAdd && g_hash_table_size(file_cache) >= CACHE_SIZE){
		while(g_hash_table_size (file_cache)>= CACHE_SIZE){
			if(file_node_end->prev != file_node_head){
				char* rem_name = file_node_end->prev->file_name;
				safe_remove(rem_name);
			}else
				break;
		}
	}
	if(!cancelAdd && g_hash_table_contains(file_table,key) && g_hash_table_size(file_cache)< CACHE_SIZE){
		StackFS_trace("Fd to be added to cache : %s", key);
		pthread_spin_lock(&spinlock);
		struct file_node* fn = (struct file_node*) malloc(sizeof(struct file_node));
		strcpy(fn->file_name, key);
		fn->prev = file_node_head;
		fn->next = file_node_head->next;
		file_node_head->next = fn;
		fn->next->prev = fn;
		g_hash_table_insert(file_cache, key , fn);
		pthread_spin_unlock(&spinlock);
	}
	if(cancelAdd){
		free(key);
		key = NULL;
	}
}

static void stackfs_ll_open(fuse_req_t req, fuse_ino_t ino,
		struct fuse_file_info *fi)
{
	int fd;
	struct stat buffer;
	generate_start_time(req);
	char * name = lo_name(req, ino);
	char * point = NULL;
	StackFS_trace("The file name: %s", name);
	string file_name (name);
        file_name = file_name+"_0" + EXT;

	char value[ATTR_SIZE];
	int ret = lgetxattr(file_name.c_str(), ATTR, value, ATTR_SIZE);
	string val(value);
	// It is assumed that the only files having size greater than 4kb will be set with 
	// the formats like FASTA, FASTQ and LEON, that serves as the additional check here. 
	// should this assumption change, the additional check has to be implemented here.
	if(ret > 0 && (val.find(FASTA)!=string::npos  || val.find(FASTQ)!=string::npos  ||
				val.find(LEON)!=string::npos )){
		string given_name(name);
		int format = 0;
		if(val.find(FASTA)!=string::npos){
			format = FASTA_F;
		}else if(val.find(FASTQ)!=string::npos){
			format = FASTQ_F;
		}	
		//We only need to decompress and save the fasta and fastq files
		if(format == FASTQ_F || format == FASTA_F) {
			fd = open(file_name.c_str(), fi->flags);
			decompress_store(name, fd, format);
		}
		else{
			//no need to decompress files with leon format
			fd = open(file_name.c_str(), fi->flags);
		}
	}else{
		// files with format none
		fd = open(name, fi->flags);
		if(stat(name, &buffer) ==0 && buffer.st_size == 0)
			created_files_map[name] = NONE_F;
	}
	generate_end_time(req);
	populate_time(req);

	StackFS_trace("Open called on name : %s and fuse inode : %llu kernel inode : %llu fd : %d", lo_name(req, ino), get_higher_fuse_inode_no(req, ino), get_lower_fuse_inode_no(req, ino), fd);
	StackFS_trace("Open name : %s and inode : %llu", lo_name(req, ino), get_lower_fuse_inode_no(req, ino));

	if (fd == -1)
		return (void) fuse_reply_err(req, errno);

	fi->fh = fd;

	fuse_reply_open(req, fi);
}

static void stackfs_ll_opendir(fuse_req_t req, fuse_ino_t ino,
		struct fuse_file_info *fi)
{
	DIR *dp;
	struct lo_dirptr *d;

	StackFS_trace("Opendir called on name : %s and inode : %llu",
			lo_name(req, ino), lo_inode(req, ino)->ino);
	generate_start_time(req);
	dp = opendir(lo_name(req, ino));
	generate_end_time(req);
	populate_time(req);

	if (dp == NULL)
		return (void) fuse_reply_err(req, errno);

	d = (struct lo_dirptr*)malloc(sizeof(struct lo_dirptr));
	d->dp = dp;
	d->offset = 0;
	d->entry = NULL;

	fi->fh = (uintptr_t) d;

	fuse_reply_open(req, fi);
}

static int findBlockId(unsigned long off, int &blockOff, vector<int> block_size){
        int blockId = 0;
	blockOff = 0;
	cout<<"Offset: "<<off<<endl;
	while(blockId < block_size.size() && off> block_size[blockId]){
		off -= block_size[blockId];
		blockId++;
	}
	blockOff = off;
	cout<<"block Id: "<<blockId<<" block Off: "<<blockOff<<endl;
	return blockId;
}


static void stackfs_ll_read(fuse_req_t req, fuse_ino_t ino, size_t size,
		off_t offset, struct fuse_file_info *fi)
{
	int res=0;
	(void) ino;
	//struct timespec start, end;
	//long time;
	//long time_sec;
	//unsigned long off_set= offset;
	StackFS_trace("StackFS Read start on inode : %llu", get_lower_fuse_inode_no(req, ino));
	if (USE_SPLICE) {
		struct fuse_bufvec buf = FUSE_BUFVEC_INIT(size);

		StackFS_trace("Splice Read name : %s, off : %lu, size : %zu",
				lo_name(req, ino), offset, size);

		generate_start_time(req);
		buf.buf[0].flags = fuse_buf_flags(FUSE_BUF_IS_FD | FUSE_BUF_FD_SEEK);
		buf.buf[0].fd = fi->fh;
		buf.buf[0].pos = offset;
		generate_end_time(req);
		populate_time(req);
		fuse_reply_data(req, &buf, FUSE_BUF_SPLICE_MOVE);
	} 
	else 
	{
		char *buf;

		StackFS_trace("Read on name : %s, Kernel inode : %llu, fuse inode : %llu, off : %lu, size : %zu",
				lo_name(req, ino), get_lower_fuse_inode_no(req, ino), get_higher_fuse_inode_no(req, ino), offset, size);
		buf = (char *)malloc(size);
		generate_start_time(req);
		char * name = lo_name(req, ino);
		int fd = 0;
		if(g_hash_table_contains(file_table,name))
		{
			struct file_pages* fp = (struct file_pages*)g_hash_table_lookup (file_table, name);
			struct file_pages* current = fp;
			struct file_pages* start = fp;
			struct file_pages* end = fp;
			struct file_pages* pre = fp;
			Leon* leon = fp->leon;
			bool block_found = false;
			vector<string>* out;
			int fromOff = 0, toOff = 0, bufOff = 0;
                        int seq_size = 0;
                        Leon* leonSeq = new Leon();
                        int count = 0;
			int last_file_size = -1;
                        char** targs;
			vector<int> block_sizes;
                        if(current->format == FASTA_F){
				getFileCount(name, block_sizes, true);
                        }else{
				getFileCount(name, block_sizes, false);
                        }
			cout<<"off_t: "<<offset<<endl;
			int fromBlock = findBlockId(offset, fromOff, block_sizes);
			int toBlock = findBlockId(offset+size, toOff, block_sizes);
			cout<<endl;
			StackFS_trace("The from block: %d to block %d", fromBlock, toBlock);
                        StackFS_trace("The fromOff: %d toOff %d", fromOff, toOff);
			if(fromBlock >= block_sizes.size()){
				fromBlock = block_sizes.size()-1;
				fromOff = block_sizes[fromBlock]-1;
			}
			if(fromOff>block_sizes[fromBlock])
				fromOff = block_sizes[fromBlock]-1;
			if(toBlock >= block_sizes.size()){
				toBlock = block_sizes.size()-1;
				toOff = block_sizes[toBlock]-1;
			}
			if(toOff>block_sizes[toBlock])
				toOff = block_sizes[toBlock]-1;
			StackFS_trace("The from block: %d to block %d", fromBlock, toBlock);
			StackFS_trace("The fromOff: %d toOff %d", fromOff, toOff);
			while(current->block_id < fromBlock){
				pre = current;
				current = current->next_page;
				if(current==NULL || current->block_id > fromBlock){
					Leon *leonRead = new Leon();
					count = 0;
					char** args;
					if(fp->format == FASTQ_F)
						args = prep_args(name, false, false, true,count);
					else 
						args = prep_args(name, false, true, true,count);
					leonRead->run(count, args);
                                        StackFS_trace("Decompressing New block id:%d", fromBlock);
					out = leonRead->executeDecompression(fromBlock);
					delete leonRead;
					struct file_pages* tmp = pre->next_page;
					current = pre;
					current->next_page = (struct file_pages*) malloc(sizeof(struct file_pages));
					current->next_page->prev_page = current;
					current->next_page->next_page = tmp;
					current->next_page->format = current->format;
                                        current->next_page->seq_size_A = current->seq_size_A;
					current->next_page->seq_size_Q = current->seq_size_Q;
					current->next_page->last_file_size_A = current->last_file_size_A;
					current->next_page->last_file_size_Q = current->last_file_size_Q;
					current = current->next_page;
					current->file_page = (char *) malloc ((*out)[0].size()+1);
					strcpy(current->file_page,(*out)[0].c_str());
					current->s_off = current->prev_page->e_off+1;
					current->e_off = current->s_off + (*out)[0].size();
					current->size = (*out)[0].size();
					current->block_id = fromBlock;
					current->leon = leon;
					current->isFastq = fp->isFastq;
					current->dirty = 0;
					current->sequenceAdded = false;
					delete out;
				}
			}
			start = current;
                        //Check if all the required blocks are present
			for(int i=fromBlock;i<=toBlock;i++)
			{
				if(current!=NULL && current->block_id == i){
					pre = current;
					current = current->next_page;
                                        StackFS_trace("Required block id:%d found in cache", i);
				}
				else
				{
					Leon *leonRead = new Leon();
					int count = 0;		
					char** args;
					if(fp->format == FASTQ_F)
						args = prep_args(name, false, false, true,count);
					else
						args = prep_args(name, false, true, true,count);
					leonRead->run(count, args);
                                        StackFS_trace("Decompressing New block id:%d", i);
					out = leonRead->executeDecompression(i);
					delete leonRead;
					struct file_pages* tmp = pre->next_page;
					current = pre;
					current->next_page = (struct file_pages*) malloc(sizeof(struct file_pages));
					current->next_page->prev_page = current;
					current->next_page->next_page = tmp;
					current->next_page->format = current->format;
                                        current->next_page->seq_size_A = current->seq_size_A;
                                        current->next_page->seq_size_Q = current->seq_size_Q;
                                        current->next_page->last_file_size_A = current->last_file_size_A;
                                        current->next_page->last_file_size_Q = current->last_file_size_Q;
					current = current->next_page;
					current->file_page = (char *) malloc ((*out)[0].size()+1);
					strcpy(current->file_page,(*out)[0].c_str());
					current->s_off = current->prev_page->e_off+1;
					current->e_off = current->s_off + (*out)[0].size();
					current->block_id = i;
					current->size = (*out)[0].size();
					current->leon = leon;
					current->isFastq = fp->isFastq;
					current->dirty = 0;
					current->sequenceAdded = false;
					pre = current;
					current = current->next_page;
					delete out;		
				}
			}
			if(start!=NULL){
				if(fromBlock == toBlock){
					if(fromOff < strlen(start->file_page)){
 						memcpy(buf, start->file_page+fromOff, toOff-fromOff);
					}
					res = toOff-fromOff;
					//StackFS_trace("Data: %s", buf);
				}else{
                                        //Get Everything from the from block
                                        if(block_sizes[fromBlock] <= strlen(start->file_page))
                                        	memcpy(buf, start->file_page+fromOff, block_sizes[fromBlock]-fromOff);
					else{
						if(fromOff < strlen(start->file_page)){
							 memcpy(buf, start->file_page+fromOff, 
								strlen(start->file_page)-fromOff);
						}
					}
					bufOff += block_sizes[fromBlock]-fromOff;
					start = start->next_page;
                                        // Get Everything in between from and to block.
					for(int i=fromBlock+1;i<toBlock && start!=NULL; i++){
 						int copy_size = min(block_sizes[i], (int)strlen(start->file_page));
						memcpy(buf+bufOff, start->file_page , copy_size); 
						bufOff+= block_sizes[i];
						start = start->next_page;
					}
                                        //Get Everything from the toBlock
					if(start!=NULL){
						if(size - bufOff < toOff){
							toOff = size - bufOff;
						}
						memcpy(buf+bufOff, start->file_page, toOff);
						bufOff += toOff;
					}
					//StackFS_trace("Data: %s", buf);
					res = bufOff;
				}
			}
		}
		else
		{	
			//clock_gettime(CLOCK_MONOTONIC, &start);
			StackFS_trace("The file name: %s not in map", name);
			res = pread(fi->fh, buf, size, offset);
		}
		generate_end_time(req);
		populate_time(req);
		if (res == -1)
			return (void) fuse_reply_err(req, errno);
		StackFS_trace("The buf size: %d, actual size: %d response:%d", strlen(buf), size, res);
		res = fuse_reply_buf(req, buf, res);
                if(buf!=NULL){
			free(buf);
			buf = NULL;
		}
	}
	StackFS_trace("StackFS Read end on inode : %llu", get_lower_fuse_inode_no(req, ino));
}

static void stackfs_ll_readdir(fuse_req_t req, fuse_ino_t ino, size_t size,
		off_t off, struct fuse_file_info *fi)
{
	struct lo_dirptr *d;
	char *buf = NULL;
	char *p = NULL;
	size_t rem;
	int err;
	(void) ino;

	StackFS_trace("Readdir called on name : %s and inode : %llu",
			lo_name(req, ino), lo_inode(req, ino)->ino);
	d = lo_dirptr(fi);
	buf = (char *)malloc(size*sizeof(char));
	unordered_set<string> file_list;
	unordered_map<string,  unordered_set<string>>::const_iterator got = 
		directory_list_map.find(lo_name(req, ino));
	if ( got != directory_list_map.end() )
		file_list = got->second;
	if (!buf)
		return (void) fuse_reply_err(req, ENOMEM);

	generate_start_time(req);
	/* If offset is not same, need to seek it */
	if (off != d->offset) {
		seekdir(d->dp, off);
		d->entry = NULL;
		d->offset = off;
	}
	p = buf;
	rem = size;
	while (1) {
		size_t entsize;
		off_t nextoff;

		if (!d->entry) {
			errno = 0;
			d->entry = readdir(d->dp);
			if (!d->entry) {
				if (errno && rem == size) {
					err = errno;
					goto error;
				}
				break;
			}
		}
		nextoff = telldir(d->dp);
		string cur_file(d->entry->d_name);
		string full_path(lo_name(req, ino));
		full_path += "/";
		full_path += cur_file;
		size_t dot = cur_file.find_last_of('.');
		char value[ATTR_SIZE];
        	int ret = lgetxattr(full_path.c_str(), ATTR, value, ATTR_SIZE);
        	string val(value);
		if( cur_file.substr(dot).compare(".qual")==0){
			string leonFile = full_path.substr(0, full_path.rfind("_"));
			leonFile += "_0.leon";
			ret = lgetxattr(leonFile.c_str(), ATTR, value, ATTR_SIZE);
			if(ret >= 0)
				val = string(value);
		}	
		if(dot != string::npos && cur_file.substr(dot).compare(EXT)==0 && val.find(LEON)==string::npos){
			if(cur_file.find("_") !=  string::npos){
				//cut the file name before the '_' in XXXX_12.leon or XXXX_0.leon
				size_t under = cur_file.rfind("_");
				string file_name  = cur_file.substr(0, under);
				if(file_list.find(file_name) == file_list.end()){
					struct stat st = {};
					st.st_ino = d->entry->d_ino;
					st.st_mode = d->entry->d_type << 12;
					entsize = fuse_add_direntry(req, p, rem,
							file_name.c_str(), &st, nextoff);
					if (entsize > rem)
						break;
					p += entsize;
					rem -= entsize;
					file_list.insert(file_name);
				}
			}
		}else if((dot == string::npos || cur_file.substr(dot).compare(".qual")!=0 ||
			 val.find(LEON)!=string::npos) && file_list.find(cur_file)==file_list.end()){
			struct stat st = {};
			st.st_ino = d->entry->d_ino;
			st.st_mode = d->entry->d_type << 12;
			entsize = fuse_add_direntry(req, p, rem,
					d->entry->d_name, &st, nextoff);
			/* The above function returns the size of the entry size even though
			 * the copy failed due to smaller buf size, so I'm checking after this
			 * function and breaking out incase we exceed the size.
			 */
			if (entsize > rem)
				break;

			p += entsize;
			rem -= entsize;
			file_list.insert(d->entry->d_name);
		}
		d->entry = NULL;
		d->offset = nextoff;
	}
	directory_list_map[lo_name(req, ino)] = file_list;
	generate_end_time(req);
	populate_time(req);
	fuse_reply_buf(req, buf, size - rem);
	free(buf);

	return;

error:
	generate_end_time(req);
	populate_time(req);
	free(buf);

	fuse_reply_err(req, err);
}

static void create_file(char* name){
	string fileName(name);
	if(created_files_map[name] == FASTQ_F)
	{
		int count = 0;
		char** args = prep_args(name,true, false, false,count);
		splitFiles(count, args, FASTQ_F);
		delete[] args;
		string newName(name);
		newName += ".fastq";
		unlink(newName.c_str());
	}
	else if(created_files_map[name] == FASTA_F)
	{
		int count = 0;
		char** args = prep_args(name,true, true, false,count);
		splitFiles(count, args, FASTA_F);
		delete[] args;
		string newName(name);
                newName += ".fasta";
		unlink(newName.c_str());
	}
	StackFS_trace("File created and compressed on name : %s ", name);
}

static void stackfs_ll_release(fuse_req_t req, fuse_ino_t ino,
		struct fuse_file_info *fi)
{
	(void) ino;

	StackFS_trace("Release called on name : %s and inode : %llu fd : %d ",
			lo_name(req, ino), lo_inode(req, ino)->ino, fi->fh);
	generate_start_time(req);
	char * name = lo_name(req, ino);
	string filename (name);
	char * point = NULL;
	unordered_map<string,int>::const_iterator got = created_files_map.find (name);
	if(g_hash_table_lookup(file_table, name)!=NULL && g_hash_table_lookup(file_cache, name)==NULL)
		add_to_cache(name, fi->fh);
	else if( got != created_files_map.end())
	{
		char* value = (char*)malloc(ATTR_SIZE);
		filename += "_0.leon";
		if(created_files_map[name] == FASTA_F){
			create_file(name);
			strcpy(value, FASTA);
			if(lsetxattr(filename.c_str(), ATTR, value, ATTR_SIZE, 0)==-1)
				StackFS_trace("Error in setxattr");
			lsetxattr(filename.c_str(), ORIG_ATTR, value, ATTR_SIZE, 0);
		}
		else if(created_files_map[name] == FASTQ_F){
			create_file(name);
			strcpy(value, FASTQ);
			if(lsetxattr(filename.c_str(), ATTR, value, ATTR_SIZE, 0)==-1)
                                StackFS_trace("Error in setxattr");
			lsetxattr(filename.c_str(), ORIG_ATTR, value, ATTR_SIZE, 0);
		}
		else if(created_files_map[name] == NONE_F){
                        try{
				StackFS_trace("File name %s is set to None", name);
				Leon* leon = new Leon();
				int format = 0;
				int count = 0;
				char** args;
				args = prep_args(name,false, false, true,count);
				leon->run(count, args);
				if(leon->checkLeonFile(name, format)){
					strcpy(value, LEON);
					if(lsetxattr(name, ATTR, value, ATTR_SIZE-1, 0)==-1)
						StackFS_trace("Error in setxattr");
					if(format == FASTA_F)
						strcpy(value, FASTA);
					else
						strcpy(value, FASTQ);
					lsetxattr(name, ORIG_ATTR, value, ATTR_SIZE-1, 0);
				}
				delete leon;
			}catch (...)  {
        			StackFS_trace("Exception Thrown Not a leon file");
    			}
		}
		if(value!=NULL){
			free(value);
		}
		created_files_map.erase(name);
	}
	close(fi->fh);
	generate_end_time(req);
	populate_time(req);

	fuse_reply_err(req, 0);
}

static void stackfs_ll_releasedir(fuse_req_t req, fuse_ino_t ino,
		struct fuse_file_info *fi)
{
	struct lo_dirptr *d;
	(void) ino;

	StackFS_trace("Releasedir called on name : %s and inode : %llu",
			lo_name(req, ino), lo_inode(req, ino)->ino);
	d = lo_dirptr(fi);
	directory_list_map.erase(lo_name(req, ino));
	generate_start_time(req);
	closedir(d->dp);
	generate_end_time(req);
	populate_time(req);
	free(d);
	fuse_reply_err(req, 0);
}


//Check for third line of every four lines(3, 7, 11, ...) to verify only a '+' is present
static bool detect_fastq(const char* buf, size_t size){
	stringstream ss(buf);
	string line;
	int lineno = 1;
	int linelen = 0;
	while(getline(ss, line)){
		/*if(lineno%4 == 2)
			linelen = line.size();
		*/if(lineno%4 == 3 && line.compare("+")!=0){
			return false;
		}
		/*if(lineno%4 == 0 && linelen!=line.size()){
			 StackFS_trace("linelen: %d and line size : %d and \nline = %s",
                        		linelen, line.size(), line.c_str());
			return false;
		}*/
		lineno++;
	}
	if(lineno<4)
		return false;
	return true;
}

//The first line should start from '>' or ';' and the first line cannot be an empty line,
//The line after the empty line should start with '>' or ';'. All the non empty line starting
//without the '<' or ';' should be less than 80 characters in length.
static bool detect_fasta(const char* buf, size_t size){
	stringstream ss(buf);
        string line;
        int lineno = 1;
        int linelen = 0;
	bool emptyline = false;
        while(getline(ss, line)){
		if(lineno ==1 && line.size()==0)
			return false;
		if(lineno ==1 &&  line.at(0)!='>' && line.at(0)!=';')
			return false;
		if(emptyline && line.size()==0)
			continue;
		if(emptyline && line.size()>0 && line.at(0)!='>' && line.at(0)!=';')
			return false;
	//	if(line.size()>0 && line.at(0)!='>' && line.at(0)!=';' && line.size() > 80)
	//		return false;
		if(line.size() == 0)
			emptyline = true;
		else if(line.size()>0)
			emptyline = false;
		lineno++;
	}
	if(lineno<2)
		return false;
	return true;
}	

static bool detect_leon(const char* buf, size_t size){
	if(buf[0]=='L')
		return true;
	else
               	return false;
}
static void stackfs_ll_write(fuse_req_t req, fuse_ino_t ino, const char *buf,
		size_t size, off_t off, struct fuse_file_info *fi)
{
	int res;
	(void) ino;
	generate_start_time(req);
	char * name = lo_name(req, ino);
	int compressedType = 2;
	bool appended = false;
	char* fname = (char*) malloc(PATH_MAX);
        strcpy(fname, name);
        strcat(fname, "_0.leon");
	char value[ATTR_SIZE];
        int ret = lgetxattr(fname, ATTR, value, ATTR_SIZE);
	string val(value);
	free(fname);
	fname = NULL;
	unordered_map<string,int>::const_iterator got = created_files_map.find (name);
	struct stat buffer;
	// To use the write operation, the file has to be either a newly created one or not a special 
	// file or less than 4096 bytes.
	if(got!= created_files_map.end()){
		if(size>=4096 && stat(name, &buffer) ==0 && buffer.st_size == 0){
			if(detect_fastq(buf, size))
				created_files_map[name] = FASTQ_F;
			else if(detect_fasta(buf, size))
				created_files_map[name] = FASTA_F;
			//else if(detect_leon(buf, size))
			//	created_files_map[name] = LEON_F;
			else
				created_files_map[name] = NONE_F;
			StackFS_trace("File detected to be : %d", created_files_map[name]);
		}
		res = pwrite(fi->fh, buf, size, off);
	}	
	else if(ret > 0 && (val.find(FASTA)!=string::npos ||val.find(FASTQ)!=string::npos  
				|| val.find(LEON)!=string::npos ))
	{
		res = -1;
                errno = -ENOTSUP;
	}else{
		res = pwrite(fi->fh, buf, size, off);
	}
	generate_end_time(req);
	populate_time(req);
	if (res == -1)
		return (void) fuse_reply_err(req, errno);
	fuse_reply_write(req, res);
}

//static void  stackfs_ll_lseek(fuse_req_t req, fuse_ino_t ino, off_t off,struct fuse_file_info *fi){
//	StackFS_trace("Seek on name : %s, off : %lu",
//                                lo_name(req, ino), off);
//}

#if	USE_SPLICE
static void stackfs_ll_write_buf(fuse_req_t req, fuse_ino_t ino,
		struct fuse_bufvec *buf, off_t off, struct fuse_file_info *fi)
{
	int res;
	(void) ino;

	struct fuse_bufvec dst = FUSE_BUFVEC_INIT(fuse_buf_size(buf));

	StackFS_trace("Splice Write_buf on name : %s, off : %lu, size : %zu",
			lo_name(req, ino), off, buf->buf[0].size);

	generate_start_time(req);
	dst.buf[0].flags = FUSE_BUF_IS_FD | FUSE_BUF_FD_SEEK;
	dst.buf[0].fd = fi->fh;
	dst.buf[0].pos = off;
	res = fuse_buf_copy(&dst, buf, FUSE_BUF_SPLICE_NONBLOCK);
	generate_end_time(req);
	populate_time(req);
	if (res >= 0)
		fuse_reply_write(req, res);
	else
		fuse_reply_err(req, res);
}
#endif

static int remove_entry(char* fullPath, bool isFastq){
	int res =0 ;
	string file_name (fullPath);
	string real_name (file_name + "_0" + EXT);
	string qual_name (file_name + "_0" + ".qual");
	int block_id = 0;
	struct stat buffer;
	while(stat(real_name.c_str(), &buffer) == 0) {
		res = unlink(real_name.c_str());
		if(isFastq && stat(qual_name.c_str(), &buffer) == 0)
			unlink(qual_name.c_str());
		block_id++;
		real_name = file_name + "_" + to_string(block_id) + EXT;
		qual_name = file_name + "_" + to_string(block_id) + ".qual";
	}
	int count = 0;
	/*Leon *leonRemove = new Leon();
	char** args = prep_args(fullPath, false, !isFastq, true,count);
	leonRemove->run(count, args);
	leonRemove->removeEntireFileConfig();;
	delete leonRemove;*/
	if(g_hash_table_contains(file_table,fullPath)){
		struct file_pages * fp = (struct file_pages*) g_hash_table_lookup(file_table, fullPath);
		StackFS_trace("File found : %s removing %p", fullPath, fp);
		struct file_pages * current = fp;
		while(current!=NULL && current->next_page!=NULL){
			current = current->next_page;
		}
		while(current!=NULL && current!=fp){
			current = current->prev_page;
			if(current->next_page!=NULL){
				free(current->next_page->file_page);
				free(current->next_page);
				current->next_page = NULL;
			}
		}
		delete fp->leon;
		if(fp!=NULL){
			free(fp->file_page);
			free(fp);
			fp = NULL;
		}
		pthread_spin_lock(&spinlock);
		g_hash_table_remove(file_table,fullPath);
		pthread_spin_unlock(&spinlock);
	}
	if(g_hash_table_contains(file_cache, fullPath)){
		struct file_node * fn = (struct file_node*) g_hash_table_lookup(file_cache, fullPath);
		fn->prev->next = fn->next;
		fn->next->prev = fn->prev;
		free(fn);
		pthread_spin_lock(&spinlock);
		g_hash_table_remove(file_cache, fullPath);
		pthread_spin_unlock(&spinlock);
	}
	return res;
}

static int remove_file(char* fullPath, char* value){
	int res = 0;
	string val(value);
	if(val.find(FASTQ)!=string::npos){
		res = remove_entry(fullPath , true);
	}else{
		res = remove_entry(fullPath, false);
	}
	return res;
}

static void stackfs_ll_unlink(fuse_req_t req, fuse_ino_t parent,
		const char *name)
{
	int res;
	char *fullPath = NULL;

	StackFS_trace("Unlink called on name : %s, parent inode : %llu",
			name, lo_inode(req, parent)->ino);
	fullPath = (char *)malloc(PATH_MAX);
	construct_full_path(req, parent, fullPath, name);
	generate_start_time(req);
	char * point = NULL;
	StackFS_trace("The file name: %s", name);
	string file_name(fullPath);
        file_name = file_name+"_0" + EXT;
        char value[ATTR_SIZE];
        int ret = lgetxattr(file_name.c_str(), ORIG_ATTR, value, ATTR_SIZE);
	string val(value);
        // It is assumed that the only files having size greater than 4kb will be set with
        // the formats like FASTA, FASTQ and LEON, that serves as the additional check here.
        // should this assumption change, the additional check has to be implemented here.
        if(ret > 0 && (val.find(FASTA)!=string::npos  || val.find(FASTQ)!=string::npos  ||
                                val.find(FASTA)!=string::npos )){
		if(val.find(FASTA)!=string::npos)
			res = remove_entry(fullPath, false);
		else
			res = remove_entry(fullPath, true);
	}else{
		res = unlink(fullPath);
	}
	generate_end_time(req);
	populate_time(req);
	if (res == -1)
		fuse_reply_err(req, errno);
	else
		fuse_reply_err(req, res);
	if (fullPath)
		free(fullPath);
}
//removes the fasta/fastq files present in the directory to be removed from map
static void remove_directory_from_map(char* fullPath, bool deleteFile){
	DIR           *d;
	struct dirent *dir;
	d = opendir(fullPath);
	if (d)
	{
		while ((dir = readdir(d)) != NULL)
		{
			string entry_name;
			char* point;
			entry_name = fullPath;
			entry_name += "/";
			entry_name +=  dir->d_name;
			StackFS_trace("%s\n", entry_name.c_str());
			char* temp_name = new char[entry_name.size()+1];
			strcpy(temp_name, entry_name.c_str());
			string visible_name = entry_name;
			if((point = strrchr(temp_name,'.')) != NULL && strcmp(point, ".leon")==0){
				if(entry_name.rfind("_")!=string::npos){
					visible_name = entry_name.substr(0,  entry_name.rfind("_"));
				}else{
					visible_name = entry_name;
				}
				StackFS_trace("File found : %s", visible_name.c_str());
			}else{
				visible_name = entry_name;
			}
			if(g_hash_table_contains(file_table, visible_name.c_str())){
				struct file_pages * fp = (struct file_pages*)
					g_hash_table_lookup(file_table, visible_name.c_str());
				StackFS_trace("File found : %s removing %p", visible_name.c_str(), fp);
				struct file_pages * current = fp;
				while(current!=NULL && current->next_page!=NULL){
					current = current->next_page;
				}
				while(current!=NULL && current!=fp){
					current = current->prev_page;
					if(current->next_page!=NULL){
						free(current->next_page->file_page);
						free(current->next_page);
						current->next_page = NULL;
					}
				}
				delete fp->leon;
				if(fp!=NULL){
					free(fp->file_page);
					free(fp);
					fp = NULL;
				}
				pthread_spin_lock(&spinlock);
				g_hash_table_remove(file_table, visible_name.c_str());
				pthread_spin_unlock(&spinlock);
			}
			if(g_hash_table_contains(file_cache, visible_name.c_str())){
				struct file_node * fn = (struct file_node*)
					g_hash_table_lookup(file_cache, visible_name.c_str());
				fn->prev->next = fn->next;
				fn->next->prev = fn->prev;
				free(fn);
				pthread_spin_lock(&spinlock);
				g_hash_table_remove(file_cache, visible_name.c_str());
				pthread_spin_unlock(&spinlock);
			}
			if(deleteFile)
				unlink(entry_name.c_str());
			delete(temp_name);
		}
		closedir(d);
	}						
}

static void stackfs_ll_rmdir(fuse_req_t req, fuse_ino_t parent,
		const char *name)
{
	int res;
	char *fullPath = NULL;

	StackFS_trace("rmdir called with name : %s, parent inode : %llu",
			name, lo_inode(req, parent)->ino);
	fullPath = (char *)malloc(PATH_MAX);
	construct_full_path(req, parent, fullPath, name);
	generate_start_time(req);
	DIR           *d;
	struct dirent *dir;
	d = opendir(fullPath);
	remove_directory_from_map(fullPath, true);
	res = rmdir(fullPath);
	generate_end_time(req);
	populate_time(req);
	if (res == -1)
		fuse_reply_err(req, errno);
	else
		fuse_reply_err(req, res);

	if (fullPath)
		free(fullPath);
}

static void forget_inode(fuse_req_t req, struct lo_inode *inode,
		uint64_t nlookup)
{
	int res;

	assert(inode->nlookup >= nlookup);
	inode->nlookup -= nlookup;

	if (!inode->nlookup)
		res = delete_from_hash_table(get_lo_data(req), inode);

	(void) res;
}

static void stackfs_ll_rename(fuse_req_t req, fuse_ino_t parent, const char* from, 
		fuse_ino_t new_parent, const char* to, unsigned int flags)
{
	int res = 0;
	struct stat statbuf;
	char * f_name = lo_name(req, parent);
	char* f_fullPath = (char *)malloc(PATH_MAX);
	construct_full_path(req, parent, f_fullPath, from);
	char * t_name = lo_name(req, new_parent);
	char* t_fullPath = (char *)malloc(PATH_MAX);
	construct_full_path(req, new_parent, t_fullPath, to);
	generate_start_time(req);
	StackFS_trace("Rename : from: %s, to: %s, from Name:%s, to_name:%s",
			f_name, t_name, from, to);
	char * point = NULL;
	char* point2 = NULL;
	char value[ATTR_SIZE];
	string f_file_name (f_fullPath);
        string t_file_name (t_fullPath);
	string f_real_name (f_file_name + "_0" + EXT);
    	string t_real_name (t_file_name + "_0" + EXT);
        string f_qual_name (f_file_name + "_0" + ".qual");
        string t_qual_name (t_file_name + "_0" + ".qual");
        int ret = lgetxattr(f_real_name.c_str(), ORIG_ATTR, value, ATTR_SIZE);
	string val(value);
	if(ret > 0 && (val.find(FASTA)!=string::npos  || val.find(FASTQ)!=string::npos  ||
				val.find(LEON)!=string::npos )){
		bool isFastq = false;
		int block_id = 0;
		struct stat buffer;
		if(val.find(FASTQ)!=string::npos)
			isFastq = true;
		while(stat(f_real_name.c_str(), &buffer) == 0) {
			res = rename(f_real_name.c_str(), t_real_name.c_str());
			if(isFastq && stat(f_qual_name.c_str(), &buffer) == 0){
				res = rename(f_qual_name.c_str(), t_qual_name.c_str());
			}
			block_id++;
			f_real_name = f_file_name + "_" + to_string(block_id) + EXT;
			t_real_name = t_file_name + "_" + to_string(block_id) + EXT;
			f_qual_name = f_file_name + "_" + to_string(block_id) + ".qual";
			t_qual_name = t_file_name + "_" + to_string(block_id) + ".qual";
		}
		int count = 0;
		/*Leon *leonRemove = new Leon();
		char** args = prep_args(f_fullPath, false, !isFastq, true,count);
		leonRemove->run(count, args);
		leonRemove->renameConfig(to, t_fullPath);
		delete leonRemove;*/
		if(g_hash_table_contains(file_table,f_fullPath)){
			struct file_pages * fp = (struct file_pages*) 
				g_hash_table_lookup(file_table, f_fullPath);
			StackFS_trace("File found : %s removing %p", f_fullPath, fp);
			struct file_pages * current = fp;
			while(current!=NULL && current->next_page!=NULL){
				current = current->next_page;
			}
			while(current!=NULL && current!=fp){
				current = current->prev_page;
				if(current->next_page!=NULL){
					free(current->next_page->file_page);
					free(current->next_page);
					current->next_page = NULL;
				}
			}
			delete fp->leon;
			if(fp!=NULL){
				free(fp->file_page);
				free(fp);
				fp = NULL;
			}
			pthread_spin_lock(&spinlock);
			g_hash_table_remove(file_table,f_fullPath);
			pthread_spin_unlock(&spinlock);
		}
		if(g_hash_table_contains(file_cache, f_fullPath)){
			struct file_node * fn = (struct file_node*) 
				g_hash_table_lookup(file_cache, f_fullPath);
			fn->prev->next = fn->next;
			fn->next->prev = fn->prev;
			free(fn);
			pthread_spin_lock(&spinlock);
			g_hash_table_remove(file_cache, f_fullPath);
			pthread_spin_unlock(&spinlock);
		}
	}else{
		res = rename(f_fullPath, t_fullPath);	
	}
	generate_end_time(req);
	populate_time(req);
	if (res == -1)
		fuse_reply_err(req, errno);
	else
		fuse_reply_err(req, res);
}	

static void stackfs_ll_forget(fuse_req_t req, fuse_ino_t ino, uint64_t nlookup)
{
	struct lo_inode *inode = lo_inode(req, ino);

	generate_start_time(req);
	StackFS_trace("Forget name : %s, inode : %llu and lookup count : %llu",
			inode->name, inode->ino, nlookup);
	forget_inode(req, inode, nlookup);
	generate_end_time(req);
	populate_time(req);

	fuse_reply_none(req);
}

static void stackfs_ll_forget_multi(fuse_req_t req, size_t count,
		struct fuse_forget_data *forgets)
{
	size_t i;
	struct lo_inode *inode;
	fuse_ino_t ino;
	uint64_t nlookup;

	generate_start_time(req);
	StackFS_trace("Batch Forget count : %zu", count);
	for (i = 0; i < count; i++) {
		ino = forgets[i].ino;
		nlookup = forgets[i].nlookup;
		inode = lo_inode(req, ino);

		//StackFS_trace("Forget %zu name : %s, lookup count : %llu",
		//				i, inode->name, nlookup);
		forget_inode(req, inode, nlookup);
	}
	generate_end_time(req);
	populate_time(req);

	fuse_reply_none(req);
}

static void stackfs_ll_flush(fuse_req_t req, fuse_ino_t ino,
		struct fuse_file_info *fi)
{
	int err;

	StackFS_trace("Flush called on name : %s and inode : %llu",
			lo_name(req, ino), lo_inode(req, ino)->ino);
	generate_start_time(req);
	err = 0;
	generate_end_time(req);
	populate_time(req);
	fuse_reply_err(req, err);
}

static void stackfs_ll_statfs(fuse_req_t req, fuse_ino_t ino)
{
	int res;
	struct statvfs buf;

	if (ino) {
		StackFS_trace("Statfs called with name : %s, and inode : %llu",
				lo_name(req, ino), lo_inode(req, ino)->ino);
		memset(&buf, 0, sizeof(buf));
		generate_start_time(req);
		res = statvfs(lo_name(req, ino), &buf);
		generate_end_time(req);
		populate_time(req);
		StackFS_trace("Statfs block size: %d and block count %d",
				buf.f_bsize, buf.f_bsize);
	}

	if (!res)
		fuse_reply_statfs(req, &buf);
	else
		fuse_reply_err(req, res);
}

static void stackfs_ll_fsync(fuse_req_t req, fuse_ino_t ino, int datasync,
		struct fuse_file_info *fi)
{
	int res;

	StackFS_trace("Fsync on name : %s, inode : %llu, datasync : %d",
			lo_name(req, ino), lo_inode(req, ino)->ino, datasync);
	generate_start_time(req);
	if (datasync)
		res = fdatasync(fi->fh);
	else
		res = fsync(fi->fh);
	generate_end_time(req);
	populate_time(req);

	fuse_reply_err(req, res);
}

static void stackfs_ll_getxattr(fuse_req_t req, fuse_ino_t ino,
		const char *name, size_t size)
{
	int res;
	struct stat buf;
	char* fullPath = lo_name(req, ino);
	//StackFS_trace("Function Trace : Getxattr %s", fullPath);
	if (size) {
		char *value = (char *) malloc(size);
		generate_start_time(req);
		char* fname = (char*) malloc(PATH_MAX);
	        strcpy(fname, fullPath);
        	strcat(fname, "_0.leon");
		if(stat(fname, &buf)==0)
			res = lgetxattr(fname, name, value, size);
		else
			res = lgetxattr(lo_name(req, ino), name, value, size);
		generate_end_time(req);
		populate_time(req);
		if (res > 0)
			fuse_reply_buf(req, value, res);
		else
			fuse_reply_err(req, errno);

		free(value);
	} else {
		generate_start_time(req);
		char* fname = (char*) malloc(PATH_MAX);
                strcpy(fname, fullPath);
                strcat(fname, "_0.leon");
                if(stat(fname, &buf)==0)
                        res = lgetxattr(fname, name, NULL, 0);
                else
			res = lgetxattr(lo_name(req, ino), name, NULL, 0);
		generate_end_time(req);
		populate_time(req);
		if (res >= 0)
			fuse_reply_xattr(req, res);
		else
			fuse_reply_err(req, errno);
	}
}

static void stackfs_ll_setxattr(fuse_req_t req, fuse_ino_t ino,
                const char *name, const char* value, size_t size, int flags)
{
        int res;
	char* fullPath = lo_name(req, ino);
	struct stat buf;
	generate_start_time(req);
	char* fname = (char*) malloc(PATH_MAX);
       	strcpy(fname, fullPath);
        strcat(fname, "_0.leon");
	StackFS_trace("Function Trace : Setxattr %s", fname);
       	if(stat(fname, &buf)==0){
		res = -1;
		errno = ENOTSUP;
		if(strcmp(name, ATTR)==0){
       	 		char vl[ATTR_SIZE];
       	 		int ret = lgetxattr(fname, ATTR, vl, ATTR_SIZE);
			char orig_val[ATTR_SIZE];
			int orig_ret = lgetxattr(fname, ORIG_ATTR, orig_val, ATTR_SIZE);
			StackFS_trace("Function Trace : 1Setxattr %s, %s, %s, %d\n", orig_val, vl, value, size);
			string val(vl);
			string orig_vl(orig_val);
			string newVal(value);
			//possible transitions FASTA->LEON, FASTA->FASTQ(if the original file is 
			//FASTQ, LEON->FASTA, LEON->FASTQ(if the original file is fastq, FASTQ->FASTA,
			//FASTQ->LEON.
			if(ret > 0 && (val.find(FASTA)!=string::npos  || val.find(FASTQ)!=string::npos ||
					val.find(LEON)!=string::npos )) {
				if(val.find(FASTA)!=string::npos && ((size ==4 && strncmp(value, LEON, size)==0) ||
						(size == 5 && strncmp(value, FASTA, size)==0) || 
					(size == 5 && orig_vl.find(FASTQ)!=string::npos  && 
					strncmp(value, FASTQ, size)==0))){
					res = lsetxattr(fname, name, value, size, flags);
				}else if(val.find(FASTQ)!=string::npos  && ((size ==4 && strncmp(value, LEON, size)==0) ||
                                                (size ==5 && strncmp(value, FASTA, size)==0) || 
						(size==5 && strncmp(value, FASTQ, size)==0))){	
        				res = lsetxattr(fname, name, value, size, flags);
					StackFS_trace("Function Trace : Setxattr %s, %s, %s", orig_val, vl, value);
				}else if(val.find(LEON)!=string::npos  && ((size==4 && strncmp(value, LEON, 4)==0) ||
                                                (size==5 && strncmp(value, FASTA, size)==0) || 
					(size == 5 && orig_vl.find(FASTQ)!=string::npos  
							&& strncmp(value, FASTQ, size)==0))){
					res = lsetxattr(fname, name, value, size, flags);
				}
			}
			if(res!=-1){
				if(g_hash_table_contains(file_table,fullPath)){
					struct file_pages * fp = (struct file_pages*) 
									g_hash_table_lookup(file_table, fullPath);
					StackFS_trace("File found : %s removing %p", fullPath, fp);
					struct file_pages * current = fp;
					while(current!=NULL && current->next_page!=NULL){
						current = current->next_page;
					}
					while(current!=NULL && current!=fp){
						current = current->prev_page;
						if(current->next_page!=NULL){
							free(current->next_page->file_page);
							free(current->next_page);
							current->next_page = NULL;
						}
					}
					delete fp->leon;
					if(fp!=NULL){
						free(fp->file_page);
						free(fp);
						fp = NULL;
					}
					pthread_spin_lock(&spinlock);
					g_hash_table_remove(file_table,fullPath);
					pthread_spin_unlock(&spinlock);
				}
				if(g_hash_table_contains(file_cache, fullPath)){
					struct file_node * fn = (struct file_node*) g_hash_table_lookup(file_cache, fullPath);
					fn->prev->next = fn->next;
					fn->next->prev = fn->prev;
					free(fn);
					pthread_spin_lock(&spinlock);
					g_hash_table_remove(file_cache, fullPath);
					pthread_spin_unlock(&spinlock);
				}		
			}
		}else if( strcmp(name, ORIG_ATTR)!=0 && strcmp(name, FASTA_SIZE_ATTR)!=0 
							&& strcmp(name, FASTQ_SIZE_ATTR)!=0){
			res = lsetxattr(fname, name, value, size, flags);
		}
	}else{
		if(strcmp(name, ATTR)==0 || strcmp(name, ORIG_ATTR)==0 
			|| strcmp(name, FASTA_SIZE_ATTR)!=0 || strcmp(name, FASTQ_SIZE_ATTR)!=0 ){
			errno = ENOTSUP;
			res = -1;
		}else 
			res = lsetxattr(lo_name(req, ino), name, value, size, flags);
	}
	generate_end_time(req);
	populate_time(req);
	free(fname);
	fname = NULL;
	if (res >= 0)
		fuse_reply_buf(req, value, res);
	else
		fuse_reply_err(req, errno);
}

static void stackfs_ll_listxattr(fuse_req_t req, fuse_ino_t ino, size_t size)
{
	int res;
        struct stat buf;
        char* fullPath = lo_name(req, ino);
        StackFS_trace("Function Trace : Getxattr %s", fullPath);
        if (size) {
                char *value = (char *) malloc(size);
                generate_start_time(req);
                char* fname = (char*) malloc(PATH_MAX);
                strcpy(fname, fullPath);
                strcat(fname, "_0.leon");
                if(stat(fname, &buf)==0)
                        res = llistxattr(fname, value, size);
                else
                        res = llistxattr(lo_name(req, ino), value, size);
                generate_end_time(req);
                populate_time(req);
                if (res > 0)
                        fuse_reply_buf(req, value, res);
                else
                        fuse_reply_err(req, errno);

                free(value);
        } else {
                generate_start_time(req);
                char* fname = (char*) malloc(PATH_MAX);
                strcpy(fname, fullPath);
                strcat(fname, "_0.leon");
                if(stat(fname, &buf)==0)
                        res = llistxattr(fname, NULL, 0);
                else
                        res = llistxattr(lo_name(req, ino), NULL, 0);
                generate_end_time(req);
		populate_time(req);
		if (res >= 0)
			fuse_reply_xattr(req, res);
		else
			fuse_reply_err(req, errno);
	}
}

static struct fuse_lowlevel_ops hello_ll_oper = {};
/*	lookup	:stackfs_ll_lookup,
forget  :stackfs_ll_forget,
getattr	:stackfs_ll_getattr,
setattr :stackfs_ll_setattr,
mkdir   :stackfs_ll_mkdir,
unlink  :stackfs_ll_unlink,
rmdleName.c_str());ir   :stackfs_ll_rmdir,
open    :stackfs_ll_open,
read    :stackfs_ll_read,
write   :stackfs_ll_write,
flush   :stackfs_ll_flush,
release :stackfs_ll_release,
fsync   :stackfs_ll_fsync,
opendir :stackfs_ll_opendir,
readdir :stackfs_ll_readdir,
releasedir:stackfs_ll_releasedir,
statfs	:stackfs_ll_statfs,
#if	TESTING_XATTR
getxattr:	stackfs_ll_getxattr,
#endif
create  :       stackfs_ll_create,
#if	USE_SPLICE
write_buf:	stackfs_ll_write_buf,
#endif
forget_multi:   stackfs_ll_forget_multi
};*/

struct stackFS_info {
	char	*rootDir;
	char	*statsDir;/* Path to copy any statistics details */
	double	attr_valid;/* Time in secs for attribute validation */
	int	is_help;
	int	tracing;
};

#define STACKFS_OPT(t, p) { t, offsetof(struct stackFS_info, p), 1 }

static const struct fuse_opt stackfs_opts[] = {
	STACKFS_OPT("-r %s", rootDir),
	STACKFS_OPT("--rootdir=%s", rootDir),
	STACKFS_OPT("--statsdir=%s", statsDir),
	STACKFS_OPT("--attrval=%lf", attr_valid),
	FUSE_OPT_KEY("--tracing", 1),
	FUSE_OPT_KEY("-h", 0),
	FUSE_OPT_KEY("--help", 0),
	FUSE_OPT_END
};

static int stackfs_process_arg(void *data, const char *arg,
		int key, struct fuse_args *outargs)
{
	struct stackFS_info *s_info =  (struct stackFS_info *)data;

	(void)outargs;
	(void)arg;

	switch (key) {
		case 0:
			s_info->is_help = 1;
			return 0;
		case 1:
			s_info->tracing	= 1;
			return 0;
		default:
			return 1;
	}
}

int main(int argc, char **argv)
{
	int res = 0, err = 0;
	char *rootDir = NULL;
	char *statsDir = NULL;
	char *resolved_statsDir = NULL;
	char *resolved_rootdir_path = NULL;
	int multithreaded;
	hello_ll_oper.lookup  = stackfs_ll_lookup;
	hello_ll_oper.forget  = stackfs_ll_forget;
	hello_ll_oper.getattr = stackfs_ll_getattr;
	hello_ll_oper.setattr = stackfs_ll_setattr;
	hello_ll_oper.mkdir   = stackfs_ll_mkdir;
	hello_ll_oper.unlink  = stackfs_ll_unlink;
	hello_ll_oper.rmdir   = stackfs_ll_rmdir;
	hello_ll_oper.rename  = stackfs_ll_rename,
	hello_ll_oper.open    = stackfs_ll_open;
	hello_ll_oper.read    = stackfs_ll_read;
	hello_ll_oper.write   = stackfs_ll_write;
	hello_ll_oper.flush   = stackfs_ll_flush;
	hello_ll_oper.release = stackfs_ll_release;
	hello_ll_oper.fsync   = stackfs_ll_fsync;
	hello_ll_oper.opendir = stackfs_ll_opendir;
	hello_ll_oper.readdir = stackfs_ll_readdir;
	hello_ll_oper.releasedir = stackfs_ll_releasedir;
	hello_ll_oper.statfs  = stackfs_ll_statfs;
//#if     TESTING_XATTR
	hello_ll_oper.getxattr = stackfs_ll_getxattr;
	hello_ll_oper.setxattr = stackfs_ll_setxattr;
	hello_ll_oper.listxattr = stackfs_ll_listxattr;	
//#endif
	hello_ll_oper.create = stackfs_ll_create;
#if     USE_SPLICE
	hello_ll_oper.write_buf = stackfs_ll_write_buf;
#endif
	hello_ll_oper.forget_multi = stackfs_ll_forget_multi;
	file_table = g_hash_table_new(g_str_hash, g_str_equal);
	file_cache = g_hash_table_new(g_str_hash, g_str_equal);
	file_node_head = (struct file_node*) malloc(sizeof(struct file_node));
	file_node_head->prev = NULL;
	file_node_end = (struct file_node*) malloc(sizeof(struct file_node));
	file_node_end->next = NULL;
	file_node_head->next = file_node_end;
	file_node_end->prev = file_node_head;
	//file_queue = g_queue_new();
	struct fuse_args args = FUSE_ARGS_INIT(argc, argv);
	/*Default attr valid time is 1 sec*/
	struct stackFS_info s_info = {NULL, NULL, 1.0, 0, 0};
	struct lo_data *lo = NULL;
	res = fuse_opt_parse(&args, &s_info, stackfs_opts, stackfs_process_arg);

	if (res) {
		printf("Failed to parse arguments\n");
		return -1;
	}

	if (s_info.is_help) {
		print_usage();
		return 0;
	}

	if (!s_info.rootDir) {
		printf("Root Directory is mandatory\n");
		print_usage();
		return -1;
	}

	if (s_info.statsDir) {
		statsDir = s_info.statsDir;
		resolved_statsDir = realpath(statsDir, NULL);
		if (resolved_statsDir == NULL) {
			printf("There is a problem in resolving the stats ");
			printf("Directory passed %s\n", statsDir);
			perror("Error");
			res = -1;
			goto out1;
		}
	}

	rootDir = s_info.rootDir;

	if (rootDir) {
		lo = (struct lo_data *) calloc(1, sizeof(struct lo_data));
		if (!lo) {
			fprintf(stderr, "fuse: memory allocation failed\n");
			res = -1;
			goto out2; /* free the resolved_statsDir */
		}
		resolved_rootdir_path = realpath(rootDir, NULL);
		if (!resolved_rootdir_path) {
			printf("There is a problem in resolving the root ");
			printf("Directory Passed %s\n", rootDir);
			perror("Error");
			res = -1;
			goto out3; /* free both resolved_statsDir, lo */
		}
		if (res == 0) {
			(lo->root).name = resolved_rootdir_path;
			(lo->root).ino = FUSE_ROOT_ID;
			(lo->root).nlookup = 2;
			(lo->root).next = (lo->root).prev = NULL;
			lo->attr_valid = s_info.attr_valid;
			/* Initialise the hash table and assign */
			res = hash_table_init(&lo->hash_table);
			if (res == -1)
				goto out4;
			/* Initialise the spin lock for table */
			pthread_spin_init(&(lo->spinlock), 0);
		}
	} else {
		res = -1;
		goto out2;
	}

	struct fuse_chan *ch;
	char *mountpoint;

	res = fuse_parse_cmdline(&args, &mountpoint, &multithreaded, NULL);

	/* Initialise the spinlock before the logfile creation */
	pthread_spin_init(&spinlock, 0);

	//if (s_info.tracing) {
	err = log_open(resolved_statsDir);
	if (err)
		printf("No log file created(but not a fatle error, ");
	printf("so proceeding)\n");
	//} else
	//	printf("No tracing\n");

	printf("Multi Threaded : %d\n", multithreaded);

	if (res != -1) {
		ch = fuse_mount(mountpoint, &args);
		if (ch) {
			struct fuse_session *se;

			printf("Mounted Successfully\n");
			se = fuse_lowlevel_new(&args, &hello_ll_oper,
					sizeof(hello_ll_oper), lo);
			if (se) {
				if (fuse_set_signal_handlers(se) != -1) {
					fuse_session_add_chan(se, ch);
					if (resolved_statsDir)
						fuse_session_add_statsDir(se,
								resolved_statsDir);
					//if (multithreaded)
					//	err = fuse_session_loop_mt(se);
					//else
						err = fuse_session_loop(se);
					(void) err;

					fuse_remove_signal_handlers(se);
					fuse_session_remove_statsDir(se);
					fuse_session_remove_chan(ch);
				}
				fuse_session_destroy(se);
			}
			StackFS_trace("Function Trace : Unmount");
			fuse_unmount(mountpoint, ch);
		}
	}

	/* free the arguments */
	fuse_opt_free_args(&args);

	/* destroy the lock protecting the hash table */
	pthread_spin_destroy(&(lo->spinlock));

	/* free up the hash table */
	free_hash_table(lo);

	/* destroy the hash table */
	hash_table_destroy(&lo->hash_table);

	/* destroy the lock protecting the log file */
	pthread_spin_destroy(&spinlock);

	/* close the log file (if any) */
	log_close();

out4:
	if (resolved_rootdir_path)
		free(resolved_rootdir_path);

out3:
	if (lo)
		free(lo);

out2:
	if (resolved_statsDir)
		free(resolved_statsDir);

out1:
	if(file_node_head !=NULL)
		free(file_node_head);
	if(file_node_end !=NULL)
		free(file_node_end);
	return res;

}
