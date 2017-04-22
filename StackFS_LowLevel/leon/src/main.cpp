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
#define MAX_SIZE 3145728
#define CACHE_SIZE 1
#define EXT ".leon"
#define CONFIG "/.biofs.dirinfo"
#define TMP_FILE "/.tmp_file"
#define CONFIG_LEN 15
#define TXT_EXT_LEN 4
#define TXT_EXT ".txt"
FILE *logfile;
#define TESTING_XATTR 0
#define USE_SPLICE 0
#define PG_SIZE 4096
#define TRACE_FILE "/trace_stackfs.log"
#define TRACE_FILE_LEN 18
pthread_spinlock_t spinlock; /* Protecting the above spin lock */
char banner[4096];
GHashTable* file_table;
GHashTable* file_cache;
struct file_node{
	char file_name[PATH_MAX];
	Leon *leon;
	int total_size;
	bool isFastq;
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
        int dirty;
        bool sequenceAdded;
	struct file_pages* prev_page;
        struct file_pages* next_page;
};


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

void log_close(void)
{

	if (logfile)
		fclose(logfile);
}

static char** prep_args(char * name, bool compress, bool isFastq, int blockId){
        char** argv = new char*[7];
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
	strcpy(argv[4], "-b");
	argv[5] = new char[33];
	snprintf(argv[5], sizeof(argv[5]), "%d", blockId);
	if(isFastq){
		argv[6] = new char[8];
		strcpy(argv[6],"-noqual");
	}
        /*argv[4] = new char[10];
        strcpy(argv[4], "-nb-cores");
        argv[5] = new char[2];
        strcpy(argv[5], "1");*/
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
	time = print_timer();
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
        StackFS_trace("The file name: %s is looked up", fullPath);
        if((point = strrchr(fullPath,'.')) != NULL && (strcmp(point, ".fasta")==0||strcmp(point, ".fastq")==0
                        ||strcmp(point, ".fa")==0||strcmp(point, ".fq")==0)) {
                struct stat buffer;
                bool isFastq = true;
		int len = strlen(fullPath)+EXT_LEN + 3;
		file_name = (char*) malloc(PATH_MAX);
		strcpy(file_name, fullPath);
		strcat(file_name, "_0");
		strcat(file_name, EXT);
		//given_name = given_name+"_0" + EXT;
		//given
		file_name[len-1] = '\0';
                /*int file_name_len = strlen(fullPath)-strlen(point);
                file_name = (char *) malloc(file_name_len+EXT_LEN+1);
                strncpy(file_name, fullPath, file_name_len);
                strncpy(file_name+file_name_len, EXT, EXT_LEN);
                file_name[file_name_len+EXT_LEN] = '\0';
                StackFS_trace("The file name: %s", file_name);*/
                if(stat(file_name, &buffer) == 0) {
                        res = stat(file_name, &e.attr);
			Leon* leon = new Leon();
                        leon->run(6, prep_args(fullPath,false, false, 0));
			/*char * config = createFile(fullPath, CONFIG);
	                string input(fullPath);
                	size_t file_base = input.find_last_of('/');
                	input = input.substr(file_base);
                	char* test = strdup(input.c_str());*/
			char * config = createFile(file_name, CONFIG);
                        char * point = strrchr(file_name, '/');
                        char * input = (char *) malloc(PATH_MAX);
                        strncpy(input, point, strlen(point)-EXT_LEN);
                	e.attr.st_size = leon->getFileSize(config, input);
                	free(config);
                	free(input);
			/*
			char * config = createFile(file_name, CONFIG);
			char * point = strrchr(file_name, '/');
			char * input = (char *) malloc(strlen(point)-EXT_LEN+1);
			strncpy(input, point, strlen(point)-EXT_LEN);
                        e.attr.st_size = leon->getFileSize(config, input);
                        StackFS_trace("Statfs block size: %d", e.attr.st_size);
			free(config);
			free(input);*/
                        delete leon;
		}else{
			res = stat(fullPath, &e.attr);
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
	char * point = NULL;
        StackFS_trace("The file name: %s is looked up", name);
        if((point = strrchr(name,'.')) != NULL && (strcmp(point, ".fasta")==0||strcmp(point, ".fastq")==0
                        ||strcmp(point, ".fa")==0||strcmp(point, ".fq")==0)) {
                struct stat buffer;
                string given_name(name);
                bool isFastq = true;
		int len = strlen(name)+EXT_LEN + 3;
                file_name = (char*) malloc(PATH_MAX);
                strcpy(file_name, name);
                strcat(file_name, "_0");
                strcat(file_name, EXT);
                //given_name = given_name+"_0" + EXT;
                //given
                file_name[len-1] = '\0';
		if(stat(file_name, &buffer) == 0) {
                        res = stat(file_name, &buf);
                        Leon* leon = new Leon();
                        leon->run(6, prep_args(name,false, false, 0));
			char * config = createFile(file_name, CONFIG);
                        char * point = strrchr(file_name, '/');
                        char * input = (char *) malloc(PATH_MAX);
                        strncpy(input, point, strlen(point)-EXT_LEN);
                        /*char * config = createFile(name, CONFIG);
                        string input(name);
                        size_t file_base = input.find_last_of('/');
                        input = input.substr(file_base);
                        char* test = strdup(input.c_str());*/
                        buf.st_size = leon->getFileSize(config, input);
                        free(config);
                        free(input);
                        delete leon;
                /*
                int file_name_len = strlen(name)-strlen(point);
                char * file_name = (char *) malloc(file_name_len+EXT_LEN+1);
                strncpy(file_name, name, file_name_len);
                strncpy(file_name+file_name_len, EXT, EXT_LEN);
                file_name[file_name_len+EXT_LEN] = '\0';
                if(stat(file_name, &buffer) == 0 && strcmp(point,EXT) != 0) {
                        res = stat(file_name, &buf);
			StackFS_trace("The file name passed: %s", file_name);
			Leon* leon = new Leon();
                	leon->run(6, prep_args(file_name,false, false, 0));
			char * config = createFile(file_name, CONFIG);
			char * point = strrchr(file_name, '/');
                        char * input = (char *) malloc(strlen(point)-EXT_LEN+1);
                        strncpy(input, point, strlen(point)-EXT_LEN);
                        buf.st_size = leon->getFileSize(config, input);
			StackFS_trace("Statfs block size: %d", buf.st_size);
			free(config);
			free(input);
			delete leon;*/
                }else{
                        res = stat(name, &buf);
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
static void decompress_store(char * name, int fd){
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
		leon->run(6, prep_args(name,false, false, 0));
		vector<string>* out =leon->executeDecompression(0);
		current->file_page = (char *) malloc ((*out)[0].size()+1);
		//strcpy(current->file_page, out);
		strcpy(current->file_page,(*out)[0].c_str());
		offSet = strlen(current->file_page)+1;
		current->s_off = 0;
		current->e_off = (*out)[0].size();
		current->block_id = 0;
		current->leon = leon;
		current->dirty = 0;
		current->sequenceAdded = false;
		//current->isFastq = isFastq;
		char * config = createFile(name, CONFIG);
		string input(name);
		size_t file_base = input.find_last_of('/');
		input = input.substr(file_base);
                //char * point = strrchr(file_name.c_str(), '/');
                //char * input = (char *) malloc(strlen(point)-EXT_LEN+1);
                //strncpy(input, point, strlen(point)-EXT_LEN);
		char* test = strdup(input.c_str());
                leon->readConfig(config, test);
                free(config);
		free(test);
		//StackFS_trace("decompress data: %s", current->file_page);	
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
				Leon *leonWrite = new Leon();
				const char* tmp_name = createFile(name, TMP_FILE);
				while(current!=NULL && current->next_page!=NULL){
					if(current->dirty != 0){
						if(current->isFastq){
							leonWrite->run(7, prep_args(name,true, 
								true, current->block_id));
						}else{
							leonWrite->run(6, prep_args(name,true,
                                                                false, current->block_id));
						}
						int fd = open(tmp_name,O_RDWR|O_TRUNC|O_CREAT);
						pwrite(fd, current->file_page, strlen(current->file_page), 0);
						leonWrite->executeCompression(current->block_id, tmp_name);
					}
					StackFS_trace("The flushed page: %s", current->file_page);
					current = current->next_page;
				}
				if(current!=NULL){
					if(current->dirty!=0){
						if(current->isFastq){
                                                        leonWrite->run(7, prep_args(name,true,
                                                                true, current->block_id));
                                                }else{
                                                        leonWrite->run(6, prep_args(name,true,
                                                                false, current->block_id));
                                                }
						int fd = open(tmp_name,O_RDWR|O_TRUNC|O_CREAT);
						pwrite(fd, current->file_page, strlen(current->file_page), 0);
						leonWrite->executeCompression(current->block_id, tmp_name);
					}
					StackFS_trace("The flushed page: %s", current->file_page);
					current = current->next_page;
				}
				delete leonWrite;
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
		compress_save(key);
		while(current!=NULL && current->next_page!=NULL){
			current = current->next_page;
		}
		while(current!=NULL && current!=fp){
			current = current->prev_page;
			if(current->next_page!=NULL){
				free(current->next_page);
				current->next_page = NULL;
			}
		}
		if(fp!=NULL){
			free(fp);
			fp = NULL;
		}
		pthread_spin_lock(&spinlock);
		g_hash_table_remove(file_table,key);
		g_hash_table_remove(file_cache, key);
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
	if(g_hash_table_size(file_cache) >= CACHE_SIZE){
		while(g_hash_table_size (file_cache)>= CACHE_SIZE){
			if(file_node_end->prev != file_node_head){
				char* rem_name = file_node_end->prev->file_name;
				safe_remove(rem_name);
				if(file_node_end->prev!=NULL){
					file_node_end->prev = file_node_end->prev->prev;
					free(file_node_end->prev->next);
					file_node_end->prev->next = file_node_end;
				}
			}else
				break;
		}
	}
	if(g_hash_table_contains(file_table,key) && g_hash_table_size(file_cache)< CACHE_SIZE){
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
}

static void stackfs_ll_open(fuse_req_t req, fuse_ino_t ino,
					struct fuse_file_info *fi)
{
	int fd;
	generate_start_time(req);
	char * name = lo_name(req, ino);
	char * point = NULL;
	StackFS_trace("The file name: %s", name);
	if((point = strrchr(name,'.')) != NULL && 
			(strcmp(point, ".fasta")==0||strcmp(point, ".fastq")==0
			||strcmp(point, ".fa")==0||strcmp(point, ".fq")==0)) {
		struct stat buffer;
		string given_name(name);
		bool isFastq = true;
		if(given_name.rfind(".fasta")!=string::npos){
			string base_name = given_name+"_0" + EXT;
			if(stat(base_name.c_str(), &buffer)!=0){
				size_t loc = given_name.rfind(".fasta");
				given_name.assign(given_name.substr(0, loc));
				given_name += ".fastq";
			}else
				isFastq = false;
		}else if(given_name.rfind(".fa")!=string::npos){
			string base_name = given_name + "_0" + EXT;
			if(stat(base_name.c_str(), &buffer)!=0){
				size_t loc = given_name.rfind(".fa");
				given_name.assign(given_name.substr(0, loc));
                                given_name += ".fq";
			}else
				isFastq = false;
		}	
		char* short_name = strdup(given_name.c_str());
		string file_name (given_name);
		file_name = file_name+"_0" + EXT; 
		/*int file_name_len = strlen(name)-strlen(point);
		char * file_name = (char *) malloc(file_name_len+EXT_LEN+1);
		strncpy(file_name, name, file_name_len);
		strncpy(file_name+file_name_len, EXT, EXT_LEN);
		file_name[file_name_len+EXT_LEN] = '\0';
		StackFS_trace("The file name: %s", file_name);*/
		if(stat(file_name.c_str(), &buffer) == 0) {
			fd = open(file_name.c_str(), fi->flags);
			decompress_store(short_name, fd);
		}
		else{
			fd = open(name, fi->flags);
		}
		free(short_name);
	}
	else{
		fd = open(name, fi->flags);
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
	
char* getOriginalFilename(char* name){
	char * point = NULL;
        if((point = strrchr(name,'.')) != NULL &&
                        (strcmp(point, ".fasta")==0||strcmp(point, ".fastq")==0
                        ||strcmp(point, ".fa")==0||strcmp(point, ".fq")==0)) {
                struct stat buffer;
                string given_name(name);
                bool isFastq = true;
                if(given_name.rfind(".fasta")!=string::npos){
                        string base_name = given_name+"_0" + EXT;
                        if(stat(base_name.c_str(), &buffer)!=0){
                                size_t loc = given_name.rfind(".fasta");
                                given_name.assign(given_name.substr(0, loc));
                                given_name += ".fastq";
                        }else
                                isFastq = false;
                }else if(given_name.rfind(".fa")!=string::npos){
                        string base_name = given_name + "_0" + EXT;
                        if(stat(base_name.c_str(), &buffer)!=0){
                                size_t loc = given_name.rfind(".fa");
                                given_name.assign(given_name.substr(0, loc));
                                given_name += ".fq";
                        }else
                                isFastq = false;
                }
                char* short_name = strdup(given_name.c_str());
                string file_name (given_name);
                file_name = file_name+"_0" + EXT;
		if(stat(file_name.c_str(), &buffer) == 0) {
			name = short_name;
                	return name;
                }
                else{
		       free(short_name);
                       return name;
                }
                //free(short_name);
        }
        else{
                return name;
        }
		
}

static void stackfs_ll_read(fuse_req_t req, fuse_ino_t ino, size_t size,
				off_t offset, struct fuse_file_info *fi)
{
	int res;
	(void) ino;
	//struct timespec start, end;
	//long time;
	//long time_sec;

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
	} else {
		char *buf;
		
		StackFS_trace("Read on name : %s, Kernel inode : %llu, fuse inode : %llu, off : %lu, size : %zu",
					lo_name(req, ino), get_lower_fuse_inode_no(req, ino), get_higher_fuse_inode_no(req, ino), offset, size);
		buf = (char *)malloc(size);
		generate_start_time(req);
		char * name = lo_name(req, ino);
		int fd = 0;
		//name = getOriginalFilename(name);
		if(g_hash_table_contains(file_table,name)){//||g_hash_table_contains(file_cache,name)) {
			/*if(g_hash_table_contains(file_cache,name)){
				fd = decompress(name);
			}*/
			struct file_pages* fp = (struct file_pages*)g_hash_table_lookup (file_table, name);
			struct file_pages* current = fp;
			struct file_pages* start = fp;
			struct file_pages* end = fp;
			struct file_pages* pre = fp;
			Leon* leon = fp->leon;
			bool block_found = false;
                	int page_no = offset/PG_SIZE;
			int page_off = offset%PG_SIZE;
			int end_off = page_off+size;
			vector<string>* out;
			int fromOff = 0, toOff = 0, bufOff = 0;
			int end_page = (offset+size)/PG_SIZE;
                        int fromBlock = leon->findBlockId(offset, fromOff);
                        int toBlock = leon->findBlockId(offset+size, toOff);
			vector<int>* block_sizes = leon->getBlockSizes();
			if(block_sizes!=NULL && fromBlock >= block_sizes->size()){
				fromBlock = block_sizes->size()-1;
				fromOff = (*block_sizes)[fromBlock]-1;
			}if(fromOff>(*block_sizes)[fromBlock])
				fromOff = (*block_sizes)[fromBlock]-1;
			if(block_sizes!=NULL && toBlock >= block_sizes->size()){
				toBlock = block_sizes->size()-1;
				toOff = (*block_sizes)[toBlock]-1;
			}if(toOff>(*block_sizes)[toBlock])
				toOff = (*block_sizes)[toBlock]-1;
			StackFS_trace("The from block: %d to block %d", fromBlock, toBlock);
			StackFS_trace("The fromOff: %d toOff %d", fromOff, toOff);
			while(current->block_id < fromBlock){
				//pthread_spin_lock(&spinlock);
				pre = current;
                                current = current->next_page;
				if(current==NULL || current->block_id > fromBlock){
					Leon *leonRead = new Leon();
                                        leonRead->run(6, prep_args(name,false,false,fromBlock));
                                        out = leonRead->executeDecompression(fromBlock);
                                        //const char* ret = leonRead->executeDecompression(fromBlock);
					delete leonRead;
                                        struct file_pages* tmp = pre->next_page;
                                        current = pre;
                                        current->next_page = (struct file_pages*) malloc(sizeof(struct file_pages));
                                        current->next_page->prev_page = current;
                                        current->next_page->next_page = tmp;
                                        current = current->next_page;
					//current->file_page = (char *) malloc (strlen(ret)+1);
					//strcpy(current->file_page,ret);
                                        current->file_page = (char *) malloc ((*out)[0].size()+1);
                                        strcpy(current->file_page,(*out)[0].c_str());
                                        int page_size = strlen(current->file_page)+1;
                                        current->s_off = page_size - fromOff;
                                        current->e_off = current->s_off + (*out)[0].size();
					//current->e_off = current->s_off + strlen(ret);
                                        current->block_id = fromBlock;
                                        current->leon = leon;
                                        current->isFastq = fp->isFastq;
                                        current->dirty = 0;
                                        current->sequenceAdded = false;
                                        delete out;
				}
				pthread_spin_unlock(&spinlock);
			}
			start = current;
			for(int i=fromBlock;i<=toBlock;i++){
				 pthread_spin_lock(&spinlock);
				if(current!=NULL && current->block_id == i){
					pre = current;
					current = current->next_page;
				}else{
					Leon *leonRead = new Leon();
                                        leonRead->run(6, prep_args(name,false,false,i));
                                	out = leonRead->executeDecompression(i);
					//const char* ret = leonRead->executeDecompression(i);
                                	delete leonRead;
					struct file_pages* tmp = pre->next_page;
					current = pre;
                                        current->next_page = (struct file_pages*) malloc(sizeof(struct file_pages));
                                        current->next_page->prev_page = current;
                                        current->next_page->next_page = tmp;
                                        current = current->next_page;
					//current->file_page = (char *) malloc (strlen(ret)+1);
                                        //strcpy(current->file_page,ret);
                                        current->file_page = (char *) malloc ((*out)[0].size()+1);
                                        strcpy(current->file_page,(*out)[0].c_str());
                                        int page_size = strlen(current->file_page)+1;
                                        current->s_off = page_size - fromOff;
                                        current->e_off = current->s_off + (*out)[0].size();
					//current->e_off = current->s_off + strlen(ret);
                                        current->block_id = i;
                                        current->leon = leon;
                                        current->isFastq = fp->isFastq;
                                        current->dirty = 0;
                                        current->sequenceAdded = false;
					pre = current;
					current = current->next_page;
					delete out;		
				}
				 pthread_spin_unlock(&spinlock);
			}
			/*bool change_made = false;
			current = start;
			//check for mismatch in block sizes, offsets from fromBlock to toBlock
			for(int i = fromBlock; i<=toBlock;i++){
				if( i==current->block_id && 
					(*block_sizes)[i]!= strlen(current->file_page)){
					(*block_sizes)[i] = strlen(current->file_page);
					change_made = true;
				}
				current = current->next_page;
			}
			if(change_made){
				 pthread_spin_lock(&spinlock);
				leon->saveConfig();
				fromBlock = leon->findBlockId(offset, fromOff);
				toBlock = leon->findBlockId(offset+size, toOff);
				if(block_sizes!=NULL && fromBlock >= block_sizes->size()){
					fromBlock = block_sizes->size()-1;
					fromOff = (*block_sizes)[fromBlock]-1;
				}if(fromOff>(*block_sizes)[fromBlock])
				fromOff = (*block_sizes)[fromBlock]-1;
				if(block_sizes!=NULL && toBlock >= block_sizes->size()){
					toBlock = block_sizes->size()-1;
					toOff = (*block_sizes)[toBlock]-1;
				}if(toOff>(*block_sizes)[toBlock])
				toOff = (*block_sizes)[toBlock]-1;
				current = start;
				for(int i=fromBlock;i<=toBlock;i++){
                                if(current!=NULL && current->block_id == i){
                                        pre = current;
                                        current = current->next_page;
                                }else{
                                        Leon *leonRead = new Leon();
                                        if(fp->isFastq)
                                                leonRead->run(7, prep_args(name,false,true,i));
                                        else
                                                leonRead->run(6, prep_args(name,false,false,i));
                                        out = leonRead->executeDecompression(i);
                                        delete leonRead;
                                        struct file_pages* tmp = pre->next_page;
                                        current = pre;
                                        current->next_page = (struct file_pages*) malloc(sizeof(struct file_pages));
                                        current->next_page->prev_page = current;
                                        current->next_page->next_page = tmp;
                                        current = current->next_page;
                                        current->file_page = (char *) malloc ((*out)[0].size()+1);
                                        strcpy(current->file_page,(*out)[0].c_str());
                                        int page_size = strlen(current->file_page)+1;
                                        current->s_off = page_size - fromOff;
                                        current->e_off = current->s_off + (*out)[0].size();
                                        current->block_id = i;
                                        current->leon = leon;
                                        current->isFastq = fp->isFastq;
                                        current->dirty = 0;
                                        current->sequenceAdded = false;
                                        pre = current;
                                        current = current->next_page;
                                        delete out;
                                }//else
				}//for
				 pthread_spin_unlock(&spinlock);
			}*/
			if(start!=NULL){
				if(fromBlock == toBlock){
					 memcpy(buf, start->file_page+fromOff, toOff-fromOff);
					 res = toOff-fromOff;
					 //StackFS_trace("Data: %s", buf);
				}else{
					memcpy(buf, start->file_page+fromOff, strlen(start->file_page)-fromOff);
					bufOff += strlen(start->file_page)-fromOff;
					start = start->next_page;
					for(int i=fromBlock+1;i<toBlock && start!=NULL; i++){
						memcpy(buf+bufOff, start->file_page ,strlen(start->file_page)); 
						bufOff+= strlen(start->file_page);
						start = start->next_page;
					}
					if(start!=NULL){
						memcpy(buf+bufOff, start->file_page, toOff);
						bufOff += toOff;
					}
					//StackFS_trace("Data: %s", buf);
					res = bufOff;
				}
			}
        	}else{	
		//clock_gettime(CLOCK_MONOTONIC, &start);
			StackFS_trace("The file name: %s not in map", name);
			res = pread(fi->fh, buf, size, offset);
		}
		//clock_gettime(CLOCK_MONOTONIC, &end);
		generate_end_time(req);
		populate_time(req);
		//time_sec = end.tv_sec - start.tv_sec;
		//time = end.tv_nsec - start.tv_nsec;
		//time_sec *= 1000000000;
		//time += time_sec;
		//StackFS_trace("Read inode : %llu off : %lu size : %zu diff : %llu", get_lower_fuse_inode_no(req, ino), offset, size, time);
		if (res == -1)
			return (void) fuse_reply_err(req, errno);
		res = fuse_reply_buf(req, buf, res);
		free(buf);
		buf = NULL;
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
	GHashTable* directory_list = g_hash_table_new(g_str_hash, g_str_equal);
	unordered_set<string> file_list;
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
		size_t dot = cur_file.find_last_of('.');
        	if(dot != string::npos && cur_file.substr(dot).compare(EXT)==0){
			if(cur_file.find(".fq") !=  string::npos || cur_file.find(".fastq") !=  string::npos){
				//cut the file name after the XXXX.fastq or XXXX.fq
				size_t under = cur_file.rfind(".fq");
				if(under==string::npos){
					under = cur_file.rfind(".fastq");
					under+=6;
				}else{
					under+=3;
				}
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

			}else if(cur_file.find(".fa") !=  string::npos || cur_file.find(".fasta") !=  string::npos)
			{
				// cut the given file name till 'XXXX.fast' or 'XXXX.f' and append 'a' afterwards
				/*size_t under = cur_file.rfind(".fq");
                                if(under==string::npos){
                                        under = cur_file.rfind(".fastq");
					if(under==string::npos){*/
						size_t under = cur_file.rfind(".fa");
						if(under==string::npos){
							under = cur_file.rfind(".fasta");
							under+=6;
						}else{
							under +=3;
						}
					/*}else{
                                        	under+=5;
					}
                                }else{
                                        under+=2;
				}*/
				string file_name = cur_file.substr(0, under);
                                //string file_name  = file_base_name + "a";
				//after the file will of the form "XXXX.fa" if it was "XXXX.fq" or 
				// "XXXX.fasta" if it was "XXXX.fastq"
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
		}else if((dot == string::npos || cur_file.substr(dot).compare(".qual")!=0)
					&& file_list.find(cur_file)==file_list.end()){
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
		/*
		char* point = NULL;
        	StackFS_trace("The file name: %s", d->entry->d_name);
        	if((point = strrchr(d->entry->d_name,'.')) != NULL && strcmp(point, EXT)==0) {
			int file_name_len = strlen(d->entry->d_name)-strlen(point);
                	char * file_name = (char *) malloc(file_name_len+TXT_EXT_LEN+1);
                	strncpy(file_name, d->entry->d_name, file_name_len);
                	strncpy(file_name+file_name_len, TXT_EXT, TXT_EXT_LEN);
                	file_name[file_name_len+TXT_EXT_LEN] = '\0';
                	StackFS_trace("The file name: %s", file_name);
			if(g_hash_table_lookup(directory_list, file_name)==NULL){
				struct stat st = {};
                        	st.st_ino = d->entry->d_ino;
                        	st.st_mode = d->entry->d_type << 12;
                		entsize = fuse_add_direntry(req, p, rem,
                                        	file_name, &st, nextoff);
				if (entsize > rem)
		                        break;

		                p += entsize;
                		rem -= entsize;
				g_hash_table_insert(directory_list, file_name, d->entry->d_name);
			}
		}
		if(g_hash_table_lookup(directory_list, d->entry->d_name)==NULL){
			struct stat st = {};
			st.st_ino = d->entry->d_ino;
			st.st_mode = d->entry->d_type << 12;
			entsize = fuse_add_direntry(req, p, rem,
					d->entry->d_name, &st, nextoff);*/
	/* The above function returns the size of the entry size even though
	* the copy failed due to smaller buf size, so I'm checking after this
	* function and breaking out incase we exceed the size.
	*/
		/*	if (entsize > rem)
			break;

			p += entsize;
			rem -= entsize;
			g_hash_table_insert(directory_list, d->entry->d_name, d->entry->d_name);
		}*/
		d->entry = NULL;
		d->offset = nextoff;
	}
	StackFS_trace("The list are %s", p);
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

static void stackfs_ll_release(fuse_req_t req, fuse_ino_t ino,
					struct fuse_file_info *fi)
{
	(void) ino;

	StackFS_trace("Release called on name : %s and inode : %llu fd : %d ",
			lo_name(req, ino), lo_inode(req, ino)->ino, fi->fh);
	generate_start_time(req);
	char * name = lo_name(req, ino);
	if(g_hash_table_lookup(file_table, name)!=NULL && g_hash_table_lookup(file_cache, name)==NULL)
		add_to_cache(name, fi->fh);
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
	generate_start_time(req);
	closedir(d->dp);
	generate_end_time(req);
	populate_time(req);
	free(d);
	fuse_reply_err(req, 0);
}

static void stackfs_ll_write(fuse_req_t req, fuse_ino_t ino, const char *buf,
			size_t size, off_t off, struct fuse_file_info *fi)
{
	int res;
	(void) ino;
	StackFS_trace("Write name : %s, inode : %llu, off : %lu, size : %zu",
			lo_name(req, ino), lo_inode(req, ino)->ino, off, size);
	generate_start_time(req);
	char * name = lo_name(req, ino);
	int compressedType = 1;
	if(g_hash_table_contains(file_table,name)) {
		struct file_pages* fp = (struct file_pages*)g_hash_table_lookup (file_table, name);
		struct file_pages* current = fp;
		struct file_pages* start = fp;
		if(compressedType == 0){
		int page_no = off/PG_SIZE;
		int page_off = off%PG_SIZE;
		int end_off = page_off+size;
		int end_page = (off+size)/PG_SIZE;
		StackFS_trace("The file name: %d %d", page_no, end_off);
		while(current!=NULL && current->next_page!=NULL && page_no>0){
			current = current->next_page;
			start = start->next_page;
			page_no--;
			end_page--;
		}
		StackFS_trace("The page_no, end_page: %d %d", page_no, end_page);
		if(current!=NULL && current->next_page == NULL && end_page >= 1){
			char * point = NULL;
			int fd = 0;
			int bytes_read = PG_SIZE;
			StackFS_trace("The file name: %s", name);
			if((point = strrchr(name,'.')) != NULL) {
				struct stat buffer;
				int file_name_len = strlen(name)-strlen(point);
				char * file_name = (char *) malloc(file_name_len+EXT_LEN+1);
				strncpy(file_name, name, file_name_len);
				strncpy(file_name+file_name_len, EXT, EXT_LEN);
				file_name[file_name_len+EXT_LEN] = '\0';
				StackFS_trace("Retrieving more data from file: %s", file_name);
				if(stat(file_name, &buffer) == 0) {
					fd = open(file_name, fi->flags);
				}
			}
			while(bytes_read == PG_SIZE && end_page!=0){
				current->next_page = (struct file_pages*) malloc(sizeof(struct file_pages));
				current->next_page->prev_page = current;
				current->next_page->next_page = NULL;
				//current->next_page->offset= current->offset;
				current = current->next_page;
				//bytes_read = pread(fd, current->file_page, PG_SIZE, current->offset);
				int off = 0;
				while(off<bytes_read){
					current->file_page[off] +=1;
					off++;
				}
				//current->offset+=bytes_read;
				end_page--;
			}
			while(end_page!=0){
				current->next_page = (struct file_pages*) malloc(sizeof(struct file_pages));
                                current->next_page->prev_page = current;
                                current->next_page->next_page = NULL;
                                //current->next_page->offset= current->offset;
                                current = current->next_page;
                                end_page--;
			}	
			while(page_no>0){				
				start = start->next_page;
				page_no--;
			}
			current->next_page = NULL;
			StackFS_trace("decompress data: %s %d", current->file_page, bytes_read);
		}
		StackFS_trace("The page_no, end_page: %d %d", page_no, end_page);
		if(start!=NULL){
			if(end_off<=PG_SIZE){
				StackFS_trace("The page before insert: %s", start->file_page);
				memcpy(start->file_page+page_off, buf, size);
				StackFS_trace("The page after insert: %s", start->file_page);
				start->dirty= 1;
				res = size;
			}else{
				StackFS_trace("The page before insert: %s", start->file_page);
				memcpy(start->file_page+page_off, buf, PG_SIZE-page_off);
				end_off -= (PG_SIZE-page_off);
				start = start->next_page;
				int buf_off = PG_SIZE-page_off;
				while(start!=NULL && end_off>=PG_SIZE){
					memcpy(start->file_page, buf+buf_off ,PG_SIZE);
					end_off -= PG_SIZE;
					buf_off+= PG_SIZE;
					start->dirty = 1;
					start = start->next_page;
				}
				if(start!=NULL){
					memcpy(start->file_page, buf+buf_off, end_off);
					buf_off += end_off;
					start->dirty = 1;
					end_off = 0;
				}
				res = buf_off;
			}
			StackFS_trace("Written data: %s", start->file_page);
		}
		}
		else if(compressedType == 1){
                	struct file_pages* end = fp;
			vector<string>* out;
			Leon* leon = fp->leon;
                        bool block_found = false;
                        int fromOff = 0, toOff = 0, bufOff = 0;
                        int fromBlock = leon->findBlockId(off, fromOff);
                        int toBlock = leon->findBlockId(off+size, toOff);
                        StackFS_trace("The from block: %d to block %d", fromBlock, toBlock);
			while(current!=NULL){
                                if(current->block_id == fromBlock){
                                        start = current;
                                        for(int i=fromBlock;i<toBlock && current!=NULL;i++){
                                                if(current->block_id == i)
                                                        current = current->next_page;
                                                else{
                                                        end = current;
                                                        break;
                                                }
                                        }
                                        if(current!=NULL && current->block_id == toBlock){
                                                block_found = true;
                                                end = current->next_page;
                                        }else if(current == NULL)
                                                end = current;
                                        break;
                                }
                                if(current->block_id>fromBlock){
                                        end = current;
                                        start = current->prev_page;
                                        break;
                                }
                                start = current;
                                current = current->next_page;
                        }
			if(start!=NULL && !block_found){
                                char* file_name = change_file_ext(name, EXT);
                                Leon *leonWrite = new Leon();
                                leonWrite->run(6, prep_args(file_name,false, false, 0));
                                out = leonWrite->executeDecompression(fromBlock);
				//const char* ret = leonWrite->executeDecompression(fromBlock);
                                delete leonWrite;
                                current = start;
                                //while(bytes_read == PG_SIZE && end_page!=0){
                                for(int i=fromBlock;i<= toBlock;i++){
                                        StackFS_trace(" populating block id: %d", current->block_id);
                                        if(current->block_id == i - 1 && end!=NULL
                                                        && end->block_id == i){
                                                current->next_page = end;
                                                current = current->next_page;
                                                end = end->next_page;
                                        }else{
                                        int index = i-fromBlock;
                                        current->next_page = (struct file_pages*) malloc(sizeof(struct file_pages));
                                        current->next_page->prev_page = current;
                                        current->next_page->next_page = NULL;
                                        current = current->next_page;
					//current->file_page = (char *) malloc (strlen(ret)+1);
                                        //strcpy(current->file_page,ret);
                                        current->file_page = (char *) malloc ((*out)[index].size()+1);
                                        strcpy(current->file_page,(*out)[index].c_str());
                                        int page_size = strlen(current->file_page)+1;
                                        current->s_off = page_size - fromOff;
                                        current->e_off = current->s_off + (*out)[index].size();
                                        current->block_id = i;
                                        current->leon = leon;
                                        current->dirty = 0;
                                        current->sequenceAdded = false;
                                        }
                                }
                                delete out;
			}
                        if(start!=NULL){
                                if(fromBlock == toBlock){
                                         memcpy(start->file_page+fromOff,buf,size);
                                         res = size;
					 start->dirty = 1;
                                         //StackFS_trace("Data: %s", buf);
                                }else{
                                        memcpy(start->file_page+fromOff, buf, strlen(start->file_page)-fromOff);
                                        start = start->next_page;
                                        bufOff += strlen(start->file_page)-fromOff;
                                        for(int i=fromBlock+1;i<toBlock && start!=NULL; i++){
                                                memcpy(start->file_page, buf+bufOff ,strlen(start->file_page));
                                                bufOff+= strlen(start->file_page);
						start->dirty = 1;
                                                start = start->next_page;
                                        }
                                        if(start!=NULL){
                                                memcpy(start->file_page,buf+bufOff, toOff);
						start->dirty =1;
                                                bufOff += toOff;
                                        }
                                        StackFS_trace("Data: %s", buf);
                                        res = bufOff;
                                }
                        }
		}
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
	res = unlink(fullPath);
	generate_end_time(req);
	populate_time(req);
	if (res == -1)
		fuse_reply_err(req, errno);
	else
		fuse_reply_err(req, res);

	if (fullPath)
		free(fullPath);
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

#if  TESTING_XATTR
static void stackfs_ll_getxattr(fuse_req_t req, fuse_ino_t ino,
					const char *name, size_t size)
{
	int res;

	StackFS_trace("Function Trace : Getxattr");
	if (size) {
		char *value = (char *) malloc(size);

		generate_start_time(req);
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
		res = lgetxattr(lo_name(req, ino), name, NULL, 0);
		generate_end_time(req);
		populate_time(req);
		if (res >= 0)
			fuse_reply_xattr(req, res);
		else
			fuse_reply_err(req, errno);
	}
}
#endif

static struct fuse_lowlevel_ops hello_ll_oper = {};
/*	lookup	:stackfs_ll_lookup,
	forget  :stackfs_ll_forget,
	getattr	:stackfs_ll_getattr,
	setattr :stackfs_ll_setattr,
	mkdir   :stackfs_ll_mkdir,
	unlink  :stackfs_ll_unlink,
	rmdir   :stackfs_ll_rmdir,
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
#if     TESTING_XATTR
        hello_ll_oper.getxattr = stackfs_ll_getxattr;
#endif
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
					if (multithreaded)
						err = fuse_session_loop_mt(se);
					else
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