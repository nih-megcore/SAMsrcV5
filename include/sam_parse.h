/* Parse SAM parameter files. */

#include "sam_param.h"

/* Parameter linked list. */

typedef struct parm {
    struct parm_table *param; /* pointer into parameter table */
    void *ptr;              /* where to stash the result of the parse */
    int set;                /* true if this parameter was specified */
    int mnum;               /* marker number or 0 */
    struct parm *next;      /* linked list pointer */
} PARM;

/* Parameters are described by entries in a table. */

typedef struct parm_table {
    int cmdline;            /* TRUE if this can be on the command line */
    char *key;              /* parameter name, aka --key */
    char *flag;             /* short (-o) version, or NULL */
    int (*parsefn)(int, char **, int, struct parm *); /* parse function */
    char *arg_help;         /* description of arguments, or NULL */
    char *help;             /* help string */
    char *env_var;          /* environment variable name or NULL */
} PARM_TABLE;

/* The %include mechanism uses this linked list. */

typedef struct inc_list {
    char *filename;         /* name of included file */
    struct inc_list *next;
} INC_LIST;

extern PARM_TABLE Parm_table[];

extern void reg_usage(char *progname, int minor, char *date);
extern void do_help();
extern void reg_parm(char *key);
extern void set_parm_arg_help(char *key, char *help);
extern void set_parm_help(char *key, char *help);
extern void reg_std_parm();
extern void do_parse_args(int argc, char **argv);
extern PARM *get_parm(char *key);
extern void parse_args(int argc, char **argv);
extern void parse_env();
extern void parse_file(char *name, int required);

extern int helpfn(int argc, char **argv, int, PARM *);
extern int boolfn(int argc, char **argv, int, PARM *);
extern int intfn(int argc, char **argv, int, PARM *);
extern int doublefn(int argc, char **argv, int, PARM *);
extern int double2fn(int argc, char **argv, int, PARM *);
extern int double3fn(int argc, char **argv, int, PARM *);
extern int trangefn(int argc, char **argv, int, PARM *);
extern int prefixfn(int argc, char **argv, int, PARM *);
extern int stringfn(int argc, char **argv, int, PARM *);
extern int includefn(int argc, char **argv, int, PARM *);
extern int markerfn(int argc, char **argv, int, PARM *);
extern int segfilefn(int argc, char **argv, int, PARM *);
extern int imgformatfn(int argc, char **argv, int, PARM *);
extern int mufn(int argc, char **argv, int, PARM *);
extern int modelfn(int argc, char **argv, int, PARM *);
extern int imgmetricfn(int argc, char **argv, int, PARM *);

extern int get_params(PARMINFO *);      // call new_params(), then do_get_params()
extern void new_params(PARMINFO *);     // fill in defaults
extern int do_get_params(PARMINFO *);   // fill in from Parm_list

extern void log_params(char *);         // write active parameters to a file
