/* SAM parameter parsing. Reads default parameters files,
with overrides from command line and environment variables. */

#include "samutil.h"
#include "sam_parse.h"
#include "sam_param.h"
#include "version.h"
#include <time.h>
#include <unistd.h>
#include <ctype.h>

/* Parameter linked list head. */

static PARM *Parm_list = NULL;
static PARM *Parm_list_tail;            /* pointer to the last element */
static int Parm_list_nmarkers = 0;      /* how many markers we've seen */

/* List of %include files. */

static INC_LIST *Inc_list = NULL;
static INC_LIST *Inc_list_tail;         /* so we can add to the end */

/* Help message. */

static char *Usage;

void reg_usage(char *progname, int minor, char *date)
{
    int i;
    char *s;

    Progname = progname;
    if ((s = rindex(progname, '/')) != NULL) {
        Progname = s + 1;
    }

    i = strlen(PRG_REV) + strlen(date) + 20;
    Usage = new_string(i);
    snprintf(Usage, i, "%s rev-%0d, %s", PRG_REV, minor, date);
}

/* Call the help function as if -h were specified. */

void do_help()
{
    helpfn(0, NULL, 0, NULL);
}

/* Create a new PARM and add it to the end of the Parm_list. */

static PARM *new_parm(PARM_TABLE *param)
{
    PARM *p;

    p = new(PARM);
    p->param = param;
    p->ptr = NULL;
    p->set = FALSE;
    p->mnum = 0;
    p->next = NULL;

    /* Add it to the end. */

    if (Parm_list == NULL) {
        Parm_list = p;
    } else {
        Parm_list_tail->next = p;
    }
    Parm_list_tail = p;

    return p;
}

/* Register a new parameter. */

void reg_parm(char *key)
{
    PARM_TABLE *param;

    /* Search the table for the key, add a new PARM to Parm_list when found. */

    for (param = Parm_table; param->key; param++) {
        if (strcmp(key, param->key) == 0) {
            (void)new_parm(param);
            return;
        }
    }
    fatalerr("unknown parameter: %s", key);
}

/* Helper to set param->arg_help. */

void set_parm_arg_help(char *key, char *help)
{
    PARM_TABLE *param;

    /* Search the table for the key. */

    for (param = Parm_table; param->key; param++) {
        if (strcmp(key, param->key) == 0) {
            param->arg_help = help;
            return;
        }
    }
    fatalerr("unknown parameter: %s", key);
}

/* Helper to set param->help. */

void set_parm_help(char *key, char *help)
{
    PARM_TABLE *param;

    /* Search the table for the key. */

    for (param = Parm_table; param->key; param++) {
        if (strcmp(key, param->key) == 0) {
            param->help = help;
            return;
        }
    }
    fatalerr("unknown parameter: %s", key);
}

/* Register standard parameters used by all programs. */

void reg_std_parm()
{
    reg_parm("help");
    reg_parm("DataSet");
#if BTI
    reg_parm("PDFName");
#endif
    reg_parm("param");
    reg_parm("verbose");
    reg_parm("%include");
}

/* Do standard parsing of arguments. */

void do_parse_args(int argc, char **argv)
{
    PARM *p;
    INC_LIST *inc;
    char *s, *home, *path;

    /* Parse the command line arguments and environment variables. */

    parse_args(argc, argv);
    parse_env();

    /* If a paramter file was specified, parse it. */

    p = get_parm("param");
    if (p->set) {
        parse_file((char *)p->ptr, TRUE);
    }

    /* Check for any included files. */

    p = get_parm("%include");
    if (p->ptr) {
        inc = (INC_LIST *)p->ptr;
        while (inc) {
            parse_file(inc->filename, TRUE);
            inc = inc->next;
        }
    }

    /* Check for optional default parameter files. */

    parse_file(".samrc", FALSE);
    parse_file(".coregrc", FALSE);

    home = getenv("HOME");
    if (home) {
        path = new_array(char, strlen(home) + 10);
        s = strecpy(path, home);
        strcpy(s, "/.samrc");
        parse_file(path, FALSE);
        s = strecpy(path, home);
        strcpy(s, "/.coregrc");
        parse_file(path, FALSE);
        free(path);
    }
}

/* This just checks to see if the parameter name is known. Allow
lowercase and abbreviations. */

static int valid_parm(char *key)
{
    int n, i;
    PARM_TABLE *p;

    /* Search the table for the key. */

    i = 0;
    n = strlen(key);
    for (p = Parm_table; p->key; p++) {
        if (strncasecmp(p->key, key, n) == 0) {
            i++;
        }
    }

    if (i == 1) {
        return TRUE;
    }

    msg("unknown parameter '%s'\n", key);
    do_help();  /* doesn't return */

    return FALSE;
}

/* Given a key, find the PARM structure. Allow lowercase and
abbreviations. If do_check is set, and the parameter isn't registered,
we just ignore it. Otherwise it's a fatal error. Markers can be
repeated, so get_param() isn't used to retrieve them. */

static PARM *_get_parm(char *key, int do_check)
{
    int n, i;
    int mnum = 0;
    PARM *p;
    PARM *q = NULL;

    /* Markers are special. Convert "MarkerN" into "Marker". Ignore N; the
    count of markers is handled separately. SegFile also specifies a marker. */

    if (strncasecmp(key, "Marker", 6) == 0) {
        key = "Marker";
        mnum = ++Parm_list_nmarkers;        // mnum is IO 1
    } else if (strncasecmp(key, "SegFile", 7) == 0) {
        mnum = ++Parm_list_nmarkers;
    }

    /* Find the PARM structure (created when the parameter was registered). */

    i = 0;
    n = strlen(key);
    for (p = Parm_list; p; p = p->next) {
        if (strncasecmp(p->param->key, key, n) == 0) {
            q = p;
            i++;
        }
    }

    if (i && mnum) {

        /* We have a new (registered) marker, so create a new PARM for it. */

        q = new_parm(q->param);
        q->mnum = mnum;
        i = 1;
    }

    if (i == 0) {
        if (mnum) {
            Parm_list_nmarkers--;   // if it's an unregistered marker, don't count it,
            return NULL;            // ignore it
        }
        if (valid_parm(key)) {
            if (!do_check) {
                fatalerr("parameter '%s' not registered", key);
            }
            return NULL;
        }
    }

    if (i != 1) {
        msg("ambiguous parameter '%s'\n", key);
        do_help();  /* doesn't return */
    }

    return q;
}

/* This is used while scanning the parameters to check
whether this parameter is already set. Valid but unregistered
parameters are OK during parsing, they are just ignored. */

static PARM *check_get_parm(char *key)
{
    return _get_parm(key, 1);
}

/* Return a parameter value. It's a fatal error if the parameter
is unregistered. */

PARM *get_parm(char *key)
{
    return _get_parm(key, 0);
}

/* Parse the command line options. */

void parse_args(int argc, char **argv)
{
    int i, k, j;
    char *s;
    PARM *parm;
    PARM_TABLE *p;

    /* The parse functions return the number of arguments consumed, or zero
    on error. */

    for (i = 0; i < argc; i++) {
        s = argv[i];
        if (s[0] == '-') {
            if (s[1] == '-') {

                /* --key */

                parm = get_parm(&s[2]);
                if (parm == NULL || !parm->param->cmdline) {

                    /* Don't allow unregistered or certain other parameters
                    on the command line. */

                    msg("parameter \"%s\" is not valid on the command line\n", s + 2);
                    do_help();  /* doesn't return */
                }
                k = (*parm->param->parsefn)(argc, argv, i + 1, parm);
                if (k) {
                    i += k;
                }
            } else if (isalpha(s[1])) {

                /* search the list for s (exact match) */

                j = TRUE;
                for (parm = Parm_list; parm; parm = parm->next) {
                    p = parm->param;
                    if (p->flag && strcmp(p->flag, s + 1) == 0) {
                        k = (*p->parsefn)(argc, argv, i + 1, parm);
                        if (k) {
                            i += k;
                        }
                        j = FALSE;
                    }
                }
                if (j) {
                    fatalerr("unknown option \"%s\"", s);
                    //msg("unknown option \"%s\"\n", s);
                    //do_help();  /* doesn't return */
                }
            }
        }
    }
}

/* Parse the arguments in the string s. */

static void parse_string(char *key, char *s)
{
    int i, k, argc;
    char *t, **argv;
    PARM *p;

    /* If the key is already set, don't do anything. This allows command
    line and environment variables to override. */

    p = check_get_parm(key);
    if (p == NULL || p->set) {
        return;
    }

    /* Convert space delimited words into argv[] style pointers. First,
    count the words. */

    i = 1;                              // 1 for the parameter name
    k = 0;                              // looking for a word
    for (t = s; *t; t++) {
        if (k == 0 && !isspace(*t)) {   // found a word
            i++;
            k = 1;                      // looking for space
        }
        if (k == 1 && isspace(*t)) {    // end of word
            k = 0;
        }
    }
    argc = i;

    /* Allocate an array and point to the words. */

    argv = new_array(char *, argc + 1);
    argv[0] = key;
    t = s;
    for (i = 1; i < argc; i++) {
        argv[i] = t;
        while (*t && !isspace(*t)) {
            t++;
        }
        if (isspace(*t)) {
            *t++ = '\0';
            while (*t && isspace(*t)) {
                t++;
            }
        }
    }
    argv[argc] = NULL;

    /* Call the parse callback for this key. */

    k = (*p->param->parsefn)(argc, argv, 1, p);
    if (k) {
        if (k + 1 != argc) {
            msg("[extra arguments for %s]\n", key);
        }
    }

    free(argv);
}

/* Parse environment variables, key=value. */

void parse_env()
{
    int i, n;
    PARM *parm;
    PARM_TABLE *p;

    /* Go through the list of parameters. For each one that has an env_var,
    see if it's present in the environment. Must be an exact match, and also
    allowed on the command line. */

    for (parm = Parm_list; parm; parm = parm->next) {
        p = parm->param;
        if (p->env_var && p->cmdline) {
            n = strlen(p->env_var);
            for (i = 0; environ[i]; i++) {
                if (strncmp(p->env_var, environ[i], n) == 0) {
                    if (environ[i][n] == '=') {

                        /* @@@ There's no trimming of the environment string here. */

                        parse_string(p->key, environ[i] + n + 1);
                    }
                }
            }
        }
    }
}

/* Parse the lines of a file. */

void parse_file(char *name, int required)
{
    char *s, *t, *k, *n;
    FILE *f;
    char buf[256];

    /* Ignore missing optional files. Complain about required ones. */

    n = name;
    if (!fileexists(name)) {

        /* Try again with .param added. */

        s = new_string(strlen(name) + 7);
        strcpy(strecpy(s, name), ".param");
        if (!fileexists(s)) {
            if (required) {
                fatalerr("file not found: %s", name);
            }
            free(s);
            return;
        }
        n = s;
    }
    f = fileopen(n, "r");
    if (n != name) {
        free(n);
    }

    /* Go through the lines of the file. */

    while (fgetline(buf, sizeof(buf), f)) {

        /* Remove leading spaces. Ignore blank lines and comments. */

        s = buf;
        while (*s && isspace(*s)) {
            s++;
        }
        if (*s == '\0' || *s == '#') {
            continue;
        }

        /* Get a space delimited keyword. */

        k = s;
        while (*s && !isspace(*s)) {
            s++;
        }

        /* If there are more arguments, parse them out. */

        t = s;
        if (isspace(*s)) {
            *s++ = '\0';    /* end the keyword */

            /* Skip to the start of the arguments. */

            while (*s && isspace(*s)) {
                s++;
            }
            t = s;

            /* Find the end, delete any comment. */

            while (*s && *s != '#') {
                s++;
            }
            *s = '\0';

            /* Trim trailing spaces. */

            while (isspace(s[-1])) {
                s--;
            }
            *s = '\0';
        }

        parse_string(k, t);
    }

    fclose(f);
}

/* Parameter callbacks. First, -h. Note: these callbacks are used for both
command line and parameter file processing. On the command line, there's
just one (argc, argv) pair, and i indicates the start of the arguments for
that option. In a file, (argc, argv) are the arguments to a single
parameter, and i is 1. The callbacks return the number of arguments
consumed. */

static void msg2(char *s0)
{
    char *s, *t;

    s = copy_string(s0);
    while (1) {
        t = index(s, '\n');
        if (t == NULL) {
            msg("%s", s);
            return;
        }
        *t++ = '\0';
        msg("%s\n", s);
        msg("%-29s", "");
        s = t;
    }
}

int helpfn(int argc, char **argv, int i, PARM *unused)
{
    int m;
    char *dash;
    PARM *parm;
    PARM_TABLE *p;
    char buf[200];

    msg("\nUsage: %s [options]\t%s\n\n", Progname, Usage);

    msg("Options:\n");
    msg("If a parameter begins with '--', it is allowed on the command line,\n");
    msg("otherwise it is only allowed in a parameter file (see -m).\n");
    msg("All times are in seconds.\n\n");

    m = FALSE;  /* haven't seen a Marker */

    for (parm = Parm_list; parm; parm = parm->next) {
        p = parm->param;
        if (p->key[0] == '%') { /* magic, don't show */
            continue;
        }
        if (p->flag) {
            if (p->arg_help) {
                snprintf(buf, sizeof(buf), "-%s %s, --%s %s", p->flag, p->arg_help, p->key, p->arg_help);
            } else {
                snprintf(buf, sizeof(buf), "-%s, --%s", p->flag, p->key);
            }
        } else {
            dash = "";          /* command line not allowed */
            if (p->cmdline) {
                dash = "--";    /* command line allowed */
            }
            if (p->arg_help) {

                /* Marker is special. */

                if (strcmp(p->key, "Marker") == 0) {
                    if (!m) {
                        m = TRUE;
                        snprintf(buf, sizeof(buf), "Marker %s", p->arg_help);
                    } else {
                        continue;
                    }
                } else {
                    snprintf(buf, sizeof(buf), "%s%s %s", dash, p->key, p->arg_help);
                }
            } else {
                snprintf(buf, sizeof(buf), "%s%s", dash, p->key);
            }
        }
        msg("    %-25s", buf);
        if (strlen(buf) >= 25) {
            msg("\n%-29s", "");
        }
        msg2(p->help);
        if (p->env_var) {
            msg(" [env_var = %s]", p->env_var);
        }
        msg("\n");
    }

    exit(0);
}

/* A simple bool that takes no arguments. */

int boolfn(int argc, char **argv, int i, PARM *p)
{
    p->set = TRUE;
    return 0;
}

/* Macros to assist with counting the required numbers of arguments. */

#define REQ1() { \
    int j; \
    char *s; \
    j = FALSE; \
    if (i < argc) { \
        s = argv[i]; \
        if (s[0] == '-' && (isalpha(s[1]) || s[1] == '-')) { \
            j = TRUE; \
        } \
    } \
    if (j || i == argc) { \
        fatalerr("%s requires an argument", argv[i-1]); \
    } \
}

#define REQ2() { \
    int j; \
    char *s; \
    j = FALSE; \
    if (i < argc) { \
        s = argv[i]; \
        if (s[0] == '-' && (isalpha(s[1]) || s[1] == '-')) { \
            j = TRUE; \
        } \
    } \
    if (i+1 < argc) { \
        s = argv[i+1]; \
        if (s[0] == '-' && (isalpha(s[1]) || s[1] == '-')) { \
            j = TRUE; \
        } \
    } \
    if (j || i+1 >= argc) { \
        fatalerr("%s requires two arguments", argv[i-1]); \
    } \
}

#define REQ3() { \
    int j; \
    char *s; \
    j = FALSE; \
    if (i < argc) { \
        s = argv[i]; \
        if (s[0] == '-' && (isalpha(s[1]) || s[1] == '-')) { \
            j = TRUE; \
        } \
    } \
    if (i+2 < argc) { \
        s = argv[i+2]; \
        if (s[0] == '-' && (isalpha(s[1]) || s[1] == '-')) { \
            j = TRUE; \
        } \
    } \
    if (j || i+2 >= argc) { \
        fatalerr("%s requires three arguments", argv[i-1]); \
    } \
}

/* An int. */

int intfn(int argc, char **argv, int i, PARM *p)
{
    long l;
    char *s;

    REQ1();

    l = strtol(argv[i], &s, 10);
    if (*s != '\0') {
        fatalerr("%s: badly formed number '%s'", argv[i-1], argv[i]);
    }
    p->ptr = new(int);
    *(int *)p->ptr = (int)l;
    p->set = TRUE;

    return 1;
}

/* A double. */

int doublefn(int argc, char **argv, int i, PARM *p)
{
    double d;
    char *s;

    REQ1();

    d = strtod(argv[i], &s);
    if (*s != '\0') {
        fatalerr("%s: badly formed number '%s'", argv[i-1], argv[i]);
    }
    p->ptr = new(double);
    *(double *)p->ptr = d;
    p->set = TRUE;

    return 1;
}

/* A pair of doubles. */

int double2fn(int argc, char **argv, int i, PARM *p)
{
    int j;
    double *d;
    char *s;

    REQ2();

    d = new_array(double, 2);

    j = 0;
    d[0] = strtod(argv[i], &s);
    if (*s != '\0') {                   // @@@ this is probably wrong
        j = i;
    }
    d[1] = strtod(argv[i+1], &s);
    if (*s != '\0') {
        j = i+1;
    }
    if (j) {
        fatalerr("%s: badly formed number '%s'", argv[i-1], argv[j]);
    }

    p->ptr = d;
    p->set = TRUE;

    return 2;
}

/* 3 doubles */

int double3fn(int argc, char **argv, int i, PARM *p)
{
    int j;
    double *d;
    char *s;

    REQ3();

    d = new_array(double, 3);

    j = 0;
    d[0] = strtod(argv[i], &s);
    if (*s != '\0') {
        j = i;
    }
    d[1] = strtod(argv[i+1], &s);
    if (*s != '\0') {
        j = i+1;
    }
    d[2] = strtod(argv[i+2], &s);
    if (*s != '\0') {
        j = i+2;
    }
    if (j) {
        fatalerr("%s: badly formed number '%s'", argv[i-1], argv[j]);
    }

    p->ptr = d;
    p->set = TRUE;

    return 3;
}

/* A time range. */

int trangefn(int argc, char **argv, int i, PARM *p)
{
    int j;
    TIMERANGE *r;
    char *s;

    REQ2();

    r = new(TIMERANGE);

    j = 0;
    r->t0 = strtod(argv[i], &s);
    if (*s != '\0') {
        j = i;
    }
    r->t1 = strtod(argv[i+1], &s);
    if (*s != '\0') {
        j = i+1;
    }
    if (j) {
        fatalerr("%s: badly formed number '%s'", argv[i-1], argv[j]);
    }

    p->ptr = r;
    p->set = TRUE;

    return 2;
}

/* The prefix may be specified as a number or as a character delimiter. */

int prefixfn(int argc, char **argv, int i, PARM *p)
{
    long l;
    char *s;

    REQ1();

    // if it's a number, use that

    l = strtol(argv[i], &s, 10);
    if (*s == '\0') {
        p->ptr = new(int);
        *(int *)p->ptr = (int)l;
        p->mnum = TRUE;             // flag that it's a number
    } else {

        // not a number; if it's a single character, use that

        if (strlen(argv[i]) == 1) {
            p->ptr = copy_string(argv[i]);
        } else {
            fatalerr("%s: badly formed Prefix '%s'", argv[i-1], argv[i]);
        }
    }
    p->set = TRUE;

    return 1;
}

/* A string argument cannot have spaces in it. */

int stringfn(int argc, char **argv, int i, PARM *p)
{
    REQ1();

    p->ptr = copy_string(argv[i]);
    p->set = TRUE;
    return 1;
}

/* %include makes a list of files to parse later. */

int includefn(int argc, char **argv, int i, PARM *p)
{
    INC_LIST *inc;

    REQ1();

    inc = new(INC_LIST);
    inc->filename = copy_string(argv[i]);
    inc->next = NULL;

    /* Add it to the end. */

    if (Inc_list == NULL) {
        Inc_list = inc;
    } else {
        Inc_list_tail->next = inc;
    }
    Inc_list_tail = inc;

    p->ptr = Inc_list;  // so we know there are some things in the list

    // don't set p->set so we can come here again
    return 1;
}

/* Parse MarkerN lines. */

int markerfn(int argc, char **argv, int i, PARM *p)
{
    int n;
    char *s;
    MARKINFO *m;

    if (argc < 5) {
        fatalerr("missing arguments for \"%s\"", argv[0]);
    }

    m = new(MARKINFO);
    m->MarkName = copy_string(argv[1]);
    m->MarkName2 = NULL;
    m->MarkStart = strtod(argv[2], &s);
    if (*s != '\0') {
        fatalerr("badly formed T0 for \"%s\"", argv[0]);
    }
    m->MarkEnd = strtod(argv[3], &s);
    if (*s != '\0') {
        fatalerr("badly formed T1 for \"%s\"", argv[0]);
    }
    n = strlen(argv[4]);
    if (strncasecmp(argv[4], "TRUE", n) == 0) {
        m->Sum = TRUE;
    } else if (strncasecmp(argv[4], "FALSE", n) == 0) {
        m->Sum = FALSE;
    } else {
        fatalerr("SUMFLAG must be TRUE or FALSE for \"%s\"", argv[0]);
    }
    if (argc == 6) {
        m->MarkName2 = copy_string(argv[5]);
    }
    m->SegFile = FALSE;
    m->FileName = NULL;     /* not a segfile */

    p->ptr = m;
    p->set = TRUE;

    return argc - 1;
}

/* SegFile name file sumflag. */

int segfilefn(int argc, char **argv, int i, PARM *p)
{
    int n;
    MARKINFO *m;

    if (argc != 4) {
        fatalerr("usage: %s markname filename sumflag", argv[0]);
    }

    m = new(MARKINFO);
    m->MarkName = copy_string(argv[1]);
    m->MarkName2 = NULL;
    m->MarkStart = 0;
    m->MarkEnd = 0;
    n = strlen(argv[3]);
    if (strncasecmp(argv[3], "TRUE", n) == 0) {
        m->Sum = TRUE;
    } else if (strncasecmp(argv[3], "FALSE", n) == 0) {
        m->Sum = FALSE;
    } else {
        fatalerr("SUMFLAG must be TRUE or FALSE for \"%s\"", argv[0]);
    }
    m->SegFile = TRUE;
    m->FileName = copy_string(argv[2]);

    p->ptr = m;
    p->set = TRUE;

    return argc - 1;
}

/* ImageFormat takes a string and a double. */

int imgformatfn(int argc, char **argv, int i, PARM *p)
{
    double d;
    char *s;
    IMGFORMAT *imf;

    if (argc > 1) {
        imf = new(IMGFORMAT);
        imf->fmt = copy_string(argv[i]);
        if (strcasecmp(imf->fmt, "TLRC") == 0) {    // abbr not allowed
            if (argc == 3) {
                d = strtod(argv[i+1], &s);
                if (*s != '\0') {
                    fatalerr("%s: badly formed number '%s'", argv[i-1], argv[i+1]);
                }
                imf->res = d;
            } else {
                fatalerr("image resolution required for ImageFormat TLRC");
            }
        } else if (strcasecmp(imf->fmt, "ORIG") != 0 && strcasecmp(imf->fmt, "ORTHO") != 0) {
            fatalerr("ImageFormat must be ORIG or TLRC");
        }
        p->ptr = imf;
        p->set = TRUE;
    }

    return argc - 1;
}

/* Mu value, for regularization. */

int mufn(int argc, char **argv, int i, PARM *p)
{
    double d;
    char *s, *t;
    MUINFO *m;

    REQ1();

    m = new(MUINFO);
    m->op = ADD_MU;             /* default is add */

    /* See if there is an initial '*' or '+'. */

    s = argv[i];
    if (*s == '*') {
        m->op = MPY_MU;
        s++;
    } else if (*s == '+') {
        m->op = ADD_MU;
        s++;
    }

    d = strtod(s, &t);
    if (*t != '\0') {
        fatalerr("%s: badly formed number '%s'", argv[i-1], argv[i]);
    }

    m->mu = d;
    p->ptr = m;
    p->set = TRUE;

    return 1;
}

/* Model */

int modelfn(int argc, char **argv, int i, PARM *p)
{
    int n;
    int j = argc - 1;
    char *arg;
    MODELINFO *m;

    arg = argv[i];
    n = strlen(arg);

    m = new(MODELINFO);

    if (strncasecmp(arg, "SingleSphere", n) == 0) {
        if (argc != 5) {
            fatalerr("improper sphere specification");
        }
        m->model = SSPHERE;
        m->sphere[0] = strtod(argv[i+1], NULL); // @@@ should add error checks
        m->sphere[1] = strtod(argv[i+2], NULL);
        m->sphere[2] = strtod(argv[i+3], NULL);
    } else if (strncasecmp(arg, "MultiSphere", n) == 0) {
        m->model = MSPHERE;
    } else if (strncasecmp(arg, "Nolte", n) == 0) {
        m->model = NOLTE;
        m->order = -1;
        if (argc >= 3) {
            m->order = strtod(argv[i+1], NULL);
        }
    } else {
        fatalerr("Model must be SingleSphere x y z (cm), MultiSphere, or Nolte");
    }

    p->ptr = m;
    p->set = TRUE;

    return j;
}

int imgmetricfn(int argc, char **argv, int i, PARM *p)
{
    int n;
    int j = argc - 1;
    char *arg;
    METRICINFO *m;

    arg = argv[i];
    n = strlen(arg);

    m = new(METRICINFO);

    if (strncasecmp(arg, "Signal", n) == 0 || strncasecmp(arg, "Moment", n) == 0) {
        m->metric = MOMENT;
    } else if (strncasecmp(arg, "Power", n) == 0) {
        m->metric = POWER;
    } else if (strncasecmp(arg, "Hilbert", n) == 0) {
        m->metric = HILBERT;
    } else if (strncasecmp(arg, "Kurtosis", n) == 0) {
        m->metric = KURTOSIS;
    } else if (strncasecmp(arg, "MutualInfo", n) == 0) {
        if (argc != 4) {
            fatalerr("%s takes 2 arguments", argv[i]);
        }
        m->metric = MUT_INFO;
        m->lags = strtod(argv[i+1], NULL); // @@@ should add error checks
        m->dims = atoi(argv[i+2]);
// #if 0
    } else if (strncasecmp(arg, "RankVectorEntropy", n) == 0) {
        if (argc != 4) {
            fatalerr("%s takes 2 arguments", argv[i]);
        }
        m->metric = RV_ENTROPY;
        m->tau = strtod(argv[i+1], NULL); // @@@ should add error checks
        m->dims = atoi(argv[i+2]);
    } else if (strncasecmp(arg, "SpectEntropy", n) == 0) {
        m->metric = S_ENTROPY;
    } else if (strncasecmp(arg, "TransferEntropy", n) == 0) {
        if (argc != 4) {
            fatalerr("%s takes 2 arguments", argv[i]);
        }
        m->metric = ST_ENTROPY;
        m->lags = strtod(argv[i+1], NULL); // @@@ should add error checks
        m->dims = atoi(argv[i+2]);
    } else if (strncasecmp(arg, "ConditionalEntropy", n) == 0) {
        if (argc != 5) {
            fatalerr("%s takes 3 arguments", argv[i]);
        }
        m->metric = RVC_ENTROPY;
        m->tau = strtod(argv[i+1], NULL); // @@@ should add error checks
        m->lags = strtod(argv[i+2], NULL);
        m->dims = atoi(argv[i+3]);
// #endif
    } else {
        msg("unknown image metric '%s'\n", argv[i]);
        fatalerr("metrics are: Signal, Power, Hilbert, Kurtosis, MutualInfo, or TransferEntropy\n");
//        fatalerr(""RankVectorEntropy, SpectEntropy, TransferEntropy, or ConditionalEntropy");
    }

    p->ptr = m;
    p->set = TRUE;

    return j;
}

/* Fill in a PARMINFO structure from the PARM list. */

int get_params(PARMINFO *Params)
{
    /* Initialize the Params struct with default values first. */

    new_params(Params);
    return do_get_params(Params);
}

int do_get_params(PARMINFO *Params)
{
    int i, j;
    double d;
    double *v;
    char *key, *s, *t;
    PARM *p;
    PARM_TABLE *param;
    MARKINFO *mark;
    IMGFORMAT *imf;
    MUINFO *mu;
    MODELINFO *model;
    METRICINFO *metric;

    /* We already know how many markers there are. */

    Params->NumMark = i = Parm_list_nmarkers;
    if (i > 0) {
        Params->Marker = new_array(MARKINFO, i);
        for (j = 0; j < i; j++) {
            Params->Marker[j].MarkName = NULL;
        }
    }

    /* Check for errors as we go through the PARM list. */

    for (p = Parm_list; p; p = p->next) {
        if (p->set) {
            param = p->param;
            key = param->key;
            if (strcmp(key, "param") == 0) {
                s = (char *)p->ptr;

                /* Copy parameter file basename without .param for use in filenames. */

                t = rindex(s, '/');
                if (t != NULL) {
                    s = t + 1;
                }
                t = rindex(s, '.');
                if (t != NULL && strcmp(t, ".param") == 0) {
                    *t = '\0';
                }
                Params->ParmName = s;
            } else if (strcmp(key, "DataSet") == 0) {
                s = (char *)p->ptr;

                /* Remove any trailing slash. */

                t = rindex(s, '/');
                if (t != NULL && t[1] == '\0') {
                    *t = '\0';
                }

                if (!direxists(s)) {

                    /* Try again with .ds added. */

                    t = new_string(strlen(s) + 3);
                    strcpy(strecpy(t, s), ".ds");
                    if (!direxists(t)) {
                        fatalerr("can't access dataset '%s'", s);
                    }
                    s = t;
                    p->ptr = t;
                }
                Params->DataSetName = s;
#if BTI
            } else if (strcmp(key, "PDFName") == 0) {
                s = (char *)p->ptr;
                if (!fileexists(s)) {
                    fatalerr("can't access pdf file '%s'", s);
                }
#endif
            } else if (strcmp(key, "Marker") == 0 || strcmp(key, "SegFile") == 0) {
                mark = (MARKINFO *)p->ptr;
                i = p->mnum - 1;    // mnum is IO 1
                Params->Marker[i].MarkName = mark->MarkName;
                Params->Marker[i].MarkName2 = mark->MarkName2;
                Params->Marker[i].MarkStart = mark->MarkStart;
                Params->Marker[i].MarkEnd = mark->MarkEnd;
                Params->Marker[i].Sum = mark->Sum;
                Params->Marker[i].SegFile = mark->SegFile;
                Params->Marker[i].FileName = mark->FileName;
            } else if (strcmp(key, "Baseline") == 0) {
                v = (double *)p->ptr;
                if (v[0] >= v[1]) {
                    fatalerr("Baseline window start time must be less than end time");
                }
                Params->BaseStart = v[0];
                Params->BaseEnd = v[1];
            } else if (strcmp(key, "DataSegment") == 0) {
                v = (double *)p->ptr;
                if (v[0] >= v[1]) {
                    fatalerr("DataSegment start time must be less than end time");
                }
                Params->DataSegStart = v[0];
                Params->DataSegEnd = v[1];
            } else if (strcmp(key, "SignSegment") == 0) {
                v = (double *)p->ptr;
                if (v[0] >= v[1]) {
                    fatalerr("SignSegment start time must be less than end time");
                }
                Params->SignSegStart = v[0];
                Params->SignSegEnd = v[1];
            } else if (strcmp(key, "FilterType") == 0) {
                s = (char *)p->ptr;
                if (strcmp(s, "IIR") == 0) {
                    Params->FilterType = IIR;
                } else if (strcmp(s, "FFT") == 0) {
                    Params->FilterType = FFT;
                } else {
                    fatalerr("FilterType must be either 'FFT' or 'IIR'");
                }
            } else if (strcmp(key, "CovBand") == 0) {
                v = (double *)p->ptr;
                if (v[0] >= v[1]) {
                    fatalerr("CovBand lowpass frequency must be greater than highpass frequency");
                }
                Params->CovHP = v[0];
                Params->CovLP = v[1];
            } else if (strcmp(key, "ImageBand") == 0) {
                v = (double *)p->ptr;
                if (v[0] >= v[1]) {
                    fatalerr("ImageBand lowpass frequency must be greater than highpass frequency");
                }
                Params->ImageHP = v[0];
                Params->ImageLP = v[1];
            } else if (strcmp(key, "OrientBand") == 0) {
                v = (double *)p->ptr;
                if (v[0] >= v[1]) {
                    fatalerr("OrientBand lowpass frequency must be greater than highpass frequency");
                }
                Params->OrientHP = v[0];
                Params->OrientLP = v[1];
            } else if (strcmp(key, "NoiseBand") == 0) {
                v = (double *)p->ptr;
                if (v[0] >= v[1]) {
                    fatalerr("NoiseBand lowpass frequency must be greater than highpass frequency");
                }
                Params->NoiseHP = v[0];
                Params->NoiseLP = v[1];
            } else if (strcmp(key, "SmoothBand") == 0) {
                v = (double *)p->ptr;
                if (v[0] >= v[1]) {
                    fatalerr("SmoothBand lowpass frequency must be greater than highpass frequency");
                }
                Params->SmoothHP = v[0];
                Params->SmoothLP = v[1];
            } else if (strcmp(key, "Notch") == 0) {
                Params->Notch = TRUE;
            } else if (strcmp(key, "Hz") == 0) {
                d = *(double *)p->ptr;
                Params->Hz = d;
            } else if (strcmp(key, "XBounds") == 0) {
                v = (double *)p->ptr;
                if (v[0] > v[1]) {
                    fatalerr("starting XBound must be <= ending XBound");
                }
                Params->SAMStart[X_] = .01 * v[0];
                Params->SAMEnd[X_] = .01 * v[1];
            } else if (strcmp(key, "YBounds") == 0) {
                v = (double *)p->ptr;
                if (v[0] > v[1]) {
                    fatalerr("starting YBound must be <= ending YBound");
                }
                Params->SAMStart[Y_] = .01 * v[0];
                Params->SAMEnd[Y_] = .01 * v[1];
            } else if (strcmp(key, "ZBounds") == 0) {
                v = (double *)p->ptr;
                if (v[0] > v[1]) {
                    fatalerr("starting ZBound must be <= ending ZBound");
                }
                Params->SAMStart[Z_] = .01 * v[0];
                Params->SAMEnd[Z_] = .01 * v[1];
            } else if (strcmp(key, "ImageStep") == 0) {
                Params->SAMStep = .01 * *(double *)p->ptr;
            } else if (strcmp(key, "TargetName") == 0) {
                s = (char *)p->ptr;
                Params->Target = s;
                t = rindex(s, '/');
                if (t != NULL) {
                    s = t + 1;
                }
                Params->TargetName = s;
            } else if (strcmp(key, "ROIList") == 0) {
                s = (char *)p->ptr;
                Params->ROIList = s;
                t = rindex(s, '/');
                if (t != NULL) {
                    s = t + 1;
                }
                Params->ROIName = s;
            } else if (strcmp(key, "Extent") == 0) {
                d = *(double *)p->ptr;
                Params->Extent = d * .001;  // convert mm to m
            } else if (strcmp(key, "Interpolate") == 0) {
                d = *(double *)p->ptr;
                Params->Interpolate = d;
            } else if (strcmp(key, "AtlasName") == 0) {
                Params->Atlas = s = (char *)p->ptr;
                // get the basename, without extension, for aux filenames
                t = rindex(s, '/');
                if (t != NULL) {
                    s = t + 1;
                }
                t = rindex(s, '.');
                if (t != NULL) {
                    i = t - s;
                } else {
                    i = strlen(s);
                }
                Params->AtlasName = t = new_array(char, i + 1);
                strncpy(t, s, i);
                t[i] = '\0';
            } else if (strcmp(key, "MRIDirectory") == 0) {
                s = (char *)p->ptr;
                if (access(s, F_OK) != 0) {
                    fatalerr("can't access MRIDirectory '%s'", s);
                }
                Params->MRIDirectory = s;
            } else if (strcmp(key, "ImageDirectory") == 0) {
                Params->DirName = (char *)p->ptr;
            } else if (strcmp(key, "MaskName") == 0) {
                Params->Mask = (char *)p->ptr;
            } else if (strcmp(key, "HullName") == 0) {
                Params->HullName = (char *)p->ptr;
            } else if (strcmp(key, "MRIPattern") == 0) {
                Params->MRIPattern = (char *)p->ptr;
            } else if (strcmp(key, "PrefixLength") == 0) {
                if (!p->mnum) { // use a delimiter character
                    Params->PrefixChar = *(char *)p->ptr;
                    Params->NumPrefix = -1;
                } else {
                    Params->NumPrefix = *(int *)p->ptr;
                }
            } else if (strcmp(key, "ImageFormat") == 0) {
                imf = (IMGFORMAT *)p->ptr;
                s = imf->fmt;
                i = strlen(s);
                if (strncasecmp(s, "ORIG", i) == 0 || strncasecmp(s, "ORTHO", i) == 0) {
                    Params->ImageFormat = ORIG;
                } else if (strncasecmp(s, "TLRC", i) == 0) {
                    Params->ImageFormat = TLRC;
                    Params->ImageRes = imf->res;
                }
            } else if (strcmp(key, "CovType") == 0) {
                s = (char *)p->ptr;
                i = strlen(s);
                if (strncasecmp(s, "GLOBAL", i) == 0) {
                    Params->CovType = GLOBAL_;
                } else if (strncasecmp(s, "SUM", i) == 0) {
                    Params->CovType = SUM_;
                } else if (strncasecmp(s, "ALL", i) == 0) {
                    Params->CovType = ALL_;
                } else {
                    fatalerr("CovType must be GLOBAL, SUM, or ALL");
                }
            } else if (strcmp(key, "Mu") == 0) {
                mu = (MUINFO *)p->ptr;
                Params->Mu = mu->mu;
                Params->Operator = mu->op;
            } else if (strcmp(key, "Model") == 0) {
                model = (MODELINFO *)p->ptr;
                Params->Model = model->model;
                if (model->model == SSPHERE) {
                    Params->Sp[0] = model->sphere[0] * .01;
                    Params->Sp[1] = model->sphere[1] * .01;
                    Params->Sp[2] = model->sphere[2] * .01;
                } else if (model->model == NOLTE) {
                    if (model->order > 0) {
                        Params->Order = model->order;
                    }
                }
            } else if (strcmp(key, "ImageMetric") == 0) {
                metric = (METRICINFO *)p->ptr;
                Params->ImageMetric = metric->metric;
                if (metric->metric == MUT_INFO) {
                    Params->Lags = metric->lags;
                    Params->Dims = metric->dims;
// #if 0
                } else if (metric->metric == RV_ENTROPY) {
                    Params->Tau = metric->tau;
                    if (metric->tau < 0.) {
                        fatalerr("RVE Tau must be >= 0");
                    }
                } else if (metric->metric == ST_ENTROPY) {
                    Params->Lags = metric->lags;
                    Params->Dims = metric->dims;
                } else if (metric->metric == RVC_ENTROPY) {
                    Params->Tau = metric->tau;
                    Params->Lags = metric->lags;
                    Params->Dims = metric->dims;
// #endif
                }
            } else if (strcmp(key, "Order") == 0) {
                i = *(int *)p->ptr;
                Params->Order = i;
            } else if (strcmp(key, "NumVerts") == 0) {
                i = *(int *)p->ptr;
                Params->MaxVertex = i;
            } else if (strcmp(key, "MinZ") == 0) {
                d = *(double *)p->ptr;
                Params->MinZ = d * .01;    // convert cm to m
            } else if (strcmp(key, "MinSpan") == 0) {
                d = *(double *)p->ptr;
                Params->MinSpan = d;
            } else if (strcmp(key, "MaxAngle") == 0) {
                d = *(double *)p->ptr;
                Params->MaxAngle = d;
            } else if (strcmp(key, "MaxDot") == 0) {
                d = *(double *)p->ptr;
                if (d < 0. || d > 1.) {
                    fatalerr("MaxDot threshold must be between 0 and 1");
                }
                Params->MaxDot = d;
            } else if (strcmp(key, "MinSNR") == 0) {
                d = *(double *)p->ptr;
                if (d < 0. || d > 1.) {
                    fatalerr("MinSNR must be between 0 and 1");
                }
                Params->MinSNR = d;
            } else if (strcmp(key, "DeltaTranslate") == 0) {
                d = *(double *)p->ptr;
                Params->DeltaTrans = d * .001;          // convert mm to metres
            } else if (strcmp(key, "DeltaRotate") == 0) {
                d = *(double *)p->ptr;
                Params->DeltaRotat = d * M_PI / 180.;   // convert degrees to radians
            } else if (strcmp(key, "InitDisplacement") == 0) {
                d = *(double *)p->ptr;
                Params->InitDisp = d * .001;            // convert mm/mrad to m/rad
            } else if (strcmp(key, "Damping") == 0) {   // sequential line-search damping
                d = *(double *)p->ptr;
                if (d <= 0. || d > 1) {
                    fatalerr("Damping must be > 0 and <= 1");
                }
                Params->Damping = d;
            } else if (strcmp(key, "LineSteps") == 0) { // sequential line-step steps (+/-)
                i = *(int *)p->ptr;
                if (i < 2) {
                    fatalerr("LineSteps must be > 1");
                }
                Params->LineSteps = i;
            } else if (strcmp(key, "EndDisplacement") == 0) {
                d = *(double *)p->ptr;
                Params->EndDisp = d * .001;             // convert mm/mrad to m/rad
            } else if (strcmp(key, "TimeStep") == 0) {
                d = *(double *)p->ptr;
                Params->TimeStep = d;
            } else if (strcmp(key, "TimeInt") == 0) {
                d = *(double *)p->ptr;
                Params->TimeInt = d;
            } else if (strcmp(key, "RemoveBaseline") == 0) {
                s = (char *)p->ptr;
                i = strlen(s);
                if (strncasecmp(s, "NONE", i) == 0) {
                    Params->RemoveBaseline = NONE;
                } else if (strncasecmp(s, "BYVOXEL", i) == 0) {
                    Params->RemoveBaseline = BYVOXEL;
                } else if (strncasecmp(s, "GLOBAL", i) == 0) {
                    Params->RemoveBaseline = GLOBAL;
                } else {
                    fatalerr("RemoveBaseline must be NONE, BYVOXEL, or GLOBAL");
                }
            } else if (strcmp(key, "BlockName") == 0) {
                s = (char *)p->ptr;
                Params->BlockName = s;
            }
        }
    }

    // If the prefix was specified using a delimiter, find it in the dataset name.

    s = Params->DataSetName;
    p = check_get_parm("PrefixLength");
    if (s && p && Params->NumPrefix < 0) {
        t = rindex(s, '/');
        if (t) {
            s = t + 1;
        }
        t = index(s, Params->PrefixChar);
        if (t == NULL) {
            // not found, use entire setname
            t = index(s, '.');
        }
        Params->NumPrefix = t - s;      // set the length
    }

    // ImageBand and OrientBand default to CovBand

    if (Params->ImageHP == -999.) {
        Params->ImageHP = Params->CovHP;
        Params->ImageLP = Params->CovLP;
    }

    if (Params->OrientHP == -999.) {
        Params->OrientHP = Params->CovHP;
        Params->OrientLP = Params->CovLP;
    }

    return 0;
}

/* Record the current set of parameters in a file. */

void log_params(char *dirname)
{
    int i;
    double d, *v;
    char *key, *name, *s;
    PARM *p;
    PARM_TABLE *param;
    MARKINFO *mark;
    IMGFORMAT *imf;
    MUINFO *mu;
    MODELINFO *model;
    METRICINFO *metric;
    FILE *f;
    time_t t;
    struct tm *tmp;
    char buf[100];

    name = new_string(strlen(dirname) + strlen(Progname) + 7);
    s = strecpy(name, dirname);
    s = strecpy(s, "/");
    s = strecpy(s, Progname);
    s = strecpy(s, ".param");
    f = fileopen(name, "a");
    free(name);

    fprintf(f, "\n# %s\n", Usage);
    t = time(NULL);
    tmp = localtime(&t);
    if (tmp && strftime(buf, sizeof(buf), "%Y-%m-%d %T", tmp)) {
        fprintf(f, "# %s\n", buf);
    }

    for (p = Parm_list; p; p = p->next) {
        if (p->set) {
            param = p->param;
            key = param->key;
            if (strcmp(key, "Marker") == 0) {
                mark = (MARKINFO *)p->ptr;
                fprintf(f, "Marker ");
                fprintf(f, "%s ", mark->MarkName);
                fprintf(f, "%g %g ", mark->MarkStart, mark->MarkEnd);
                fprintf(f, "%s", mark->Sum ? "True" : "False");
                if (mark->MarkName2) {
                    fprintf(f, " %s", mark->MarkName2);
                }
                fprintf(f, "\n");
            } else if (strcmp(key, "SegFile") == 0) {
                mark = (MARKINFO *)p->ptr;
                fprintf(f, "SegFile ");
                fprintf(f, "%s ", mark->MarkName);
                fprintf(f, "%s ", mark->FileName);
                fprintf(f, "%s", mark->Sum ? "True" : "False");
                fprintf(f, "\n");
            } else if (strcmp(key, "ImageFormat") == 0) {
                imf = (IMGFORMAT *)p->ptr;
                s = imf->fmt;
                fprintf(f, "%s %s", key, s);
                i = strlen(s);
                if (strncasecmp(s, "TLRC", i) == 0) {
                    fprintf(f, " %g", imf->res);
                }
                fprintf(f, "\n");
            } else if (strcmp(key, "Mu") == 0) {
                mu = (MUINFO *)p->ptr;
                fprintf(f, "%s %s%g\n", key, mu->op == MPY_MU ? "*" : "", mu->mu);
            } else if (strcmp(key, "PropMu") == 0) {
                mu = (MUINFO *)p->ptr;
                fprintf(f, "%s %g\n", key, mu->mu);
            } else if (strcmp(key, "Model") == 0) {
                model = (MODELINFO *)p->ptr;
                if (model->model == SSPHERE) {
                    fprintf(f, "%s SingleSphere %g %g %g\n", key,
                            model->sphere[0], model->sphere[1], model->sphere[2]);
                } else if (model->model == MSPHERE) {
                    fprintf(f, "%s MultiSphere\n", key);
                } else {
                    fprintf(f, "%s Nolte", key);
                    if (model->order > 0) {
                        fprintf(f, " %d", model->order);
                    }
                    fprintf(f, "\n");
                }
            } else if (strcmp(key, "ImageMetric") == 0) {
                metric = (METRICINFO *)p->ptr;
                if (metric->metric == MOMENT) {
                    fprintf(f, "%s Signal\n", key);
                } else if (metric->metric == POWER) {
                    fprintf(f, "%s Power\n", key);
                } else if (metric->metric == HILBERT) {
                    fprintf(f, "%s Hilbert\n", key);
// #if 0
                } else if (metric->metric == RV_ENTROPY) {
                    fprintf(f, "%s RankVectorEntropy %g %d\n", key, metric->tau, metric->dims);
                } else if (metric->metric == S_ENTROPY) {
                    fprintf(f, "%s SpectEntropy\n", key);
                } else if (metric->metric == ST_ENTROPY) {
                    fprintf(f, "%s TransferEntropy %g %d\n", key, metric->lags, metric->dims);
                } else if (metric->metric == RVC_ENTROPY) {
                    fprintf(f, "%s ConditionalEntropy %g %g %d\n", key, metric->tau, metric->lags, metric->dims);
// #endif
                } else if (metric->metric == KURTOSIS) {
                    fprintf(f, "%s Kurtosis\n", key);
                } else if (metric->metric == MUT_INFO) {
                    fprintf(f, "%s MutualInfo %g %d\n", key, metric->lags, metric->dims);
                } else {
                    fprintf(f, "%s unknown\n", key);
                }
            } else if (strcmp(key, "MinZ") == 0) {
                d = *(double *)p->ptr;
                fprintf(f, "%s %g\n", key, d * 100.);
            } else if (strcmp(key, "InitDisplacement") == 0) {
                d = *(double *)p->ptr;
                fprintf(f, "%s %g\n", key, d * 100.);
            } else if (strcmp(key, "EndDisplacement") == 0) {
                d = *(double *)p->ptr;
                fprintf(f, "%s %g\n", key, d * 100.);
            } else if (param->parsefn == stringfn) {
                s = (char *)p->ptr;
                fprintf(f, "%s %s\n", key, s);
            } else if (param->parsefn == intfn) {
                i = *(int *)p->ptr;
                fprintf(f, "%s %d\n", key, i);
            } else if (param->parsefn == doublefn) {
                d = *(double *)p->ptr;
                fprintf(f, "%s %g\n", key, d);
            } else if (param->parsefn == double2fn) {
                v = (double *)p->ptr;
                fprintf(f, "%s %g %g\n", key, v[0], v[1]);
            } else if (param->parsefn == double3fn) {
                v = (double *)p->ptr;
                fprintf(f, "%s %g %g %g\n", key, v[0], v[1], v[2]);
            } else if (param->parsefn == boolfn) {
                fprintf(f, "%s\n", key);
            } else if (param->parsefn == prefixfn) {
                if (p->mnum == 0) {
                    fprintf(f, "%s %c\n", key, *(char *)p->ptr);
                } else {
                    fprintf(f, "%s %d\n", key, *(int *)p->ptr);
                }
            }
        }
    }

    fclose(f);
}
