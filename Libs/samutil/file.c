/* File utilities. */

#include "samutil.h"
#include <ctype.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* Open a file or exit on failure. */

FILE *fileopen(char *name, char *mode)
{
    FILE *f;

    if ((f = fopen(name, mode)) == NULL) {
        fatalerr("can't %s '%s'", *mode == 'r' ? "open" : "write", name);
    }
    return f;
}

/* Check to see if a file exists before we open it. */

int fileexists(char *name)
{
    struct stat buf;

    if (stat(name, &buf) < 0 || (buf.st_mode & S_IFMT) != S_IFREG) {
        return FALSE;
    }
    return TRUE;
}

/* Check to see if a directory exists. */

int direxists(char *name)
{
    struct stat buf;

    if (stat(name, &buf) < 0 || (buf.st_mode & S_IFMT) != S_IFDIR) {
        return FALSE;
    }
    return TRUE;
}

/* Create a directory and any missing parents. */

static int is_dir_separator(char c)
{
#ifdef _WIN32
    return c == '/' || c == '\\';
#else
    return c == '/';
#endif
}

int makedirs(char *name)
{
    char *copy, *s, *start;
    char separator;
    int saved_errno;

    if (name == NULL || *name == '\0') {
        errno = EINVAL;
        return -1;
    }
    if (direxists(name)) {
        return 0;
    }

    copy = copy_string(name);
    start = copy;
#ifdef _WIN32
    if (isalpha((unsigned char)copy[0]) && copy[1] == ':') {
        start = copy + 2;
    } else if (is_dir_separator(copy[0]) && is_dir_separator(copy[1])) {
        start = copy + 2;
        while (*start && !is_dir_separator(*start))
            start++;
        if (*start)
            start++;
        while (*start && !is_dir_separator(*start))
            start++;
    }
#endif

    for (s = start; *s; s++) {
        if (!is_dir_separator(*s)) {
            continue;
        }
        if (s == copy || s[-1] == ':' || is_dir_separator(s[-1])) {
            continue;
        }
        separator = *s;
        *s = '\0';
        if (mkdir(copy, 0777) == -1 && errno != EEXIST) {
            saved_errno = errno;
            free(copy);
            errno = saved_errno;
            return -1;
        }
        if (!direxists(copy)) {
            free(copy);
            errno = ENOTDIR;
            return -1;
        }
        *s = separator;
    }

    if (mkdir(copy, 0777) == -1 && errno != EEXIST) {
        saved_errno = errno;
        free(copy);
        errno = saved_errno;
        return -1;
    }
    if (!direxists(copy)) {
        free(copy);
        errno = ENOTDIR;
        return -1;
    }
    free(copy);
    return 0;
}

/* Extract the exact multi-character SAM directory options used by legacy
 * getopt-based programs, removing them before getopt sees their prefixes. */

void parse_samdir_args(int *argc, char **argv, char **input, char **output)
{
    char **target;
    int i, j;

    for (i = 1; i < *argc;) {
        target = NULL;
        if (strcmp(argv[i], "-i_SAMdir") == 0) {
            target = input;
        } else if (strcmp(argv[i], "-o_SAMdir") == 0) {
            target = output;
        }
        if (target == NULL) {
            i++;
            continue;
        }
        if (i + 1 >= *argc) {
            fatalerr("%s requires an argument", argv[i]);
        }
        *target = argv[i + 1];
        for (j = i; j + 2 <= *argc; j++) {
            argv[j] = argv[j + 2];
        }
        *argc -= 2;
    }
}

/* Like fgets() but remove the trailing newline. Return a pointer to
the nul at the end of the string or NULL on error. */

char *fgetline(char *buf, int maxlen, FILE *infile)
{
    int i, c;
    char *s;

    s = buf;
    i = maxlen - 1;
    if (i < 0) {
        return NULL;
    }
    if (i == 0) {
        *s = '\0';
        return s;
    }
    c = getc(infile);
    if (c == EOF) {
        return NULL;
    }
    while (c != '\n' && c != EOF) {
        *s++ = c;
        if (--i <= 0) {
            break;
        }
        c = getc(infile);
    }
    *s = '\0';
    return s;
}
