/** @file mygetline.cpp*/

#include "mygetline.h"

char * fgetln(FILE *fp, size_t *len) {
    static char *buf = NULL;
    static size_t bufsiz = 0;
    char *ptr;


    if (buf == NULL) {
        bufsiz = BUFSIZ;
        if ((buf = (char*) malloc(bufsiz)) == NULL)
            return NULL;
    }

    if (fgets(buf, bufsiz, fp) == NULL)
        return NULL;

    *len = 0;
    while ((ptr = strchr(&buf[*len], '\n')) == NULL) {
        size_t nbufsiz = bufsiz + BUFSIZ;
        char *nbuf = (char*) realloc(buf, nbufsiz);

        if (nbuf == NULL) {
            int oerrno = errno;
            free(buf);
            errno = oerrno;
            buf = NULL;
            return NULL;
        } else
            buf = nbuf;

        if (fgets(&buf[bufsiz], BUFSIZ, fp) == NULL) {
            buf[bufsiz] = '\0';
            *len = strlen(buf);
            return buf;
        }

        *len = bufsiz;
        bufsiz = nbufsiz;
    }

    *len = (ptr - buf) + 1;
    return buf;
}
