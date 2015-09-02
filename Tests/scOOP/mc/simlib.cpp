#include "simlib.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

void trim (char *str) {
    ltrim (str);
    rtrim (str);
}


void ltrim (char *str) {
    char *tr;
    int c;

    if (!str) return;
    tr = strdup (str);
    c  = 0;
    while ((tr[c] == ' ') || (tr[c] == '\n') || (tr[c] == '\t') || (tr[c] == '\r')) c++;
    strcpy (str,tr+c);
    free (tr);
}




void rtrim (char *str) {
    int nul;

    if (!str) return;
    nul = strlen(str)-1;
    while ((nul > 0) && ((str[nul] == ' ') || (str[nul] == '\t') || (str[nul] == '\n')  || (str[nul] == '\r') )) {
        str[nul] = '\0';
        nul--;
    }
}




void strip_comment (char *line) {
    char *c;

    if (!line) return;
    /* search for a comment mark and replace it by a zero */
    if ((c = strchr(line,COMMENTSIGN)) != NULL) (*c) = 0;
}




int continuing(char *s) {
    int sl;

    rtrim(s);
    sl = strlen(s);
    if ((sl > 0) && (s[sl-1] == CONTINUE)) {
        s[sl-1] = 0;
        return 1; /*true*/
    } else return 0; /*false*/
}




void upstring (char *str) {
    int i;
    for (i=0; (i < (int)strlen(str)); i++) str[i] = toupper(str[i]);
}




void beforecommand(char *str,char *pline,char commandc) {
    char *dummy;

    strcpy(str,pline);
    if ((dummy = strchr (str,commandc)) != NULL) (*dummy) = 0;
    trim (str);
}




void aftercommand(char *str, char *pline,char commandc) {
    char *dummy;
    int i;

    strcpy(str,pline);
    if ((dummy = strchr (str,commandc)) != NULL) {
        i=0;
        while( (*dummy) != str[i]) {
            str[i] = ' ';
            i++;
        }
        str[i] = ' ';
    }
    trim (str);
}




char *fgets2(char *line, int n, FILE *stream) {
    char *c;

    if (fgets(line,n,stream)==NULL) {
        return NULL;
    }
    if ((c=strchr(line,'\n'))!=NULL)
        *c=0;

    return line;
} // issues with compiling, from IOoperations


int longarrayCpy (long **target, long **source, long targetlength, long sourcelength) {

    /*if ( (*target) != NULL) {
        if ( targetlength == sourcelength ) {
            memcpy((*target),(*source), sizeof(long)*(sourcelength));
            return 0;
        } else {
            free(*target);
        }
    }*/
    if ( (*target) != NULL)
        (*target) = (long*) realloc((*target), sizeof(long)*(sourcelength));
    else
        (*target) = (long int*) malloc( sizeof(long)*(sourcelength));
    memcpy((*target),(*source), sizeof(long)*(sourcelength));

    return 0;
}



