/** @file simlib.h*/

#ifndef SIMLIB_H
#define SIMLIB_H

#include "../structures/macros.h"

#include <stdio.h>

/**
 * @brief trim string from left and right
 * @param str
 */
void trim (char *str);

/**
 * @brief trim string from left
 * @param str
 */
void ltrim (char *str);

/**
 * @brief trim string from right
 * @param str
 */
void rtrim (char *str);

/**
 * @brief removes comments
 * @param line
 */
void strip_comment (char *line);

/**
 * @brief test is there is still something left in string
 * @param s
 * @return
 */
int continuing(char *s);

/**
 * @brief make string uppercase
 * @param str
 */
void upstring (char *str);

/**
 * @brief return string that goes before comand character
 * @param str
 * @param pline
 * @param commandc
 */
void beforecommand(char *str,char *pline,char commandc);

/**
 * @brief return string that goes after command character
 * @param str
 * @param pline
 * @param commandc
 */
void aftercommand(char *str, char *pline,char commandc);

/**
 * @brief reads a string from stream of max length n
 * @param line
 * @param n
 * @param stream
 * @return
 */
char *fgets2(char *line, int n, FILE *stream);

/**
 * @brief longarray_cpy
 * @param target
 * @param source
 * @param targetlength
 * @param sourcelength
 * @return
 */
int longarrayCpy (long **target, long **source, long targetlength, long sourcelength);



#endif // SIMLIB_H
