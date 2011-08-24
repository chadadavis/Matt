/* util.h -- A couple miscellaneous string-related utility functions.
 *
 * Author: Matt Menke (2007)
 *
 * Matt is licensed under the GNU public license version 2.0.
 *
 * If you would like to license Matt in an enviroment where the GNU public license is
 * unacceptable (such as inclusion in a non-GPL software package) comercial Matt
 * licensing is available through the MIT and Tufts offices of Technology Transfer.
 * Contact betawrap@csail.mit.edu or cowen@cs.tufts.edu for more information.
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(WIN32) || defined(__WIN32__)
#define DIR_DELIMITER '\\'
#else
#define DIR_DELIMITER '/'
#endif

/* Needed in a couple situations to let me use fastcall calling convention with MSCV */
#if defined(WIN32)
#define CDECL __cdecl
#else
#define CDECL
#endif

#define UCASE(x) ((x)>='a' && (x) <= 'z' ? (x)-'a'+'A' : (x))
#define LCASE(x) ((x)>='A' && (x) <= 'Z' ? (x)-'A'+'a' : (x))

/* Note:  s must be at least of length len, not including
 * null terminator.  Returns 0 on failure.
 */
int ReadDouble(char *s, int len, double *out);
int ReadFloat(char *s, int len, float *out);
int ReadInt(char *s, int len, int *out);

void *memdup (const void *mem, int size);

/* Checks if all of s1 matches the first strlen(s1) characters
 * of s2, and is not case sensitive.  If there were a case-insensitive
 * form of strstr, would just be str(s2,s1)==s1.
 */
int stristart(const char *s1, const char *s2);

#define WHITESPACE(x) ((x)==' ' ||(x)=='\t' ||(x)=='\n' ||(x)=='\r')

/* Used to cleanup old status display line and create a new one. */
void DisplayProgress(int *backspace, char *progressTemplate, int count);
void CleanupAndBackspace(int len);

#endif
