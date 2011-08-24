/* util.c -- A couple miscellaneous string-related utility functions.
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

#include "util.h"
#include "FileReader.h"


int ReadDouble(char *s, int len, double *d) {
	int success;
	char c = s[len];
	s[len] = 0;
	success = (1 == sscanf(s, "%lf", d));
	if (!success) *d = 0;
	s[len] = c;
	return success;
}

int ReadFloat(char *s, int len, float *f) {
	int success;
	char c = s[len];
	s[len] = 0;
	success = (1 == sscanf(s, "%f", f));
	if (!success) *f = 0;
	s[len] = c;
	return success;
}

void *memdup (const void *mem, int size) {
	void *out = malloc(size);
	if (!out) return 0;
	memcpy(out, mem, size);
	return out;
}

char ** parseArgs(char *s, int *argc, char split) {
	int args=0;
	int i=0;
	char **argv = (char**) malloc(sizeof(char*) * 32);
	if (s[0]) {
		argv[0]=s;
		args++;
	}
	while (s[i]) {
		if (s[i]==split) {
			s[i]=0;
			if (args%32==0) {
				argv = (char**) realloc(argv, sizeof(char*) * (args+32));
			}
			argv[args]=&s[i+1];
			args++;
		}
		i++;
	}
	*argc = args;
	return argv;
}

int stristart(const char *s1, const char *s2) {
	int i=0;
	while (s1[i] && UCASE(s1[i]) == UCASE(s2[i])) i++;
	if (s1[i]==0) {
		return 0;
	}
	else return 2*(UCASE(s1[i]) < UCASE(s2[i]))-1;
}


struct LINE_READER {
	FileReader *file;
	char *buffer;
	int bufferSize;
	int lineStart;
	/* Note:  -1 indicates still reading data, so buffer is full up to the end.
	 * Otherwise, has last valid value.
	 */
	int bufferEnd;
	int ignoreEmptyLines;
};

int ReadInt(char *s, int len, int *out) {
	int success;
	char c = s[len];
	s[len] = 0;
	success = (1 == sscanf(s, "%i", out));
	if (!success) *out = 0;
	s[len] = c;
	return success;
}


/* Curiously, seems to flicker less (At least on Windows) if done 3 characters at a time,
 * rather than with a single buffer. */
void CleanupAndBackspace(int len) {
	while (len) {
		fprintf(stderr, "\b \b");
		len--;
	}
}

void DisplayProgress(int *backspace, char *progressTemplate, int count) {
	int oldBackspace = *backspace;
	char temp[100];
	memset(temp, '\b', oldBackspace*sizeof(char));
	fwrite(temp, sizeof(char), oldBackspace, stderr);
	*backspace = fprintf(stderr, progressTemplate, count);
	if (*backspace < oldBackspace) {
		int diff = (oldBackspace-*backspace)*sizeof(char);
		memset(temp, ' ', diff);
		memset(temp+diff, '\b', diff);
		fwrite(temp, 2, diff, stderr);
	}
}
