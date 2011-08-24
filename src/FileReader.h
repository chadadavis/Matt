/* Protein.h -- Prototypes to read both plain and compressed (gzip, deflate) files
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

#ifndef FILE_READER_H
#define FILE_READER_H

/* File reader. Handles raw text, gzip, and compress files.
 * Can also handle files compressed multiple times with one or
 * both algorithms. Recognizes compressed files by their header.
 */
struct FILE_READER;
typedef struct FILE_READER FileReader;

FileReader *readerOpen(char *name);
void readerClose(FileReader *FileReader);

/* Ignores empty lines. */
int readerGetLine(FileReader* in, char **line, int *lineLen);

/* Note:  Size must be <= 65536. Works much like fread(). */
int readerRead(void *out, int size, int count, FileReader *fileReader);

#endif
