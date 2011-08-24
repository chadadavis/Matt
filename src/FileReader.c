/* Protein.c -- Functions to read both plain and compressed (gzip, deflate) files.
 * Could have used zlib, but this is more fun.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FileReader.h"

FileReader *readerOpen(char *name);
void readerClose(FileReader *FileReader);


/* Large enough to handle reading next token in a .Z file with space left over.  Unlikey
 * to encounter it, as a file would have to contain 2 GB of a single character.
 */
#define BUFFER_SIZE (1<<17)

struct BUFFER {
	int used;
	int read;
	int minSize;
	unsigned char data[BUFFER_SIZE];
};
typedef struct BUFFER Buffer;

/* deflated file - read only */
struct DFILE;
typedef struct DFILE DFile;

/* gzipped file - read only */
struct GZFILE;
typedef struct GZFILE GZFile;

struct FILE_READER {
	FILE *file;
	DFile *dfile;
	GZFile *gzfile;
	Buffer buffer;
};


struct EXTENEDED_CODE {
	unsigned int lengths;
	unsigned int bits;
};

typedef struct EXTENEDED_CODE ExtendedCode;

const ExtendedCode lengthTable[29] = {{  3,0}, {  4,0}, {  5,0}, {  6,0}, {  7,0}, {  8,0},
									{  9,0}, { 10,0}, { 11,1}, { 13,1}, { 15,1}, { 17,1},
									{ 19,2}, { 23,2}, { 27,2}, { 31,2}, { 35,3}, { 43,3},
									{ 51,3}, { 59,3}, { 67,4}, { 83,4}, { 99,4}, {115,4},
									{131,5}, {163,5}, {195,5}, {227,5}, {258,0}};

const ExtendedCode distTable[30] =   {{   1,0}, {   2,0}, {   3,0}, {   4,0}, {   5,1}, {   7,1},
									{   9,2}, {  13,2}, {  17,3}, {  25,3}, {  33,4}, {  49,4},
									{  65,5}, {  97,5}, { 129,6}, { 193,6}, { 257,7}, { 385,7},
									{ 513,8}, { 769,8}, {1025,9}, {1537,9}, {2049,10}, {3073,10},
									{4097,11}, {6145,11}, {8193,12}, {12289,12}, {16385,13}, {24577,13}};

const unsigned int lengthOrder[19] = {16, 17, 18,
		0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};

const ExtendedCode codeLengthTable[2] = {{3,3}, {11,7}};

struct HUFFMAN_CODES {
	/* Longer than needed, but best to be on the safe side */
	unsigned int code[1200];
};

typedef struct HUFFMAN_CODES HuffmanCodes;


struct BIT_READER {
	FileReader *fileReader;
	unsigned int bits;
	unsigned int numBits;
};
typedef struct BIT_READER BitReader;

int readBits(BitReader *reader, unsigned int numBits, unsigned int *out) {
	unsigned int res = reader->bits;
	unsigned int read = reader->numBits;
	while (read < numBits) {
		unsigned char bits;
		if (readerRead(&bits, 1, 1, reader->fileReader) <= 0) {
			numBits = read;
			break;
		}
		res += (bits << read);
		read+=8;
	}
	reader->numBits = read-numBits;
	reader->bits = res >> numBits;
	*out = (res & ((1 << numBits)-1));
	return numBits;
}

void unreadBits(BitReader *reader, unsigned int numBits, unsigned int bits) {
	reader->bits = (reader->bits << numBits) + bits;
	reader->numBits += numBits;
}

void evenBits(BitReader *reader) {
	int remove = reader->bits & 7;
	reader->numBits -= remove;
	reader->bits >>= remove;
}

struct GZFILE {
	BitReader reader;
	unsigned int type, final, dead;
	HuffmanCodes charCodes, distCodes;
};

GZFile * gzopen(FileReader *fileReader) {
	unsigned char header[10];
	int crc, extra, fname, fcomment;
	GZFile *gzfile = 0;
	if (fileReader->buffer.used < 10 ||
		fileReader->buffer.data[0] != 0x1F ||
		fileReader->buffer.data[1] != 0x8B ||
		fileReader->buffer.data[2] != 0x08 ||
		fileReader->buffer.data[3] > 32 ||
		!(gzfile = (GZFile*) calloc(1, sizeof(GZFile))) ||
		1 != readerRead(header, 10, 1, fileReader)) {
			free(gzfile);
			return 0;
	}
	crc = header[3] & 2;
	extra = header[3] & 4;
	fname = header[3] & 8;
	fcomment = header[3] & 16;
	if (extra) {
		unsigned short xlen;
		if (readerRead(&header, 2, 1, fileReader) == 1) {
			xlen = header[0] + 256*(unsigned short)header[1];
			while (xlen--) {
				readerRead(header, 1, 1, fileReader);
			}
		}
	}
	if (fname) {
		while (1 == readerRead(header, 1, 1, fileReader) && header[0]);
	}
	if (fcomment) {
		while (1 == readerRead(header, 1, 1, fileReader) && header[0]);
	}
	if (crc) {
		unsigned short crc16;
		readerRead(&crc16, 2, 1, fileReader);
	}

	gzfile->reader.fileReader = fileReader;

	return gzfile;
}

unsigned int getHuffmanCode(BitReader *reader, HuffmanCodes *hc) {
	unsigned int pos, bits;
	unsigned int read, consumed = 8;
	read = readBits(reader, 16, &bits);
	pos = hc->code[bits & 0xFF];
	while (pos < 0x8000) {
		pos = hc->code[pos + ((bits>> consumed) & 0x3)];
		consumed += 2;
	}
	consumed = pos>>16;
	if (consumed <= read) {
		unreadBits(reader, read-consumed, bits >> consumed);
		return pos & 0x7FFF;
	}
	return 0xFFFF;
}


int initHuffmanCodes(HuffmanCodes *hc, const unsigned int* lengths, const unsigned int num) {
	unsigned int numCodes[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	unsigned int codes[16][300];
	unsigned int i, j, code = 0;
	unsigned int bit = 1<<15;
	unsigned int next = 256;
	unsigned int k;

	for (i=0; i<num; i++) {
		codes[lengths[i]][numCodes[lengths[i]]++] = i;
	}
	for (i=0; i<256; i++) {
		hc->code[i] = 0x8FFF;
	}
	for (i=1; i<16; i++) {
		unsigned int numQuads = ((i+1)&~1)-1;
		for (j=0; j<numCodes[i]; j++) {
			unsigned int pos = 0;
			unsigned int d = ((code>>15)&1) + ((code>>13)&2) + ((code>>11)&4) + ((code>>9)&8) +
							 ((code>>7)&16) + ((code>>5)&32) + ((code>>3)&64) + ((code>>1)&128);
			if (i<=8) {
				int inc = 1<<i;
				int code2 =  codes[i][j] + 0x8000 + (i<<16);
				for (k = d; k < 256; k+=inc) {
					hc->code[k] = code2;
				}
				code+=bit;
				continue;
			}
			else {
				if (hc->code[d] & 0x8000) {
					hc->code[d] = next;
					hc->code[next] = 0x8FFF;
					hc->code[next+1] = 0x8FFF;
					hc->code[next+2] = 0x8FFF;
					hc->code[next+3] = 0x8FFF;
					pos = next;
					next += 4;
				}
				else {
					pos = hc->code[d];
				}
			}
			for (k=9; k<numQuads; k+=2) {
				d = ((code>>(16-k))&1) + ((code>>(14-k))&2);
				if (hc->code[d+pos] & 0x8000) {
					hc->code[d+pos] = next;
					hc->code[next] = 0x8FFF;
					hc->code[next+1] = 0x8FFF;
					hc->code[next+2] = 0x8FFF;
					hc->code[next+3] = 0x8FFF;
					pos = next;
					next += 4;
				}
				else {
					pos = hc->code[d+pos];
				}
			}

			d = ((code>>(16-k))&1);
			if (i&1) {
				hc->code[d+pos] = codes[i][j] + 0x8000 + (i<<16);
				hc->code[d+pos+2] = codes[i][j] + 0x8000 + (i<<16);
			}
			else {
				d += ((code>>(14-k))&2);
				hc->code[d+pos] = codes[i][j] + 0x8000 + (i<<16);
			}

			code+=bit;
			/* Reached only when table is overfull.  Overflow. */
			if (code==0) {
				return 0;
			}
		}
		bit>>=1;
	}
	return 1;
}


int gzread(Buffer *buffer, GZFile *gzfile) {
	int read = 0;
	unsigned int i;
	while (buffer->used + (1<<15) < BUFFER_SIZE) {
		if (gzfile->dead) break;
		if (!gzfile->type) {
			unsigned int codeLengths[350];
			unsigned int hlit, hdist, hclen, test;
			/* Needed for repeat codes. */
			buffer->minSize = 32768;
			if (readBits(&gzfile->reader, 1, &gzfile->final) != 1 || readBits(&gzfile->reader, 2, &gzfile->type) != 2 || gzfile->type == 3) {
				gzfile->dead = 1;
				return read;
			}
			if (gzfile->type == 0) {
				int len, nlen;
				int len2;
				unsigned int bytes[4];
				evenBits(&gzfile->reader);
				if (readBits(&gzfile->reader, 16, &bytes[0]) != 8 ||
					readBits(&gzfile->reader, 16, &bytes[1]) != 8 ||
					readBits(&gzfile->reader, 16, &bytes[2]) != 8 ||
					readBits(&gzfile->reader, 16, &bytes[3]) != 8) {
					gzfile->dead = 1;
					return read;
				}
				len = bytes[0] + (((int)bytes[1]) << 8);
				nlen = bytes[2] + (((int)bytes[3]) << 8);
				if (len ^ nlen ^ 0xFFFF) {
					gzfile->dead = 1;
				}
				/* Always byte-aligned after calling even bits and reading 16 bits, so don't need to use readBits here.*/
				len2 = readerRead(buffer->data+buffer->used, 1, len, gzfile->reader.fileReader);
				if (len2 != len || gzfile->final) {
					gzfile->dead = 1;
					if (len2 < 0)
						return read;
				}
				read += len2;
				continue;
			}
			else if (gzfile->type == 2) {
				HuffmanCodes lengths;
				if (readBits(&gzfile->reader, 5, &hlit) != 5 ||
					(hlit += 257) > 286 ||
					readBits(&gzfile->reader, 5, &hdist) != 5 ||
					(hdist += 1) > 30 ||
					readBits(&gzfile->reader, 4, &hclen) != 4 ||
					(hclen += 4) > 19) {
						gzfile->dead = 1;
						return read;
				}
				for (i=0; i<hclen; i++) {
					if (readBits(&gzfile->reader, 3, &codeLengths[lengthOrder[i]]) != 3) {
						gzfile->dead = 1;
						return read;
					}
				}
				for (;i<19;i++) {
					codeLengths[lengthOrder[i]] = 0;
				}

				if (!initHuffmanCodes(&lengths, codeLengths, 19)) {
					gzfile->dead = 1;
					return read;
				}
				for (i=0; i<hlit+hdist; i++) {
					codeLengths[i] = getHuffmanCode(&gzfile->reader, &lengths);
					if (codeLengths[i]>15) {
						unsigned int len;
						if (codeLengths[i]==16) {
							if (readBits(&gzfile->reader, 2, &len) != 2) {
								gzfile->dead = 1;
								return read;
							}
							/* Actually 3 + len, but first is current one */
							len += 2;
							codeLengths[i] = codeLengths[i-1];
							while (len--) {
								i++;
								codeLengths[i] = codeLengths[i-1];
							}
						}
						else if (codeLengths[i] < 20) {
							if (readBits(&gzfile->reader, codeLengthTable[codeLengths[i]-17].bits, &len) != (int)codeLengthTable[codeLengths[i]-17].bits) {
								gzfile->dead = 1;
								return read;
							}
							/* -1 is for buttent entry. */
							len += -1 + codeLengthTable[codeLengths[i]-17].lengths;
							codeLengths[i]=0;
							while (len--) {
								i++;
								codeLengths[i] = 0;
							}
						}
						else {
							gzfile->dead = 1;
							return read;
						}
					}
				}
				test = 0;
				for (i=0; i<hlit; i++) {
					if (codeLengths[i])
						test += (1<<16) >> codeLengths[i];
				}
			}
			else {
				hlit = 288;
				hdist = 30;
				for (i=0; i<144; i++)
					codeLengths[i] = 8;
				for (   ; i<256; i++)
					codeLengths[i] = 9;
				for (   ; i<280; i++)
					codeLengths[i] = 7;
				for (   ; i<288; i++)
					codeLengths[i] = 8;

				for (i=0; i<30; i++)
					codeLengths[hlit+i] = 5;
			}
			if (!initHuffmanCodes(&gzfile->charCodes, codeLengths, hlit) || !initHuffmanCodes(&gzfile->distCodes, codeLengths+hlit, hdist)) {
				gzfile->dead = 1;
				return read;
			}
		}
		while (buffer->used + (1<<15) < BUFFER_SIZE) {
			unsigned int c = getHuffmanCode(&gzfile->reader, &gzfile->charCodes);
			if (c < 256) {
				read++;
				buffer->data[buffer->used++] = (unsigned char) c;
			}
			else if (c>256) {
				c-=257;
				if (c < 29) {
					unsigned int len, dist, distEntry;
					if (readBits(&gzfile->reader, lengthTable[c].bits, &len) != (int)lengthTable[c].bits) {
						gzfile->dead = 1;
						return read;
					}
					len += lengthTable[c].lengths;
					distEntry = getHuffmanCode(&gzfile->reader, &gzfile->distCodes);
					if (distEntry >= 30 || readBits(&gzfile->reader, distTable[distEntry].bits, &dist) != (int)distTable[distEntry].bits) {
						gzfile->dead = 1;
						return read;
					}
					dist += distTable[distEntry].lengths;
					if ((int)dist >= buffer->used) {
						gzfile->dead = 1;
						return read;
					}
					for (i=0; i<len; i++) {
						buffer->data[buffer->used] = buffer->data[buffer->used-dist];
						buffer->used++;
					}
					read += len;
				}
				else {
					/* Error */
					gzfile->dead = 1;
					return read;
				}
			}
			else {
				/* End of block. */
				gzfile->dead = gzfile->final;
				if (gzfile->final) {
					gzfile->final = gzfile->final;
				}
				gzfile->type = 0;
				break;
			}
		}
	}
	return read;
}

void gzclose(GZFile *gzfile) {
	if (gzfile) {
		readerClose(gzfile->reader.fileReader);
		free(gzfile);
	}
}



struct compressTable {
	unsigned int size;
	unsigned short prefix[0x10000];
	char suffix[0x10000];
};
typedef struct compressTable CompressTable;

struct DFILE {
	int blockMode;
	int maxBits;
	int bits;
	int prev;
	int end;
	int count;
	BitReader reader;
	CompressTable table;
};

DFile *dopen(FileReader *fileReader) {
	DFile *dfile;
	int maxBits;
	unsigned char header[3];
	if (!fileReader || fileReader->buffer.used <= 2 ||
		fileReader->buffer.data[0] != 0x1F || fileReader->buffer.data[1] != 0x9D ||
		fileReader->buffer.data[2] & 0x60) return 0;
	maxBits = fileReader->buffer.data[2] & 0x1F;
	if (maxBits < 8 || maxBits > 16) return 0;

	readerRead(header, 1, 3, fileReader);

	dfile = (DFile*) calloc(sizeof(DFile), 1);

	dfile->reader.fileReader = fileReader;
	dfile->blockMode = header[2] & 0x80;
	dfile->maxBits = maxBits;
	dfile->bits = 9;
	dfile->table.size = 256;
	if (dfile->blockMode) {
		dfile->table.size = 257;
	}
	return dfile;
}
void dclose(DFile *dfile) {
	if (dfile) {
		readerClose(dfile->reader.fileReader);
		free(dfile);
	}
}

int dread(Buffer *buffer, DFile *dfile) {
	unsigned int code;
	int read = 0;
	while (buffer->used < BUFFER_SIZE-(1<<16)) {
		if (readBits(&dfile->reader, dfile->bits, &code) != dfile->bits) {
			dfile->end = 1;
			break;
		}
		dfile->count ++;
		if (code == 256 && dfile->blockMode) {
			dfile->count &= 7;
			if (dfile->count) {
				while (dfile->count < 8) {
					readBits(&dfile->reader, dfile->bits, &code);
					dfile->count++;
				}
				dfile->count = 0;
			}
			dfile->table.size = 257;
			dfile->bits = 9;
		}
		else {
			int size = 1;
			int i = code;
			if (code >= dfile->table.size) {
				dfile->end = 1;
				break;
			}
			while (i > 255) {
				size ++;
				i = dfile->table.prefix[i];
			}
			i = code;
			buffer->used += size;
			size = 1;
			while (i > 255) {
				buffer->data[buffer->used - size] = dfile->table.suffix[i];
				i = dfile->table.prefix[i];
				size++;
			}
			buffer->data[buffer->used - size] = (char)i;
			read += size;
			if (dfile->table.size-1 < (unsigned int)(1<<dfile->maxBits)) {
				dfile->table.suffix[dfile->table.size-1] = i; /* first bit output. */
				if (dfile->table.size < (unsigned int)(1<<dfile->maxBits)) {
					dfile->table.prefix[dfile->table.size] = code;
					dfile->table.suffix[dfile->table.size] = i; /* first bit output. */
				}
				dfile->table.size++;
				if (((dfile->table.size-1) & (1<<dfile->bits)) && dfile->bits != dfile->maxBits) {
					dfile->bits++;
				}
			}
			dfile->prev = code;
		}
	}
	return read;
}



FileReader *readerOpenSub(FILE *file, int format) {
	FileReader *fileReader;
	if (!file) return 0;
	fileReader = (FileReader*) malloc(sizeof(FileReader));
	if (!fileReader) {
		fclose(file);
		return 0;
	}
	fileReader->dfile = 0;
	fileReader->gzfile = 0;
	fileReader->file = file;
	fileReader->buffer.used = 0;
	fileReader->buffer.read = 0;
	fileReader->buffer.minSize = 0;

	while (1) {
		DFile *dfile;
		GZFile *gzfile;
		char temp[10];
		if (readerRead(temp, sizeof(temp), 1, fileReader) != 1) break;
		fileReader->buffer.read-=sizeof(temp);
		if ((!format || format == 1) &&
			(dfile = dopen(fileReader))) {
					FileReader *fileReader2 = (FileReader*) malloc(sizeof(FileReader));
					format = 0;
					if (!fileReader2) {
						dclose(dfile);
						return 0;
					}
					fileReader2->dfile = dfile;
					fileReader2->gzfile = 0;
					fileReader2->buffer.used = 0;
					fileReader2->buffer.read = 0;
					fileReader2->file = 0;
					fileReader2->buffer.minSize = 0;
					fileReader = fileReader2;
		}
		else if ((!format || format == 2) &&
				 (gzfile = gzopen(fileReader))) {
					FileReader *fileReader2 = (FileReader*) malloc(sizeof(FileReader));
					format = 0;
					if (!fileReader2) {
						dclose(dfile);
						return 0;
					}
					fileReader2->dfile = 0;
					fileReader2->gzfile = gzfile;
					fileReader2->buffer.used = 0;
					fileReader2->buffer.read = 0;
					fileReader2->file = 0;
					fileReader2->buffer.minSize = 0;
					fileReader = fileReader2;
		}
		else break;
	}
	if (format) {
		readerClose(fileReader);
		return 0;
	}
	return fileReader;
}


FileReader *readerOpen(char *name) {
	FileReader *fileReader = readerOpenSub(fopen(name, "rb"), 0);
	if (!fileReader) {
		char *temp = (char*) malloc(sizeof(char)*(strlen(name)+10));
		sprintf(temp, "%s.Z", name);
		if (!(fileReader = readerOpenSub(fopen(temp, "rb"), 1))) {
			sprintf(temp, "%s.gz", name);
			fileReader = readerOpenSub(fopen(temp, "rb"), 2);
		}
		free(temp);
	}

	return fileReader;
}

void readerClose(FileReader *fileReader) {
	if (fileReader) {
		if (fileReader->file) fclose(fileReader->file);
		else if (fileReader->dfile) dclose(fileReader->dfile);
		else if (fileReader->gzfile) gzclose(fileReader->gzfile);
		free(fileReader);
	}
}

int readerRead(void *out, int size, int count, FileReader *fileReader) {
	int read = 0;
	while (read<count) {
		if (fileReader->buffer.used - fileReader->buffer.read >= size) {
			int copyCount = count;
			int copyBytes = count*size;
			if (copyBytes > fileReader->buffer.used - fileReader->buffer.read) {
				copyCount = (fileReader->buffer.used - fileReader->buffer.read)/size;
				copyBytes = size * copyCount;
			}
			memcpy(out, fileReader->buffer.data + fileReader->buffer.read, copyBytes);
			out = ((char*)out)+copyBytes;
			fileReader->buffer.read += copyBytes;
			read += copyCount;
			count -= copyCount;
		}
		else {
			if (fileReader->buffer.minSize <= fileReader->buffer.used - fileReader->buffer.read) {
				memmove(fileReader->buffer.data, fileReader->buffer.data + fileReader->buffer.read, fileReader->buffer.used -= fileReader->buffer.read);
				fileReader->buffer.read = 0;
			}
			else if (fileReader->buffer.minSize < fileReader->buffer.used) {
				int keep = fileReader->buffer.minSize;
				if (keep < fileReader->buffer.used - fileReader->buffer.used)
					keep = fileReader->buffer.used - fileReader->buffer.used;
				memmove(fileReader->buffer.data, fileReader->buffer.data + fileReader->buffer.used - keep, keep);
				fileReader->buffer.read -= fileReader->buffer.used - keep;
				fileReader->buffer.used = keep;
			}

			if (fileReader->file) {
				/* Buffering for fread doesn't really make much difference */
				int read = (int)fread(fileReader->buffer.data + fileReader->buffer.used, 1, BUFFER_SIZE-fileReader->buffer.used, fileReader->file);
				if (read <= 0) break;
				fileReader->buffer.used += read;
			}
			else if (fileReader->dfile) {
				if (dread(&fileReader->buffer, fileReader->dfile) == 0) break;
			}
			else {
				if (gzread(&fileReader->buffer, fileReader->gzfile) == 0) break;
			}
		}
	}
	if (read || size == 0) return read;
	return -1;
}


int readerGetLine(FileReader* in, char **line, int *lineLen) {
	int i;
	if (!*line) {
		*lineLen = 100;
		*line = (char*) malloc(sizeof(char)* (2 + *lineLen));
	}
	i=0;
	while (readerRead(&line[0][i], 1, 1, in) == 1) {
		if (i>=*lineLen) {
			*lineLen *= 2;
			*line = (char*)realloc(*line, sizeof(char)* (2+*lineLen));
		}
		if (line[0][i]=='\n' || line[0][i]=='\r') {
			if (i) break;
		}
		else i++;
	}
	line[0][i]=0;
	return i;
}
