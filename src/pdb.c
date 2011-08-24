/* pdb.c -- Contains functions for loading and manipulating pdb files.
 * Rather slow, but it handles small errors in source files reasonably
 * well.
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

#include "pdb.h"
#include "string.h"

void CleanupModel(PDBModel *model) {
	free(model->pdbName);
	if (model->chains) {
		int i;
		for (i=0; i<model->numChains; i++) {
			CleanupChain(model->chains[i]);
		}
		free(model->chains);
	}
	free(model);
}

PDBModel *CreateModel(char *name, int number) {
	PDBModel *model = (PDBModel *)calloc(1, sizeof(PDBModel));
	if (name)
		model->pdbName = strdup(name);
	model->number = number;
	return model;
}

void RotateModel(PDBModel *model, const Matrix *m) {
	int i;
	for (i=0; i<model->numChains; i++) {
		RotateChain(model->chains[i], m);
	}
}

PDBChain *PopChain(PDBModel *model, char chain) {
	int i;
	if (chain == 0) chain = ' ';
	for (i=0; i<model->numChains; i++) {
		if (model->chains[i]->chainName == chain) {
			PDBChain * temp = model->chains[i];
			memcpy(model->chains+i, model->chains+i+1, sizeof(PDBChain*) * (model->numChains-i-1));
			model->numChains --;
			return temp;
		}
	}
	return 0;
}

PDBChain *GetChain(PDBModel *model, char chain) {
	int i;
	if (chain == 0) chain = ' ';
	for (i=0; i<model->numChains; i++) {
		if (model->chains[i]->chainName == chain)
			return model->chains[i];
	}
	return 0;
}

PDBChain *AddNewChain(PDBModel *model, char chain) {
	PDBChain *out = GetChain(model, chain);
	if (out) return out;
	if (chain == 0) chain = ' ';
	model->chains = (PDBChain**) realloc(model->chains, sizeof(PDBChain*) * (model->numChains+1));
	model->chains[model->numChains] = CreateChain(model->pdbName, chain);
	model->numChains++;

	return model->chains[model->numChains-1];
}

void AddExistingChain(PDBModel *model, PDBChain *chain) {
	model->chains = (PDBChain**) realloc(model->chains, sizeof(PDBChain*) * (model->numChains+1));
	model->chains[model->numChains++] = chain;
}

void AddChain(PDBData *pdb, int model, PDBChain *chain) {
	int i;
	for (i=0; i<pdb->numModels; i++) {
		if (model == pdb->models[i]->number) break;
	}
	if (i == pdb->numModels) {
		AddModel(pdb, model);
	}
	AddExistingChain(pdb->models[i], chain);
}

PDBModel *AddModel(PDBData *pdb, int number) {
	pdb->models = (PDBModel **)realloc(pdb->models, sizeof(PDBModel*) * (pdb->numModels+1));
	pdb->models[pdb->numModels] = CreateModel(pdb->name, number);
	return pdb->models[pdb->numModels++];
}

PDBModel *FindModel(PDBData *pdb, int number) {
	int i;
	for (i=0; i<pdb->numModels; i++) {
		if (pdb->models[i]->number == number) return pdb->models[i];
	}
	return 0;
}


void RotatePDB(PDBData *pdb, const Matrix *m) {
	int i;
	for (i=0; i<pdb->numModels; i++) {
		RotateModel(pdb->models[i], m);
	}
}

void CleanupPDB(PDBData *pdb) {
	int i;
	if (!pdb) return;
	for (i=0; i<pdb->numModels; i++) {
		CleanupModel(pdb->models[i]);
	}
	free(pdb->models);
	free(pdb->name);
	free(pdb->header);
	free(pdb->cryst);
	free(pdb);
}

void CalcPDBSecondary(PDBData *pdb) {
	int m, c;
	for (m=0; m<pdb->numModels; m++) {
		for (c=0; c<pdb->models[m]->numChains; c++) {
			CalcChainSecondary(pdb->models[m]->chains[c]);
		}
	}
}

void OutputChains(PDBChain **chain, int numChains, FILE *out, int seqres, int het, int secondary, int extraInfo) {
	PDBData *pdb = (PDBData *)calloc(sizeof(PDBData), 1);
	pdb->numModels = 1;
	pdb->models = (PDBModel **)calloc(sizeof(PDBModel*), 1);
	pdb->models[0] = (PDBModel *)calloc(sizeof(PDBModel), 1);
	pdb->models[0]->number = 0;
	pdb->models[0]->numChains = numChains;
	pdb->models[0]->chains = (PDBChain**)memdup(chain, sizeof(PDBChain*)*numChains);
	Output(pdb, out, seqres, het, secondary, extraInfo);
	pdb->models[0]->numChains = 0;
	CleanupPDB(pdb);
}

void Output(PDBData *pdb, FILE *out, int seqres, int het, int secondary, int extraInfo) {
	int atomIndex = 1;
	int c, m, i;
	int numRemarks = 0;
	int numHet = 0;
	int sheetLines = 0;
	int numTurns = 0;
	int numSites = 0;
	int numXform = 0;
	int numHelices = 0;
	int numSheets = 0;
	int numConnect = 0;
	int numSeqres = 0;
	int numTer = 0;
	if (extraInfo) {
		if (pdb->header) {
			char *temp = pdb->header;
			do {
				int output = 0;
				if (!stristart("REMARK 987 CHAIN", temp) && !stristart("FROM", temp+19)) {
					output = 1;
					numRemarks++;
				}
				else if (!stristart("REMARK   1", temp) ||
						 !stristart("REMARK   2", temp) ||
						 !stristart("REMARK   3", temp)) {
							 if (extraInfo >= 2) {
								output = 1;
								numRemarks++;
							}
				}
				else if (!stristart("REMARK", temp)) {
					 if (extraInfo >= 3) {
						output = 1;
						numRemarks++;
					}
				}
				else if (extraInfo >= 2) output = 1;
				if (output) {
					fwrite(temp, 1, 80, out);
					fprintf(out, "\n");
				}
				temp += 81;
			}
			while (temp[-1]);
		}
		if (extraInfo >= 2 && pdb->models) {
			for (c=0; c<pdb->models[0]->numChains; c++) {
				PDBChain *chain = pdb->models[0]->chains[c];
				for (i=0; i<chain->numDbrefs; i++) {
					DBRef * dbref = chain->dbrefs+i;
					fprintf(out, "DBREF  %s %c %4i%c %4i%c %-54s\n", dbref->idCode, chain->chainName, chain->residues[dbref->start].resNum, chain->residues[dbref->start].resLetter, chain->residues[dbref->end].resNum, chain->residues[dbref->end].resLetter, dbref->stuff);
				}
			}
			for (c=0; c<pdb->models[0]->numChains; c++) {
				PDBChain *chain = pdb->models[0]->chains[c];
				for (i=0; i<chain->numSeqAdvs; i++) {
					SeqAdv * seqAdv = chain->seqAdvs+i;
					fprintf(out, "SEQADV %s %s %c %4i%c %-56s\n", seqAdv->idCode, seqAdv->resName, chain->chainName, chain->residues[seqAdv->res].resNum, chain->residues[seqAdv->res].resLetter, seqAdv->stuff);
				}
			}
		}
	}
	if (seqres && pdb->models) {
		for (c=0; c<pdb->models[0]->numChains; c++) {
			numSeqres += (pdb->models[0]->chains[c]->length+12)/13;
			OutputChainSeqres(pdb->models[0]->chains[c], out, het);
		}
	}
	if (extraInfo >= 2 && pdb->cryst) {
		char *c = pdb->cryst-1;
		while (c) {
			c++;
			if (!stristart("ORIGX", c) ||
				!stristart("SCALE", c) ||
				!stristart("MTRIX", c))
				numXform++;
			c = strchr(c, '\n');
		}
		fprintf(out, "%s\n", pdb->cryst);
	}
	if (secondary && pdb->models) {
		char *helixIDs=0;
		char *strandIDs=0;
		int numCisPeps = 0;
		int numSSbonds = 0;
		for (c=0; c<pdb->models[0]->numChains; c++) {
			OutputChainAlphas(pdb->models[0]->chains[c], out, het, &numHelices, &helixIDs);
		}
		free(helixIDs);
		for (c=0; c<pdb->models[0]->numChains; c++) {
			OutputChainSheets(pdb->models[0]->chains[c], out, het, &numSheets, &strandIDs, &sheetLines);
		}
		free(strandIDs);
		for (c=0; c<pdb->models[0]->numChains; c++) {
			OutputChainSSbonds(pdb->models[0]->chains[c], out, het, &numSSbonds);
		}
		for (c=0; c<pdb->models[0]->numChains; c++) {
		}
		for (m=0; m<pdb->numModels; m++) {
			for (c=0; c<pdb->models[m]->numChains; c++) {
				for (i=0; i<pdb->models[m]->chains[c]->numCisPeps; i++) {
					OutputChainCisPeps(pdb->models[m]->chains[c], out, het, &numCisPeps, pdb->models[m]->number);
				}
			}
		}
	}
	for (m=0; m<pdb->numModels; m++) {
		for (c=0; c<pdb->models[m]->numChains; c++) {
			if (!pdb->models[m]->chains[c]->length) continue;
			numTer++;
			OutputChain(pdb->models[m]->chains[c], out, seqres, het, &atomIndex);
		}
	}
	fprintf(out, "MASTER    %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i          \n",
		numRemarks%100000,         0%100000,        numHet%100000,     numHelices%100000,
		sheetLines%100000,         numTurns%100000, numSites%100000,   numXform%100000,
		(atomIndex-numTer)%100000, numTer%100000,   numConnect%100000, numSeqres%100000);
	fprintf(out, "%-80s\n", "END");
}

void OutputFasta(PDBData *pdb, FILE *out, int het) {
	int c;
	for (c=0; c<pdb->models[0]->numChains; c++) {
		OutputChainFasta(pdb->models[0]->chains[c], out, het);
	}
}





struct RAW_RESIDUE_POSITION {
	int index;
	char chain;
	char resLetter;
	char resNumber;
};
typedef struct RAW_RESIDUE_POSITION RawResiduePosition;

struct TEMP_ALPHA_HELIX {
	char id[4];
	RawResiduePosition start, end;
};
typedef struct TEMP_ALPHA_HELIX TempAlphaHelix;

char HeaderOrder[][11] = {
	"HEADER",
	"OBSLTE",
	"TITLE ",
	"CAVEAT",
	"COMPND",
	"SOURCE",
	"KEYWDS",
	"EXPDTA",
	"AUTHOR",
	"REVDAT",
	"SPRSDE",
	"JRNL  ",
	"REMARK   1",
	"REMARK   2",
	"REMARK   3",
	"REMARK",
};

char CrystOrder[][11] = {
	"CRYST1",
	"ORIGX1",
	"ORIGX2",
	"ORIGX3",
	"SCALE1",
	"SCALE2",
	"SCALE3",
	"MTRIX1",
	"MTRIX2",
	"MTRIX3",
	"TVECT "};

#define CACHED_DBREF  0
#define CACHED_SEQADV 1
#define CACHED_HELIX  2
#define CACHED_SHEET  3
#define CACHED_TURN   4
#define CACHED_SSBOND 5
#define CACHED_LINK   6
#define CACHED_CISPEP 7
#define CACHED_TYPES  8

char CachedStrings[CACHED_TYPES][8] = {
	"DBREF ",
	"SEQADV",
	"HELIX ",
	"SHEET ",
	"TURN  ",
	"SSBOND",
	"LINK  ",
	"CISPEP",
};




PDBData * LoadFileReader(FileReader *FileReader, char *suggestedName, int loadSeqres, int maxGaps, int caonly, int het) {
	PDBData *out = 0;
	char *cachedStrings[CACHED_TYPES];
	int seqres = 0;

	int modelNum = 0;

	PDBModel *model = 0;
	PDBChain *chain = 0;

	int started=0;
	int i, j;
	char *line = 0;
	int lineLen = 0;
	int modelStart = 0;
	char *string;

	int HET = 0;

	memset(cachedStrings, 0, sizeof(cachedStrings));
	if (!FileReader) return out;
	out = (PDBData*)calloc(1, sizeof(PDBData));
	if (suggestedName)
		out->name = strdup(suggestedName);
	else
		out->name = strdup("none");

	while ((i=readerGetLine(FileReader, &line, &lineLen)) > 0) {
		if (loadSeqres && stristart("SEQRES", line)==0 && i>=21) {
			seqres++;
			if (!model)
				model = AddModel(out, 0);
			line[11] = UCASE(line[11]);
			if (!chain || chain->chainName != line[11]) {
				chain = AddNewChain(model, line[11]);
			}
			/* SEQRES code which means no SEQRES data.  Useful!
			 * Currently add chain before skipping rest of the handling
			 * UNK skipped because it gives no useful information here, either.
			 */
			if (stristart(" 0", &line[8]) == 0 &&
				stristart("UNK", &line[19]) == 0) {
				continue;
			}
			if (chain->terminated) continue;
			if (i>67) i = 67;
			for (j=19; j<=i; j+=4) {
				char resName[4] = {UCASE(line[j]), UCASE(line[j+1]), UCASE(line[j+2]), 0};
				if (!strcmp(resName, "   ") || resName[0] == 0) break;
				AddResidue(chain, resName);
			}
		}
		else if (i>53 && (stristart("ATOM  ", line)==0 || (het && stristart("HETATM", line)==0 && (HET=1)))) {
			Atom *atom;
			int atomIndex;
			int foundCharge = 0;
			int altLoc;
			int temp;
			modelStart = 1;
			if (!model)
				model = AddModel(out, 0);
			if (!chain || chain->chainName != line[21])
				chain = AddNewChain(model, line[21]);
			if (chain->terminated)
				continue;
			chain->atoms = (Atom*)realloc(chain->atoms, (++(chain->numAtoms)) * sizeof(Atom));
			chain->tempAtoms |= !HET;

			atomIndex = chain->numAtoms-1;
			atom = &chain->atoms[atomIndex];
			memset(atom, 0, sizeof(Atom));

			strncpy(atom->name, &line[12],4);

			if (!ReadInt(&line[22], 4, &atom->res)) continue;

			ReadInt(&line[6], 5, &atom->atomIndex);
			atom->altLoc = line[16];
			if (atom->altLoc != ' ') {
				atom->altLoc = atom->altLoc;
			}

			if (!ReadDouble(&line[30], 8, &atom->pos.x) ||
				!ReadDouble(&line[38], 8, &atom->pos.y) ||
				!ReadDouble(&line[46], 8, &atom->pos.z)) continue;

			if (i>=66) {
				ReadFloat(&line[54], 6, &atom->occupancy);
				ReadFloat(&line[60], 6, &atom->tempFactor);
			}

			atom->resLetter = line[26];

			strncpy(atom->resName, &line[17],3);
			if (caonly && !strcmp(atom->resName, "CA")) {
				chain->numAtoms--;
				continue;
			}

			atom->charge = 0;
			if (i >= 80 && line[79] >= '0' && line[79] <= '9') {
				if (line[80] == '+') {
					atom->charge = line[79]-'0';
					foundCharge = 1;
				}
				else if (line[80] == '-') {
					atom->charge = -(line[79]-'0');
					foundCharge = 1;
				}
			}
			if (i >= 80 && (foundCharge || 
				(line[78]== ' ' && line[79]== ' ' && line[77] >= 'A' && line[77] <= 'Z' && (line[76]== ' ' || (line[76] >= 'A' && line[76] <= 'Z'))))) {
					atom->element[0] = line[76];
					atom->element[1] = line[77];
			}
			else {
				atom->element[0] = line[12];
				if (atom->element[0]<'A' || atom->element[0]>'Z') atom->element[0] = ' ';
				atom->element[1] = line[13];
			}

			altLoc = line[16];

			temp = chain->numTempResidues;
			while (temp && (atom->res != chain->tempResidues[temp-1].res ||
							atom->resLetter != chain->tempResidues[temp-1].resLetter)) {
				/* If same number appears in non-adjacent lines, assume they are different residues.
				 * Not always true, but need for HOMSTRAD.
				 */
				if (atom->res != chain->tempResidues[temp-1].res) {
					temp = 0;
					break;
				}
				temp--;
			}
			if (temp && strcmp(atom->resName, chain->tempResidues[temp-1].resName) &&
				temp!=chain->numTempResidues) temp = 0;

			if (temp == 0) {
				temp = chain->numTempResidues;
				chain->numTempResidues++;
				chain->tempResidues = (TempResidue*)realloc(chain->tempResidues, sizeof(TempResidue)*chain->numTempResidues);

				chain->tempResidues[temp].res = atom->res;
				chain->tempResidues[temp].resLetter = atom->resLetter;

				strcpy(chain->tempResidues[temp].resName, atom->resName);
				chain->tempResidues[temp].het = HET;
				chain->tempResidues[temp].finalIndex = -1;
				chain->tempResidues[temp].alphaCarbon = 0;
				temp = 1;
			}
			else {
				if (chain->tempResidues[temp-1].het && !HET) {
					chain->tempResidues[temp-1].het = 0;
					strcpy(chain->tempResidues[temp-1].resName, atom->resName);
				}
				while (chain->atoms[atomIndex].res != chain->atoms[atomIndex-1].res ||
					chain->atoms[atomIndex].resLetter != chain->atoms[atomIndex-1].resLetter) {
						Atom a = chain->atoms[atomIndex-1];
						chain->atoms[atomIndex-1] = chain->atoms[atomIndex];
						chain->atoms[atomIndex] = a;
						atomIndex--;
				}
				atom = &chain->atoms[atomIndex];
			}
			if (!strcmp(atom->name, " CA ")) {
				chain->tempResidues[temp-1].alphaCarbon = 1;
			}
			HET = 0;
		}
		else if (stristart("TER", line)==0 && i>=21) {
			PDBChain *lchain = chain;
			chain = GetChain(model, line[21]);
			if (chain != lchain) {
				/* Hack for 1dx5, which has a typo, and its many friends. */
				if (!chain || chain->terminated || chain->numAtoms == 0) chain = lchain;
			}
			if (!chain || chain->terminated) continue;
			Terminate(chain, maxGaps, 1);
		}
		else if (stristart("MODEL", line)==0) {
			int number;
			sscanf(&line[10], "%i", &number);
			if (!model || modelStart) {
				PDBModel *oldModel = 0;
				if (seqres)
					oldModel = out->models[out->numModels-1];
				model = AddModel(out, number);
				if (oldModel) {
					for (i=0; i<oldModel->numChains; i++) {
						PDBChain *oldChain = oldModel->chains[i];
						PDBChain *chain = AddNewChain(model, oldChain->chainName);
						for (j=0; j<oldChain->length; j++)
							AddResidue(chain, oldChain->residues[j].resName);
					}
				}
			}
			else {
				model->number = number;
			}
			modelStart = 1;
			chain = 0;
		}
		else if (stristart("ENDMDL", line)==0) {
			model = 0;
		}
		else {
			int h;
			for (h=0; h<CACHED_TYPES; h++) {
				if (stristart(CachedStrings[h], line)==0) {
					int len = 0, newLen;
					if (cachedStrings[h]) {
						len = (int)strlen(cachedStrings[h]);
						cachedStrings[h][len++] = '\n';
					}
					newLen = len + 81;
					cachedStrings[h] = (char*) realloc(cachedStrings[h], newLen);
					if (i >= 72) i = 72;
					memcpy(cachedStrings[h]+len, line, i);
					while (i<80) {
						cachedStrings[h][len+i] = ' ';
						i++;
					}
					cachedStrings[h][len+80] = 0;
					break;
				}
			}
			if (h<CACHED_TYPES) continue;
			for (h=0; h<sizeof(HeaderOrder)/sizeof(HeaderOrder[0]); h++) {
				if (stristart(HeaderOrder[h], line)==0) {
					int len = 0, newLen;
					if (out->header) {
						len = (int)strlen(out->header);
						out->header[len++] = '\n';
					}
					newLen = len + 81;
					out->header = (char*) realloc(out->header, newLen);
					if (i >= 72) i = 72;
					memcpy(out->header+len, line, i);
					while (i<80) {
						out->header[len+i] = ' ';
						i++;
					}
					out->header[len+80] = 0;
					while (len) {
						char temp[80];
						int h2;
						int pos = len-81;
						for (h2=0; h2<sizeof(HeaderOrder)/sizeof(HeaderOrder[0]); h2++) {
							if (stristart(HeaderOrder[h2], out->header+pos)==0) {
								break;
							}
						}
						if (h2 <= h) break;
						memcpy(temp, out->header+len, 80);
						memcpy(out->header+len, out->header+pos, 80);
						memcpy(out->header+pos, temp, 80);
						len -= 81;
					}
					break;
				}
			}
			if (h<sizeof(HeaderOrder)/sizeof(HeaderOrder[0])) continue;
			for (h=0; h<sizeof(CrystOrder)/sizeof(CrystOrder[0]); h++) {
				if (stristart(CrystOrder[h], line)==0) {
					int len = 0, newLen;
					if (out->cryst) {
						len = (int)strlen(out->cryst);
						out->cryst[len++] = '\n';
					}
					newLen = len + 81;
					out->cryst = (char*) realloc(out->cryst, newLen);
					if (i >= 72) i = 72;
					memcpy(out->cryst+len, line, i);
					while (i<80) {
						out->cryst[len+i] = ' ';
						i++;
					}
					out->cryst[len+80] = 0;
					while (len) {
						char temp[80];
						int h2;
						int pos = len-81;
						for (h2=0; h2<sizeof(CrystOrder)/sizeof(CrystOrder[0]); h2++) {
							if (stristart(CrystOrder[h2], out->cryst+pos)==0) {
								break;
							}
						}
						if (h2 <= h) break;
						memcpy(temp, out->cryst+len, 80);
						memcpy(out->cryst+len, out->cryst+pos, 80);
						memcpy(out->cryst+pos, temp, 80);
						len -= 81;
					}
					break;
				}
			}
		}
	}
	free(line);
	readerClose(FileReader);

	for (i=0; i<out->numModels; i++) {
		for (j=0; j<out->models[i]->numChains; j++) {
			char search[80];
			End(out->models[i]->chains[j], maxGaps);
			sprintf(search, "REMARK 987 CHAIN %c FROM ", out->models[i]->chains[j]->chainName);
			if (out->header) {
				char *s = out->header;
				while (s = strstr(s, search)) {
					if (s == out->header || s[-1] == '\n') {
						char *end = s+79;
						while (*end == ' ' && end > s) end --;
						s+=24;
						if (end > s) {
							char c = end[1];
							end[1] = 0;
							free(out->models[i]->chains[j]->idString);
							out->models[i]->chains[j]->idString = strdup(s);
							end[1] = c;
							break;
						}
					}
					s++;
				}
			}
		}
	}

	/* Remove old remark 987s. */
	if (out->header) {
		char * s = out->header;
		while (s = strstr(s, "REMARK 987 CHAIN ")) {
			char *s2 = s+80;
			while (*s2 == '\n' && !stristart("REMARK 987 CHAIN ", s2+1)) s2+=81;
			if (s != out->header) s--;
			else if (*s2) s2++;
			strcpy(s, s2);
		}
		if (!out->header[0]) {
			free(out->header);
			out->header = 0;
		}
		else {
			out->header = (char*) realloc(out->header, sizeof(char)*(1+strlen(out->header)));
		}
	}

	for (i=0; i<out->numModels; i++) {
		if (cachedStrings[CACHED_DBREF]) {
			for (string = cachedStrings[CACHED_DBREF]; string[-1]; string+=81) {
				int start, end;

				PDBChain * chain = GetChain(out->models[i], string[12]);

				if (chain && GetResidueFirstIndexString(chain, &start, string+14) && GetResidueFirstIndexString(chain, &end, string+20)) {
					chain->dbrefs = (DBRef*) realloc(chain->dbrefs, (1+chain->numDbrefs)*sizeof(DBRef));
					memcpy(chain->dbrefs[chain->numDbrefs].idCode, &string[7], 4);
					chain->dbrefs[chain->numDbrefs].idCode[4] = 0;

					chain->dbrefs[chain->numDbrefs].start = start;
					chain->dbrefs[chain->numDbrefs].end = end;
					memcpy(chain->dbrefs[chain->numDbrefs].stuff, string+26, 43);
					chain->dbrefs[chain->numDbrefs++].stuff[43] = 0;
				}
			}
		}
		if (cachedStrings[CACHED_SEQADV]) {
			for (string = cachedStrings[CACHED_SEQADV]; string[-1]; string+=81) {
				int res;

				PDBChain * chain = GetChain(out->models[i], string[16]);

				if (chain && GetResidueFirstIndexString(chain, &res, string+18)) {
					chain->seqAdvs = (SeqAdv*) realloc(chain->seqAdvs, (1+chain->numSeqAdvs)*sizeof(SeqAdv));
					memcpy(chain->seqAdvs[chain->numSeqAdvs].idCode, &string[7], 4);
					chain->seqAdvs[chain->numSeqAdvs].idCode[4] = 0;

					chain->seqAdvs[chain->numSeqAdvs].res = res;

					memcpy(chain->seqAdvs[chain->numSeqAdvs].resName, string+12, 3);
					chain->seqAdvs[chain->numSeqAdvs].resName[3] = 0;

					memcpy(chain->seqAdvs[chain->numSeqAdvs].stuff, string+24, 46);
					chain->seqAdvs[chain->numSeqAdvs++].stuff[46] = 0;
				}
			}
		}
		if (cachedStrings[CACHED_HELIX]) {
			for (string = cachedStrings[CACHED_HELIX]; string[-1]; string+=81) {
				int start, end;
				PDBChain * chain;
				if (string[19] != string[31]) continue;

				chain = GetChain(out->models[i], string[19]);
				if (chain && GetResidueFirstIndexString(chain, &start, string+21) && GetResidueFirstIndexString(chain, &end, string+33)) {
					AlphaHelix helix;
					helix.start = start;
					helix.length = end-start+1;
					memcpy(helix.id, string+11, 3);
					helix.id[3] = 0;
					memcpy(helix.comment, string+40, 30);
					helix.comment[30] = 0;
					string[40] = 0;
					helix.type = atoi(string+38);
					AddAlphaHelix(chain, &helix);
				}
			}
		}
		if (cachedStrings[CACHED_SHEET]) {
			char id[3] = "";
			int start = -1;
			int end = -1;
			int prevStart = -1;
			int prevEnd = -1;
			int prev, current;
			int dir, len;
			PDBChain * lastChain = 0;
			PDBChain * chain;
			BetaSheet *sheet = 0;
			for (string = cachedStrings[CACHED_SHEET]; string[-1]; string+=81) {
				if (string[21] != string[32] || string[49] != string[64]) {
					prevStart = -1;
					prevEnd = -1;
					lastChain = 0;
					sheet = 0;
					continue;
				}
				chain = GetChain(out->models[i], string[21]);
				if (!chain) continue;
				if (strncmp(string+11, id, 3) || lastChain != chain || prevStart == -1) {
					lastChain = chain;
					memcpy(id, string+11, 3);
					sheet = 0;
					if (!GetResidueFirstIndexString(chain, &prevStart, string+22) || !GetResidueFirstIndexString(chain, &prevEnd, string+34)) {
						prevStart = -1;
						prevEnd = -1;
					}
					continue;
				}
				if (string[49] != string[21] || string[39] != '1') {
					prevStart = -1;
					prevEnd = -1;
					lastChain = 0;
					sheet = 0;
					continue;
				}
				dir = atoi(string + 38);
				if (!GetResidueFirstIndexString(chain, &start, string+22) || !GetResidueFirstIndexString(chain, &end, string+33) ||
					!GetResidueFirstIndexString(chain, &current, string+50) || !GetResidueFirstIndexString(chain, &prev, string+65) ||
					(dir != 1 && dir != -1)) {
					lastChain = 0;
					sheet = 0;
					continue;
				}

				if (dir == 1) {
					if (string[42] == 'N') {
						if (current < start) current--;
						else prev++;
					}
					else {
						if (current > end) current++;
						else prev--;
					}
				}

				while (prev > prevStart && ((dir == 1 && current > start) || (dir == -1 && current < end))) {
					current -= dir;
					prev--;
				}
				len = 1;
				while (prev+len <= prevEnd && ((dir == -1 && current-len >= start) || (dir == 1 && current+len <= end))) {
					len ++;
				}
				if (len >= 2) {
					AddStrandPair(chain, prev, current, dir, len);
					if (!sheet) {
						int check = 0;
						chain->betaSheets = (BetaSheet*) realloc(chain->betaSheets, sizeof(BetaSheet) * (chain->numBetaSheets+1));
						sheet = chain->betaSheets + chain->numBetaSheets;
						memcpy(sheet->id, string+11, 3*sizeof(char));
						sheet->id[3] = 0;
						while (check < chain->numBetaSheets) {
							if (!strcmp(chain->betaSheets[check].id, sheet->id)) {
								sheet->id[0] = 0;
								break;
							}
							check++;
						}
						sheet->numStrands = 0;
						sheet->strands = 0;
						chain->numBetaSheets++;
					}
					sheet->strands = (BetaSheetPair*) realloc(sheet->strands, sizeof(BetaSheetPair) * (1+sheet->numStrands));
					sheet->strands[sheet->numStrands].direction = dir;
					sheet->strands[sheet->numStrands].res1 = prev;
					sheet->strands[sheet->numStrands].res2 = current;
					sheet->strands[sheet->numStrands].length = len;
					sheet->numStrands++;
				}
				else {
					sheet = 0;
				}

				prevStart = start;
				prevEnd = end;
			}
		}
		if (cachedStrings[CACHED_SSBOND]) {
			for (string = cachedStrings[CACHED_SSBOND]; string[-1]; string+=81) {
				int res1, res2;
				PDBChain * chain;
				if (string[15] != string[29]) continue;

				chain = GetChain(out->models[i], string[15]);
				if (chain && GetResidueFirstIndexString(chain, &res1, string+17) && GetResidueFirstIndexString(chain, &res2, string+31)) {
					chain->ssbonds = (SSbond*) realloc(chain->ssbonds, sizeof(SSbond) * (1+chain->numSSbonds));
					chain->ssbonds[chain->numSSbonds].res1 = res1;
					chain->ssbonds[chain->numSSbonds].res2 = res2;
					chain->ssbonds[chain->numSSbonds].sym1 = 0;
					chain->ssbonds[chain->numSSbonds].sym2 = 0;
					ReadInt(string+59, 6, &chain->ssbonds[chain->numSSbonds].sym1);
					ReadInt(string+66, 6, &chain->ssbonds[chain->numSSbonds].sym2);
					chain->numSSbonds++;
				}
			}
		}
	}
	if (cachedStrings[CACHED_CISPEP]) {
		for (string = cachedStrings[CACHED_CISPEP]; string[-1]; string+=81) {
			int res1, res2;
			int modelNum;
			if (!ReadInt(string+43, 3, &modelNum)) continue;
			if (string[15] != string[29]) continue;
			model = FindModel(out, modelNum);
			if (!model) continue;
			chain = GetChain(model, string[15]);
			if (chain && GetResidueFirstIndexString(chain, &res1, string+17) && GetResidueFirstIndexString(chain, &res2, string+31)) {

				chain->cisPeps = (CisPep*) realloc(chain->cisPeps, sizeof(CisPep) * (1+chain->numCisPeps));
				chain->cisPeps[chain->numCisPeps].res1 = res1;
				chain->cisPeps[chain->numCisPeps].res2 = res2;
				chain->cisPeps[chain->numCisPeps].angle = 0;
				ReadDouble(string+53, 6, &chain->cisPeps[chain->numCisPeps].angle);
				chain->numCisPeps++;
			}
		}
	}
	for (i=0; i<CACHED_TYPES; i++)
		free(cachedStrings[i]);

	return out;
}

/* Try all 12 possibilities. If path is the same pointer as name, then there's no given path. 
 * Don't care about relative/absolute paths here. */
static FileReader *TryOpen(char *path, char *name) {
	const char *searchPaths[4] = {
		"%s%c%s",
		"%s%c%s.pdb",
		"%s%c%s.ent",
		"%s%cpdb%s.ent"
	};
	FileReader *out;
	int pathLen = strlen(path);
	char *temp = (char*) malloc(pathLen + strlen(name) + 15);
	char *file = temp;
	int i;
	/* Ignore path, even though I still print it to the buffer,
	 * if it's the same as name. */
	if (name == path) {
		file += pathLen+1;
	}
	for (i=0; i<4; i++) {
		sprintf(temp, searchPaths[i], path, DIR_DELIMITER, name);
		out = readerOpen(file);
		if (out) break;
	}
	free(temp);
	return out;
}

PDBData *LoadPDBFile(char *file, int seqres, int maxGaps, int caonly, int het) {
	int i = 0, j=0;

	char *name;
	char *path;
	char *end;
	PDBData *out = 0;
	FileReader *reader = 0;
	path = strdup(file);
	name = strrchr(path, DIR_DELIMITER);
	if (name && name[1]) {
		name[0] = 0;
		name++;
	}
	else {
		name = path;
	}

	reader = TryOpen(path, name);

	if (!reader) {
		if (DIR_DELIMITER != path[0] &&
			(DIR_DELIMITER != '\\' || path[1] != ':')) {
				static char *searchPath = 0, *pdbPath = 0;
				static int haveEnvironment = 0;
				if (!haveEnvironment) {
					searchPath = getenv("MATT_SEARCH_PATH");
					pdbPath = getenv("MATT_PDB_PATH");
					haveEnvironment = 1;
				}
				if (searchPath) {
					char *search = strdup(searchPath);
					char *prefix = strtok(search, ";");
					while (prefix && !reader) {
						char *path2;
						if (path == name) {
							path2 = strdup(prefix);
						}
						else {
							path2 = (char*) malloc(strlen(prefix) + 2 + strlen(path));
							sprintf(path2, "%s%c%s", prefix, DIR_DELIMITER, path);
						}
						reader = TryOpen(path2, name);
						free(path2);
						prefix = strtok(0, ";");
					}
					free(search);
				}
				if (!reader && pdbPath) {
					int len = strlen(pdbPath) + 10 + strlen(path);
					char *path2 = (char*) malloc(2*len);
					char *path3 = path2+len;
					char subpath[5] = "";
					if (strlen(name) == 7 || strchr(name, '.') == name+7) {
						sprintf(subpath, "%c%c%c", DIR_DELIMITER, name[4], name[5]);
					}
					else if (strlen(name) >=3) {
						sprintf(subpath, "%c%c%c", DIR_DELIMITER, name[1], name[2]);
					}
					if (path == name) {
						strcpy(path2, pdbPath);
						sprintf(path3, "%s%s", pdbPath, subpath);
					}
					else {
						sprintf(path2, "%s%c%s", pdbPath, DIR_DELIMITER, path);
						sprintf(path3, "%s%s", pdbPath, subpath);
					}
					reader = TryOpen(path2, name);
					if (!reader && subpath[0]) reader = TryOpen(path3, name);
					free(path2);
				}
		}
	}

	end = strchr(name, 0);
	if (end - 5 > name) {
		/* Remove known extensions. */
		if (stristart(end-4, ".pdb") == 0)
			end[-4] = 0;
		else if (stristart(end-4, ".ent") == 0)
			end[-4] = 0;
	}

	if (reader) {
		out = LoadFileReader(reader, name, seqres, maxGaps, caonly, het);
	}
	free(path);
	return out;
}
