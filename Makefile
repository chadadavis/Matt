CC = gcc
CFLAGS = -lm -O3

Matt : src/chain.c src/chain.h src/pdb.c src/pdb.h src/Protein.c src/Protein.h src/secondary.c src/secondary.h src/FileReader.c src/FileReader.h src/util.c src/util.h src/Vector.c src/Vector.h src/MultipleAlignment.c src/MultipleAlignment.h src/Score.c src/Score.h src/OctTree.c src/OctTree.h src/RMSD.c src/RMSD.h src/AssemblyOrder.c src/AssemblyOrder.h src/Extend.c src/Extend.h src/MultipleAlignmentOutput.h src/MultipleAlignmentOutput.c
	$(CC) $(CFLAGS) -o bin/Matt src/Matt.c src/chain.c src/pdb.c src/Protein.c src/secondary.c src/FileReader.c src/util.c src/Vector.c src/MultipleAlignment.c src/Score.c src/OctTree.c src/RMSD.c src/AssemblyOrder.c src/Extend.c src/MultipleAlignmentOutput.c

clean :
	rm -f bin/Matt src/*~ bin/core bin/core.*
