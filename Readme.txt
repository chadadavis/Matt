About:

Matt is a multiple protein structure alignment program. It uses local geometry to align segments of two sets of proteins, allowing limited bends in the backbones between the segments.  If you use Matt, please cite: M. Menke, B. Berger, L. Cowen, "Matt: Local Flexibility Aids Protein Multiple Structure Alignment", 2007, preprint.

Matt is licensed under the GNU public license version 2.0. If you would like to license Matt in an environment where the GNU public license is unacceptable (such as inclusion in a non-GPL software package) commercial Matt licensing is available through the MIT and Tufts offices of Technology Transfer. Contact betawrap@csail.mit.edu or cowen@cs.tufts.edu for more information. Contact mmenke@mit.edu for issues involving the code or binaries.


Compilation:

To compile under Linux, simply type "make". Note that the makefile will build Matt without OpenMP support.  To build it with OpenMP support, gcc 4.2 or higher is required.  Just edit the makefile and add the "-fopenmp" switch to the command.  Matt has not yet been tested with OpenMP under Linux..

Microsoft Visual Studio 6.0 and 2005 project files are included.  To compile with either one, just open up the corresponding project file and compile.  By default, the Visual Studio 2005 project is set to compile with OpenMP enabled.  The Express Edition does not support OpenMP, so will return an error message.  The option to disable it can be found under Project Properties > C/C++ > Language.


Installation:

To install, simply copy the binary to the directory you want Matt to run from and type in the command to run it.  Matt does not need to be in the active directory to run properly and uses no temp files.

Under Windows, if you use the binary distribution, you'll also need the SP1 version of Microsoft's Visual Studio 2005 runtime.  If it's not already installed, and odds are it isn't, you can get it at http://207.46.19.190/downloads/details.aspx?FamilyID=200b2fd9-ae1a-4a14-984d-389c36f85647&displaylang=en


Overview:

Matt takes a set of pdb files as input.  Individual chains can optionally be specified for individual source files.  Source files can be uncompressed or compressed in gzip or compress file formats.  The ".Z" and ".gz" extensions can optionally be left off the file name, as can a terminal ".ent" and ".pdb", and the combined prefix/suffix of "pdb" and ".ent" when they occur together.  A complete list of command line options appear at the bottom of this file.


File Formats:

Up to fifteen files will be created per alignment.  Their names are <outprefix>.<extension> and <outprefix>_bent.<extension>, where extension is fasta, txt, pdb, pir, msf, spt, or spdbv.  See information on the "-f" and "-b" parameters for information on how to set which files Matt creates.  In addition, a single log file will be generated when doing multiple alignments with a single call to Matt.  By default, only fasta, txt, pdb, and spt files are output.  <outprefix> is specified as a command line option.  In all files, proteins are listed in the input order, except for the assembly order section at the top of the txt files.  By default, the "bent" files will not be created.

The fasta, pir, and msf files contain the alignment the corresponding formats, using dashes (fasta and pir) or periods (msf) to indicate gap positions.

The txt is an interleaved visual alignment of the sequences of the three structures.  It also includes the assembly order, RMSD, number of core residues, raw score, p-value (For pairwise alignments), and the reference structure.  The reference structure is the one that is untransformed in the final pass of the algorithm.  It also is only slightly transformed in the output pdb files.  Other than this and the fact that it's one of the two structures used to calculate rotation angles in the final alignment, the reference structure has no special significance.

The pdb files contain 3D atomic coordinates of the structural alignment.

The spt files are RasMol/Jmol scripts that highlight regions in the common core.  To run the scripts with Jmol, just open the PDB files and drag the script to the Jmol window.  With RasMol, open the pdb file, and type "script <filename>.spt".  The core residues from each structure will be set to a different color.  The colors repeat after 10 structures.  If you get the error "Unable to allocate shade" when running the script in RasMol, try "select all" and then "color gray" and then run the script again.

The spdbv files are Swiss PDB Viewer / DeepView scripts.  They work much like the spt files, and use roughly the same set of 10 colors in the same order.  DeepView apparently does not support scripted highlighting by insertion code, so be aware that a few extra residues may be colored in PDB files with insertion codes.  Enabling renumbering ("-r") will get around this issue.

The bent files, enabled by "-b", contain the results generated before the final extension pass, which align the unbent structures and fills in gaps.  One bent file will be created for each unbent file, in the same file format.  The bent pdb files are the output structures generated by the first phase of the algorithm.  The source structures are broken apart halfway between adjacent fragments.  The RMSD in the bent txt file is the RMSD of the alignment of deformed structures, so should not be compared directly to the RMSD of other structural alignment algorithms.


Environment Variables:

Matt supports several environment variables to make command line usage simpler.  Use of these is completely optional.


MATT_SEARCH_PATH:  If a pdb input file specified by a relative path isn't found relative to the current directory, Matt then looks for it in MATT_SEARCH_PATH.  MATT_SEARCH_PATH can be a semi-colon delimited list to search though multiple directories.  These are searched through in order of occurrence.  Note that all variations of the name of a specified file are tried in each path before moving on to the next.  That is, if you specify "1plu", Matt will use "1plu.pdb" if it occurs earlier in the search path than a file named just "1plu".

MATT_PDB_PATH:  If a pdb input file specified by a relative path isn't found relative to the current directory or in MATT_SEARCH_PATH, Matt then looks for it in MATT_PDB_PATH.  The only difference between this and MATT_SEARCH_PATH is that it will also look for a subdirectory using the second and third letters of the name you give it, so if you specify "1tsp", Matt will find "<pdb path>/ts/pdb1tsp.ent.Z", for example, which is how the files are arranged on RCSB's FTP site.  If the name starts with "pdb" and has 7 characters, or seven characters before a period, Matt will use the fifth and sixth characters instead, so "pdb1tsp" and "pdb1tsp.ent" would both find the file as well.

MATT_PARAMS:  Contains default parameters to use.  Uses the exact same syntax as the command line itself.  Any command line parameter except pdb files and list files can be put in this variable.  Anything explicitly specified on the command line will override these values.  For example, if set to "-o temp -b", Matt will use an output prefix of "temp" when none is specified and will create bent alignment files by default.  If MATT_PARAMS is set to a set of invalid parameters, Matt may refuse to run until its value is fixed.


Multithreading:

Matt uses the OpenMP multithreading extensions when built with OpenMP enabled with a compliant compiler, such as gcc 4.2 and retail versions of Microsoft Visual Studio 2005.  Note that for compatibility reasons, the included makefile will not create a binary with OpenMP support.

Matt supports two types of division of labor between threads.  When a single multiple alignment is run, threads are assigned pairwise alignment.  This method only benefits multiple alignments.  The other method will assign threads complete alignments.  The former method works best for small numbers of alignments, or when a small set of alignments take significantly longer than the others.  The latter method works best when running on a large numbers of alignments, particularly when running many alignments of small numbers of structures or running on computers with more than 2 cores.  Matt will use the first method when a single alignment is specified on the command line, and the second method when multiple alignments are specified.


Unless otherwise specified, Matt will use OpenMP's default number of threads, which is generally one thread per CPU per core.  More information on that command line option (-t) is below.  Reducing the number of threads will also decrease memory usage

Matt reports the number of threads that are created, not the number of threads that are active.  Therefore, when running a pairwise alignment on a multi-core system, it may report multiple threads, even though only one is doing any work.

Running Matt with no parameters will display version and usage information.  If built with multithreading support, the version number will be followed by "OpenMP".


Command Line:
Matt [-o <outprefix>]* [-c <cutoff>] [-t <threads>] [-[rslbdvp][01]]*
     [-f <extension list>] [file[:<chains>]]* [-L <listfile>]

Square brackets indicate optional parameters, asterisks indicate parameters that can appear more than once, and angle brackets indicate that values should be substituted for the text.  None of those symbols should actually be used when running Matt.

<chains> must be in the format:
     <chain 1 id>:<model number>,<chain 2 id>:<model number> ...
     
Leaving off the ":<model number>" indicates the single unnumbered model (model 0) in PDB files with no numbered model, or model 1 in pdb files with numbered models.  Leaving off the chain id but keeping the colon means all chains with an id, if there are any, or the single unnamed chain if there aren't any.  Note that two colons will occur in a row when using this format immediately after the file name ("1tsp::1,:2").  When multiple chains come from the same model, the chains can be listed without commas, so "1tsp:ABC" and "1tsp:ABC:0" both mean chains A, B, and C from model 0.  Optionally, a subset of a chain's residues can also be listed.  "1tsp:A(1-100)B(1-50)" means two chains, one consisting of residues 1 to 100 of chain A and another consisting of residues 1 to 50 of chain B.  When using this format, chain ids must be explicitly specified, even when there are no named chains, so quotes may be needed if using the command line.  The numbers are either 1-based residue indices or PDB ids, depending on whether or not residue renumbering is enabled (See -r.  Renumbering is disabled by default).


Command line notes:

For options that don't take a space before their parameter (r, s, l, b, d, p, and v), giving the option with no parameter is equivalent to specifying a parameter of 1.  Also, the order of parameters is irrelevant, unless you use more than one "-o" parameter, and then only the position of the "-o" relative to listed pdb files matters.  All other options apply to all runs.  "-s", for example, which affects how pdb files will be read, affects both pdb files before and after the "-s" option is specified.  You can also combine multiple options with a single hyphen.  The following lines are equivalent, when there is no MATT_PARAMS environment variable:

Matt 1plu.pdb 1tsp.pdb -b1 -s1 -d1 -c 5.0 -o alignment
Matt -b1sdc 5.0 -o alignment 1plu.pdb 1tsp.pdb
Matt -b -o alignment 1plu 1tsp

Note that the last line removes all the parameters that the earlier lines set to their default values and takes advantage of the automatic search for files with the .pdb extension.  If there are files named both "1plu" and "1plu.pdb" in the same directory, the last line will load "1plu" instead of "1plu.pdb", unlike the first two.


Primary command line parameters:

-o <outprefix>:  Specifies prefix of output file names.  Unlike most other parameters, you can use -o more than once to run multiple jobs in parallel.  Every pdb file before the second -o is put in the first alignment, every one before the third goes into the second, and so on.  Running on multiple files at once will result in one log file being generated for each alignment, in place of some of the console output to avoid interlaced error messages.  Note that running this way causes threads to be assigned complete multiple alignments, as opposed to pairwise alignments that will be assembled into a single multiple alignment.  The difference is discussed in more detail in the Multithreading section above.  If no output prefix is specified, "MattAlignment" is used by default.


-L <listfile>:  Specifies a file containing a list of pdb files to load.  Each line of the file must specify a different file.  Chains, residues, and models can optionally be specified using the same syntax as with pdb files listed on the command line.  Leading and trailing white space is ignored.  Lines starting with "#" are treated as comments and skipped.  Lines with the format "-o outprefix" can also be used to specify more than one alignment, though no other commands are allowed in a list file.  Blank lines are allowed.


Secondary command line parameters:

-b[01]:  Disables or enables creation of bent files.  Disabled by default.

-c <cutoff>:  Sets the distance cutoff value, in Angstroms, for the final pass, which fills in some of the gaps.  Cutoff can be any non-negative floating point number.  This does not affect any of the bent files.  The default value is 5.0 angstroms, which is what was used in the paper.  A value of 0 prevents the last pass from running at all.  Note that there should be a space after the c.  P-values are calculated before the final extension pass, so the cutoff does no affect reported p-values.

-r[01]:  Disables or enables renumbering of all residues in all proteins.  Each protein will start from residue 1 and all residues will be numbered consecutively.  Insertion codes will be removed.  Note that all loaded residues will be given a number, so if SEQRES entries are loaded or some residues have no alpha carbons, the first residue used in the alignment may not be residue 1.  Disabled by default.

-l[01]:  Disables or enables renaming chains in the output pdb files.  When enabled, chains are first labeled by capital letters, then numbers, then symbols, then lowercase letters, and then the pattern repeats when there are over 90 chains.  The limit is due to the fact that chains in a pdb file can only have a single character label.  Chains are numbered according to the order they're specified in the command line.  Enabled by default.

-s[01]:  Disables or enables reading SEQRES lines in source pdb files.  When enabled, the program tries to align residues in the ATOM entries to residues in the SEQRES entries.  This allows detecting gaps between residues that would otherwise be assumed to be adjacent.  Fragments cannot cross over regions with no alpha carbon coordinates.  Note that residues with ATOM entries but no alpha carbons coordinate will always be loaded.  When enabled, unalignable residues without alpha carbons will still appear in output files.  -s also affects residue renumbering if -r is set.  Enabled by default.

-v[01]:  Sets verbosity of feedback to stdout or log file.  A value of 0 will only display errors and warnings, and 1 will display a list of chains as they are loaded.  Default value is 1.

-d[01]:  Disables or enables sending current progress to stderr.  Current progress is either how many pairwise alignments are completed when running a single multiple alignment, or how many multiple alignments have been completed when running more than one multiple alignment at once.  Enabled by default.

-p[01]:  Disables or enables partial alignments.  Has no effect on core residues and bent alignments.  Partial alignment can have a very significant performance impact, particularly on alignments of large numbers of structures with few core residues.  Enabled by default.

-t <thread count>:  Sets the number of threads Matt uses.  If not specified, Matt will use OpenMP's default number of threads, which is implementation dependent, though it is generally the number of threads a system is capable of running synchronously.  When not compiled with OpenMP support, a warning will be displayed and the option will be ignored. Note that there should be a space between t and the number of threads.

-f <extension list>:  Sets the output file formats to use.  Extension list is case sensitive and must be comma-delimited and have no spaces.  "all" is equivalent to listing all supported formats.  The default value is "pdb,fasta,pdb,spt".  The only other formats currently supported are msf, pir, and spdbv.  When using -f, you must specify the default file formats again if you want them to be created, so "-f msf" would create no pdb files, for example.


Notes:

Matt will list residues with no alpha-carbon coordinates in its sequence alignments, but will never align them.  The recommended way to unambiguously figure out which atom entries in the alignment files corresponds to which entry in the created pdb files is to enable renumbering, keeping in mind Matt includes HETATOM entries with alpha carbons in the alignment.  The numbers of the PDB files will then correspond exactly to the positions of residues in the alignment files.

When run on more than 26 chains, Matt will start using symbols and then lowercase letters for chain ids, which can cause issues with programs that read PDB files.  For large enough sets of chains, Matt will start using duplicate chain names.


Algorithm changes:

The realignment phase and final RMSD alignment now use two passes to get marginally tighter alignments.  Because of this change, Matt tends to give fewer core residues and lower RMSD with the same cutoff value.

A couple changes were made so that the partial alignment algorithm adds no more core residues.  Matt 1.0's final extension phase uses the average of <alpha carbon distance> for all aligned residue pairs instead of using their RMS distance.  In addition, Matt 1.0 repeatedly runs an RMSD alignment followed by the final extension phase until no more core residues are added.  This prevents any core residues from being added in the partial alignment phase.


Partial Alignments:

Matt's partial alignments are primarily intended to make visual inspection of generated linear alignment files easier.

Matt's partial alignment algorithm starts with the final 3D super imposition, computed as discussed on the paper.  For each unaligned region between sets of core regions, the algorithm goes through pairs of sets of structures in the original assembly order, using dynamic programming to find an optimal alignment of unaligned residues using a simple scoring function.  As partial alignments don't always have a well-defined sequence order, and the optimal aligned sets of residues from a subset of structures may not be a subset of the optimal aligned sets of residues on the full set of residues, the returned result is not guaranteed to be optimal, given the scoring function.  Some effort has been made to correct for this in subsequent passes by increasing the search space slightly.

The scoring function is just the sum of <rmsd cutoff> - <alpha carbon distance> for all aligned residue pairs.  Having a purely additive function is simpler in a dynamic programming framework.
