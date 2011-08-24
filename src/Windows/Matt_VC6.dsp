# Microsoft Developer Studio Project File - Name="Matt" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Matt - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Matt_VC6.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Matt_VC6.mak" CFG="Matt - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Matt - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Matt - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Matt - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\Release"
# PROP Intermediate_Dir "..\..\Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /Gr /MD /W3 /GX /O2 /Ob2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "Matt - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\Debug"
# PROP Intermediate_Dir "..\..\Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Matt - Win32 Release"
# Name "Matt - Win32 Debug"
# Begin Source File

SOURCE=..\AssemblyOrder.c
# End Source File
# Begin Source File

SOURCE=..\AssemblyOrder.h
# End Source File
# Begin Source File

SOURCE=..\chain.c
# End Source File
# Begin Source File

SOURCE=..\chain.h
# End Source File
# Begin Source File

SOURCE=..\Extend.c
# End Source File
# Begin Source File

SOURCE=..\Extend.h
# End Source File
# Begin Source File

SOURCE=..\FileReader.c
# End Source File
# Begin Source File

SOURCE=..\FileReader.h
# End Source File
# Begin Source File

SOURCE=..\Matt.c
# End Source File
# Begin Source File

SOURCE=..\MultipleAlignment.c
# End Source File
# Begin Source File

SOURCE=..\MultipleAlignment.h
# End Source File
# Begin Source File

SOURCE=..\MultipleAlignmentOutput.c
# End Source File
# Begin Source File

SOURCE=..\MultipleAlignmentOutput.h
# End Source File
# Begin Source File

SOURCE=..\OctTree.c
# End Source File
# Begin Source File

SOURCE=..\OctTree.h
# End Source File
# Begin Source File

SOURCE=..\pdb.c
# End Source File
# Begin Source File

SOURCE=..\pdb.h
# End Source File
# Begin Source File

SOURCE=..\Protein.c
# End Source File
# Begin Source File

SOURCE=..\Protein.h
# End Source File
# Begin Source File

SOURCE=..\RMSD.c
# End Source File
# Begin Source File

SOURCE=..\RMSD.h
# End Source File
# Begin Source File

SOURCE=..\Score.c
# End Source File
# Begin Source File

SOURCE=..\Score.h
# End Source File
# Begin Source File

SOURCE=..\secondary.c
# End Source File
# Begin Source File

SOURCE=..\secondary.h
# End Source File
# Begin Source File

SOURCE=..\util.c
# End Source File
# Begin Source File

SOURCE=..\util.h
# End Source File
# Begin Source File

SOURCE=..\Vector.c
# End Source File
# Begin Source File

SOURCE=..\Vector.h
# End Source File
# End Target
# End Project
