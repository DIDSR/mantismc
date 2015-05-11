# Microsoft Developer Studio Project File - Name="mantis" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=mantis - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "mantis.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "mantis.mak" CFG="mantis - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "mantis - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "mantis - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "mantis - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt /fast
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "mantis - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /compile_only /debug:full /fltconsistency /libs:dll /nologo /traceback /warn:argument_checking /warn:nofileopt /warn:nouninitialized
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# SUBTRACT CPP /Fr
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# SUBTRACT LINK32 /profile /map /nodefaultlib

!ENDIF 

# Begin Target

# Name "mantis - Win32 Release"
# Name "mantis - Win32 Debug"
# Begin Group "inputs"

# PROP Default_Filter ""
# Begin Source File

SOURCE=".\test_pc\abs-vacuum.dm2"
# End Source File
# Begin Source File

SOURCE=.\test_pc\backing.dm2
# End Source File
# Begin Source File

SOURCE=.\test_pc\csi.dm2
# End Source File
# Begin Source File

SOURCE=.\test_pc\diode.dm2
# End Source File
# Begin Source File

SOURCE=.\test_pc\dixid00.geo
# End Source File
# Begin Source File

SOURCE=.\test_pc\gas.dm2
# End Source File
# Begin Source File

SOURCE=.\test_pc\mantis.job
# End Source File
# Begin Source File

SOURCE=.\test_pc\mantis_demo_CsI.cd2
# End Source File
# Begin Source File

SOURCE=.\test_pc\mantis_demo_CsI.dd2
# End Source File
# Begin Source File

SOURCE=.\test_pc\mantis_demo_CsI.md2
# End Source File
# Begin Source File

SOURCE=.\test_pc\mantis_demo_CsI.pen
# End Source File
# Begin Source File

SOURCE=.\test_pc\mantis_demo_CsI.sd2
# End Source File
# Begin Source File

SOURCE=.\test_pc\mantis_demo_CsI.xd2
# End Source File
# Begin Source File

SOURCE=".\test_pc\pseudo-vacuum.dm2"
# End Source File
# Begin Source File

SOURCE=.\test_pc\substrate.dm2
# End Source File
# End Group
# Begin Group "outputs"

# PROP Default_Filter ""
# End Group
# Begin Group "gnuplots"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\test_pc\tallyEnergyDepositionPulseHeightSpectrum.gpl
# End Source File
# Begin Source File

SOURCE=.\test_pc\tallyParticleCurrentSpectrum.gpl
# End Source File
# Begin Source File

SOURCE=.\test_pc\tallyParticleTrackStructure.gpl
# End Source File
# Begin Source File

SOURCE=.\test_pc\tallyPhaseSpaceFile.gpl
# End Source File
# Begin Source File

SOURCE=".\test_pc\tallySpatialDoseDistrib-1D.gpl"
# End Source File
# Begin Source File

SOURCE=".\test_pc\tallySpatialDoseDistrib-2D.gpl"
# End Source File
# End Group
# Begin Source File

SOURCE=.\DETECT2.F90
# End Source File
# Begin Source File

SOURCE=.\mantis.f
DEP_F90_MANTI=\
	".\penaux.f"\
	".\penelope.f"\
	".\pengeom.f"\
	".\penvared.f"\
	".\sourceBoxIsotropicGaussSpectrum.f"\
	".\sourcePhaseSpaceFile.f"\
	".\tallyDetect2.f"\
	".\tallyEnergyDepositionPulseHeightSpectrum.f"\
	".\tallyParticleCurrentSpectrum.f"\
	".\tallyParticleTrackStructure.f"\
	".\tallyPhaseSpaceFile.f"\
	".\tallySpatialDoseDistrib_c1.f"\
	".\tallySpatialDoseDistrib_c2.f"\
	".\tallySpatialDoseDistrib_c3.f"\
	".\timing.f"\
	
NODEP_F90_MANTI=\
	".\Release\mtpar.mod"\
	

!IF  "$(CFG)" == "mantis - Win32 Release"

# ADD F90 /inline:size

!ELSEIF  "$(CFG)" == "mantis - Win32 Debug"

!ENDIF 

# End Source File
# End Target
# End Project
