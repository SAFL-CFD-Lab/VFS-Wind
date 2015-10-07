#!MC 1400
# Created by Tecplot 360 build 14.0.0.25097
$!VarSet |MFBD| = './'
$!READDATASET  '"../TSR5.0/ave/Turbine_AL06_00" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = XYLINE
  VARNAMELIST = '"time" "angle" "angvel_axis" "Force_axis" "Torque_fluid" "Torque_generator" "Uref" "Ud" "TSR"'
$!READDATASET  '"../TSR8.0/Turbine_AL06_00" '
  READDATAOPTION = APPEND
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = XYLINE
  VARNAMELIST = '"time" "angle" "angvel_axis" "Force_axis" "Torque_fluid" "Torque_generator" "Uref" "Ud" "TSR"'
$!ALTERDATA 
  EQUATION = '{CP}=v3*v5/(0.5*3.14*0.5*0.5*(1**3))'
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer4'
  COMMAND = 'Integrate [1-2] VariableOption=\'Average\' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar=9 Absolute=\'F\' ExcludeBlanked=\'F\' XVariable=1 YVariable=2 ZVariable=3 IntegrateOver=\'Cells\' IntegrateBy=\'Zones\' IRange={MIN =1 MAX = 0 SKIP = 1} JRange={MIN =1 MAX = 0 SKIP = 1} KRange={MIN =1 MAX = 0 SKIP = 1} PlotResults=\'F\' PlotAs=\'Result\' TimeMin=0 TimeMax=0'
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer4'
  COMMAND = 'Integrate [1-2] VariableOption=\'Average\' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar=9 Absolute=\'F\' ExcludeBlanked=\'F\' XVariable=1 YVariable=2 ZVariable=3 IntegrateOver=\'Cells\' IntegrateBy=\'Zones\' IRange={MIN =1 MAX = 0 SKIP = 1} JRange={MIN =1 MAX = 0 SKIP = 1} KRange={MIN =1 MAX = 0 SKIP = 1} PlotResults=\'F\' PlotAs=\'Result\' TimeMin=0 TimeMax=0'
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer4'
  COMMAND = 'SaveIntegrationResults FileName=\'./TSR.txt\''
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer4'
  COMMAND = 'Integrate [1-2] VariableOption=\'Average\' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar=10 Absolute=\'F\' ExcludeBlanked=\'F\' XVariable=1 YVariable=2 ZVariable=3 IntegrateOver=\'Cells\' IntegrateBy=\'Zones\' IRange={MIN =1 MAX = 0 SKIP = 1} JRange={MIN =1 MAX = 0 SKIP = 1} KRange={MIN =1 MAX = 0 SKIP = 1} PlotResults=\'F\' PlotAs=\'Result\' TimeMin=0 TimeMax=0'
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer4'
  COMMAND = 'SaveIntegrationResults FileName=\'./CP.txt\''
$!RemoveVar |MFBD|
