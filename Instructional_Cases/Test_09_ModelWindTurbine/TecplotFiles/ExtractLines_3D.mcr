#!MC 1400
# Created by Tecplot 360 build 14.0.0.25097
$!VarSet |MFBD| = '../TSR4/ave'
$!READDATASET  '"|MFBD|/Result002000-avg.plt" '
#$!READDATASET  '"/safl/projects/fsotiropoulos/xyang/channel_inflow_0-3H1/Result004200-avg.plt" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"X" "Y" "Z" "U" "V" "W" "uu" "vv" "ww" "uv" "vw" "uw" "Nv"'

$!VarSet |X0|=0.0
$!VarSet |Y0|=0.0867
$!VarSet |Z0|=0.0

$!VarSet |Y1|=0.01
$!VarSet |Y2|=0.3
$!VarSet |X1|=-0.6
$!VarSet |X2|=0.6
$!VarSet |Z1|=0
$!VarSet |Z2|=3

$!VarSet |NX|=241
$!VarSet |NY|=30
$!VarSet |NZ|=21
#$!VarSet |NZ|=1

$!VarSet |dZ|=0.15
$!VarSet |dZ|/=1.5

$!VarSet |Case|='TSR4'

$!VarSet |Z|=|Z0|

$!LOOP |NZ|


$!EXTRACTFROMPOLYLINE 
  EXTRACTTHROUGHVOLUME = YES
  EXTRACTLINEPOINTSONLY = NO
  INCLUDEDISTANCEVAR = NO
  NUMPTS = |NY|
  EXTRACTTOFILE = YES
  FNAME = '|MFBD|/Profile_x|X0%3.3f|_z|Z%3.3f|_|Case|.dat'
  RAWDATA
2
|X0| |Y1| |Z|
|X0| |Y2| |Z|


$!VarSet |Z|+=|dZ|
$!ENDLOOP

$!RemoveVar |MFBD|
