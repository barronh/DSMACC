TUV4.5
last build Sep 2007 by S. Madronich

To run in UNIX/Linux, unpack the source code zip or tar file, then:
compile using "make" 
run using "tuv"

To run in Windows (pre-compiled executable):
After unpacking the zip file, double-click on the icon Run_tuv.bat and follow
the on-screen instructions.

TABLE OF CONTENTS:
******************  1.   SELECTING AN INPUT FILE:
******************  2.   CHANGING THE VALUES
******************  3.   SAVING THE NEW INPUT FILE
******************  4.   OUTPUT CONTROL
******************  5.   FREQUENTLY ASKED QUESTIONS
******************  6.   SOURCE FILES AND SUBROUTINES
******************  7.   DATA FILES


******************  1.   SELECTING AN INPUT FILE:

Input files are used to modify the conditions for the radiation 
calculations, and to select different desired outputs. Several sample
input files are provided.  They can be modified by the user, saved
(with a different name), and re-called.  The default input files are:

defin1 = suggested for calculations of UV spectral irradiances 
at the Earth's surface. The wavelength range is 280-420 nm in 1 nm steps, 
and the outputs are spectral irradiance and several weighted irradiances:
UV-B (280-315 nm), UV-B*(280-320 nm), UV-A (315-400 nm), DNA damage,
CIE erythema, and the UV-Index. The radiation is calculated with a 
2-stream delta-Eddington code, and values are the sum of the direct sun
and down-welling diffuse radiation (dirsun = difdn = 1, difup = 0).

defin2 = suggested for calculations of spectral actinic fluxes and photolysis
rate coefficients (J-values), at 0.5 km above sea level.  The wavelength 
range is 121-735 nm in 156 non-equally spaced intervals (standard WMO grid).
Wavelengths shorter than ca. 290 nm are important if the output altitude
(zout) is set to the stratosphere or higher. J-values for several different 
reactions are given.  The radiation is calculated with a 4-stream discrete
ordinates code, and values are the sum of all directions: Direct sun, 
down-welling diffuse, and up-welling diffuse (dirsun = difdn = difup = 1).

defin3 = Output for use with the NCAR Master Mechanism. 

defin4 = Illustrates the possible outputs that can be obtained by use of the
logical switches (lirrad, laflux, lmmech, lrates, ljvals all set to .TRUE.) and
by using non-zero values for the integers isfix, ijfix, iwfix, itfix, izfix.
NOTE:  The tables of values (e.g. spectral irradiance vs. altitude and wavelength)
can get quite large, therefore only very coarse resolution is specified for this
illustration of possible outputs.  
CAUTION: The actual values of radiative quantities are not very accurate when
using such low wavelength and altitude resolution.

Sample output, obtained with each of these default input files (defin1-4) can
be found in the subdirectory SAMPLES.

******************  2.  CHANGING THE VALUES

The values of various input and control variables can be changed.  
By typing ?variable, you will get a few lines of help on that variable.
To change a value, type the name of the variable shown in the table. 
(NOTE:  The variable names are case-sensitive)

If the variable requires a number, you will be prompted to enter the new value.

If the variable requires a character string (e.g. a filename), you will be 
prompted to enter the new value.

If the variable is of type LOGICAL, it will switch (toggle) between True and
False.

If the variable is nms, you will be taken to a new menu of available
weighting functions (action spectra).  In the new menu, type the number of the
weighting function to switch  on/off (True/False).

If the variable is mnj, you will be taken to a new menu of available
photolysis reactions. In the new menu, type the number of the
photolysis reaction to switch  on/off (True/False).

When done with the changes, press <enter> to continue.

******************  3.  SAVING THE NEW INPUT FILE

You will be asked if you want to save your modifications to the input file.
You can give any name to the new input file (default = usrinp).
Next time you run the program, you can read this input file rather than
starting from one of the default inputs.

******************  4.  OUTPUT CONTROL

Output will be written to files in the directory ONE LEVEL ABOVE the
directory containing the executable.  In windows, this means the output
will appear in the directory that contains the file Run_tuv.bat.

You can change the name of the output file by changing the variable name
"outfil" in the main tabular menu.  The default is "usrout"

The program will add the extension .txt to the output name, e.g.
    usrout.txt

You can also send the output to the screen instead of a file, but this
is not recommended if printing large tables.

Output options are enabled by switching (True/False) several variables:
lirrad, laflux, lmmech, lrates, ljvals.

Detailed output tables can also be created by setting non-zero integer values 
for isfix, ijfix, iwfix, itfix, izfix.  These create "slices" through three-
dimensional matrices of results.  For example, the spectral irradiance is
a function of 3 coordinates:  time (t), wavelength (w), and altitude (z).
If you set itfix = 2, the output will have a table of spectral irradiance
as a function of wavelength and altitude, for the 2nd time step.  If you
set iwfix = 10, the output will be the spectral irradiance at the 10th 
wavelength (as set with wstart, wstop, nwint), as a function of altitude and
time.

Sample output files are provided in subdirectory SAMPLE.  The file
SAMPLE/usrout4 shows the various available outputs.

The code will also output a file 'tuvlog.txt' which records the inputs used
in the calculations.

******************  5.  FREQUENTLY ASKED QUESTIONS

a)  What is the difference between zstart and zout? 
The variable zstart is the elevation, in km above sea level,  of the atmosphere.
If zstart = 0, the surface is at sea level.
On the other hand, zout is the altitude (km, above sea level) for which you want
output.  For example, suppose you are flying in an airplane 500 m above the
ground in Boulder, Colorado.  Boulder is located 1.7 km above sea level, so
zstart = 1.7, while zout = 2.2 (1.7 + 0.5 km).

b)  What is the purpose of zaird, ztemp?
For some photolysis reactions, the absorption cross sections and quantum
yields may depend on temperature and/or pressure.  If the local (at z = zout)
temperature and pressure are known, and if they are different from the
US Standard Atmosphere, they can be put in here.  These values do not affect
the radiation field.

c)  What is the purpose of the variables dirsun, difdn, and difup?
They allow the user to compute separately the direct solar beam component
(dirsun), the diffuse down-welling or scattered radiation (difdn), or the
diffuse up-welling radiation (difup, which includes reflections from the 
surface and, if zout > zstart, from the atmosphere below).  Examples:
- down-welling irradiance or actinic flux:  dirsun = difdn = 1, difup = 0
- total actinic flux: dirsun = difdn = difup = 1
- net irradiance:  dirsun = difdn = 1, difup = -1

e)  Can I calculate the diffuse/direct ratio?
Not in a single calculation.  You must do two calculations:
the first with dirsun = 1, difdn = difup = 0
the second with  dirsun = 0, difdn = 1, difup = 0
Then you can take the ratio of the results.

f)  How can I get the extraterrestrial solar spectral irradiance or 
actinic flux?
Set 
      zout = zstop 
      lzenith = T  
      tstart = 0
      nt = 1
      dirsun = 1
      difdn = 0
      difup = 0
Then for spectral irradiance use
     lirrad = T  (for W m-2 nm-1)
or for actinic flux use:
     laflux = T ( for quanta cm-2 s-1 nm-1)

g) How can I calculate atmospheric transmission down to a particular level?
Calculate the spectral irradiance or spectral actinic flux at desired
altitude (zout), then separately calculate the extraterrestrial spectral
irradiance or actinic flux (see FAQ-f), then divide the results.

h) Are the wavelengths defined IN-VACUUM or IN-AIR?
The TUV model works strictly with IN-VACUUM wavelengths.  This includes
all spectral data, such as extraterrestrial spectral irradiance, 
absorption cross sections, quantum yields, and action spectra.
However:  Some users may want spectral irradiances or spectral actinic
fluxes on a wavelength scale defined IN-AIR. In this case, you must edit
the main program (TUV.f) and set the logical variable lrefr=.TRUE.
If lrefr=.TRUE., 
- the specified wavelength grid will be assumed to be IN-AIR
- the calculations will be done IN-VACUUM (the code will shift the grid)
- the results will be reported on the IN-AIR grid (the code will shift back).
The wavelength shift will be calculated with an index of refraction appropriate
for the main output altitude, zout.  Note that the wavelength scale at other
altitudes will not be correct and should not be used.
Control of this option is not available from the interactive input table.

******************  6.   SOURCE FILES AND SUBROUTINES

The source code (Fortran 77) is contained in several files, each of which
may have a number of related subroutines:

TUV.f = main program 

params = include file, which sets dimensioning parameters and some constants

functs.f 
* This file contains the following user-defined fortran functions:
*     fery
*     fo3qy
*     fsum
*     futr

grids.f 
* This file contains the following subroutine, related to setting up
* grids for numerical calculations:
*     gridw
*     gridz
*     gridt
*     gridck

la_srb.f 
* This file contains the following subroutines, related to the calculation
* of radiation at Lyman-alpha and Schumann-Runge wavelengths:
*     la_srb
*     lymana
*     schum
*     effxs
*     calc_params
*     init_xs
*     sjo2   
* and the following function
*     chebev

numer.f
* This file contains the following subroutines, related to interpolations
* of input data, addition of points to arrays, and zeroing of arrays:
*     inter1
*     inter2
*     inter3
*     inter4
*     addpnt
*     zero1
*     zero2

odo3.f
* Compute ozone optical depths.

odrl.f
* Compute Rayleigh optical depths.

orbit.f
* This file contains the following subroutines, related to the orbit and
* rotation of the Earth:
*     calend
*     sunae

qys.f
* This file contains subroutines used for calculation of quantum yields for 
* various photoreactions:
*     qyacet - q.y. for acetone, based on Blitz et al. (2004)

rdetfl.f
* This file contains the following subroutines, related to reading the
* extraterrestrial spectral irradiances:
*     rdetfl
*     read1
*     read2

rdinp.f
* This file contains the following subroutines, related to reading
* simple input parameters from an input file, and interactive control. 
*     rdinp
*     write1
*     readin
*     chkval
*     newval
*     gethlp
*     select
*     atrim

rdxs.f
* This file contains the following subroutines, related to reading the
* absorption cross sections of atmospheric gases:
*     rdno2xs
*     rdo2xs
*     rdo3xs
*       o3xs_mm
*       o3xs_mal
*       o3xs_bass
*     rdso2xs

rtrans.f
* This file contains the following subroutines, related to the solution of
* the equation of radiative transfer in multiple homogeneous layers.
*     rtlink
*     ps2str
*        tridag
*     psndo
*        asymtx
*        chekin
*        fluxes
*        lepoly
*        pravin
*        prtinp
*        prtint
*        qgausn
*        setdis
*        setmtx
*        soleig
*        solve0
*        surfac
*        solvec
*        upbeam
*        zeroal
*        zeroit
*        errmsg
*        sgbco
*        sgbfa
*        sgbsl
*        sgeco
*        sgefa
*        sgesl
*        saxpy
*        sscal
*        sswap
*        t665d
*        t665r
* and the functions
*        dref
*        ratio
*        wrtbad
*        wrtdim
*        tstbad
*        sasum
*        sdot
*        isamax
*        d1mach
*        r1mach

rxn.f
* This file contains the following subroutines, related to reading/loading
* the product (cross section) x (quantum yield) for photo-reactions:
*     r01 through r48
*     r101 through r115

savout.f
* This file contains the following subroutines, related to saving and writing
* some specific outputs:
*     saver1
*     saver2
*     outpt1
*     outpt2

setaer.f
* specify aerosol vertical profile and wavelength-dependent optical 
* properties

setalb.f
* specify surface reflectivity (albedo)

setcld.f
* specify cloud vertical profile and wavelength-dependent optical 
* properties

setno2.f
* specify NO2 vertical profile and wavelength-dependent optical depths 

setno2.f
* specify SO2 vertical profile and wavelength-dependent optical depths 

seto2.f
*  specify O2 vertical profile and wavelength-dependent optical depths
*  optical depths will in Lyman-alpha and Schumann-Runge bands will be
*  overwritten in subroutine la_srb.f

setsnw.f
* specify snowpack depth and wavelength-dependent optical properties
* Snowpack parameters are set manually in setsnw, not interactively

setso2.f
* specify SO2 vertical profile and wavelength-dependent optical depths 

sphers.f
* This file contains the following subroutines, related to the
* spherical geometry of the Earth's atmosphere
*     sphers
*     airmas

swbiol.f
* This file contains the following subroutines, related to specifying 
* biological spectral weighting functions:
*     swbiol

swchem.f
* This file contains the following subroutines, related to specifying 
* chemical spectral weighting functions (cross sections x quantum yields)
*     swphys

swphys.f
* This file contains the following subroutines, related to specifying 
* physical spectral weighting functions:
*     swphys

vpair.f
* Specify vertical profile of air density (molec cm-3).

vpo3.f
* Specify vertical profile of ozone (molec cm-3).

vptmp.f
* Specify vertical profile of temperature (Kelvin).

wshift.f
* Shift wavelengths between vacuum and air scales.


******************  7.   DATA FILES

The main directory contains the following data files:
*   defin1,defin2,defin3,defin4 = default input files
*   helpin = data file for interactive help

Most other data is contained in the following directories
DATAE1 = data related to the atmospheric environment (e.g O3 profile)
DATAS1 = spectral weighting functions (e.g. biological action spectra)
DATAJ1 = cross sections and quantum yields for photolysis reactions

********************** END of README.txt ********************************
