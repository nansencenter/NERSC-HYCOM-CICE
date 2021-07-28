# How to use the sections and transport tools

## Preparations
Whether you are extracting vertical sections of fields using **m2section** or computing transports using **m2transport** or **m2transport2** you need to define the start and end point of the sections in a file called *sections.in*. An example of such a file can be found in the folder **NERSC-HYCOM-CICE/hycom/MSCPROGS/Input/**. With *sections.in* defined, it was formerly necessary to first run `section_intersect`, whose purpose it is to set up the grid points which are needed to calculate transports and to interpolate section details. However, this is automatically called from the respective section tools executables and no longer needs to be done manually.

## Extracting sections
First define which variables to extracted in *extract.sec*.  An example of such a file can be found in the folder **NERSC-HYCOM-CICE/hycom/MSCPROGS/Input/**. Then run: 

`  m2sections [list of .b hycom files to precess] `

## Scalar transport

The “scalartransport.in” file can be used to make the routine m2transport2 calculate transport for additional tracer fields,

such as for nutrients in biological models. The format is :

’temp ’ 6. 3987000  
’saln ’ 0. 1.

The first column denotes the variable to calculate the scalar transport for.

The second column is an offset value which is subtracted from the variable when the transport is calculated.

Finally, the third column is a scale factor, which is multiplied with the final scalar transport estimate.

The example on the first line above will calculate temperature transport relative to 6 degrees Celsius,

with a scale factor of 3987000 ≈ cp × ρ, where cp is the heat capacity of sea water and ρ is water density.

Line 1 would therefore give a heat transport estimate relative to 6 degrees celsius.

## Usage, requirements & outputs of the section and transport tools
Here a brief explanation of each tool, its usage, options, required input files and outputs will be described.

### m2transport
**Description:**
   This script will calculate barotropic transports across a section. It will also try to calculate ice transports. As it only needs to read the 2D field barotropic velocity, it is really fast. You will need to specify the sections in the file *sections.in*.

**Usage:**

`m2transport [-nosec] <files>`,

   where the optional argument:
   `-nosec`    specifies that section positions don't need to be calculated again, convenient on multiple passes.

**Required files:**
- regional.grid.a
- regional.grid.b
- regional.depth.a
- regional.depth.b
- sections.in
   These files need to be in your current working directory when executing the `m2transport` routine.

**Outputs:**
- transports.nc       # these are the transports outputs
- section001.dat      # not sure if this is supposed to be the actual data
- transports001.dat   # not sure what this is, looks like grid indices and 1-or-0 vales in two columns
- tst001.nc           # this looks like a bathymetry mask for the section(s)
- .zoneinfo           # just says "ZONE>1"
   Some of these outputs require clarification from someone who is more familiar with them.

### m2transport2
**Description:**
   This script will calculate transports across a section, but with a more fine-grained approach than `m2transport`. You will need to specify the sections in the file *sections.in*. You will also have to define how to calculate transports in *transports.in*.

**Usage:**

`m2transport2 [-nosec] file1 file2 ...`,

   where the optional argument:
   `-nosec`    specifies that section positions don't need to be calculated again, convenient on multiple passes.

**Required files:**
- regional.grid.a
- regional.grid.b
- regional.depth.a
- regional.depth.b
- sections.in
- transport.in
   These files need to be in your current working directory when executing the `m2transport` routine.

**Outputs:**
- transports2.nc      # these are the transports outputs
- section001.dat      # not sure if this is supposed to be the actual data
- transports001.dat   # not sure what this is, looks like grid indices and 1-or-0 vales in two columns
- tst001.nc           # this looks like a bathymetry mask for the section(s)
   Some of these outputs require clarification from someone who is more familiar with them.

> **An important note:**
> With `m2transport2` you can also calculate scalar transports in addition to volume transports. The routine `m2transport2` will look for a file *scalartransports.in* which can be used to calculate these transports, such as heat transport. However, when I, Harry, tried to calculate scalar transports with the *scalartransport.in* file present no additional output compared to standard m2transport2 is created. This behaviour occurs even when the scalar transport flag is turned on in the compilation makefile. Feel free to add instructions here if you know how to overcome this issue.
