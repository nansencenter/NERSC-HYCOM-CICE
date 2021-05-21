# How to use the sections and transport tools

## Preparations
Wether you are extraction vertical sections of fields using **m2section** or computing transports using **m2transport** or **m2transport2** you need to define the start and end point of the sections in a files called **sections.in**.  An example of such a file can be found in the folder **NERSC-HYCOM-CICE/hycom/MSCPROGS/Input/**. When sections.in is defined, run **section_intersect**, whose purpose is to set up the grid points which are needed to calculate transports and to interpolate section details.

## Extracting sections
First define which varibales to extract in **extract.sec**.  An example of such a file can be found in the folder **NERSC-HYCOM-CICE/hycom/MSCPROGS/Input/**. Then run: 

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
