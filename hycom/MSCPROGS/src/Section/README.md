## How to use the Sections and transport tools
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
