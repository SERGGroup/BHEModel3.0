This folder is used with activity coefficient models beginning with Chapter 11.

Note that unifac.m is used with the files in either of the subfolders UnifacVLE or UnifacLLE.

All activity models expect the composition and parameter values to be passed to the function.
Some routines (e.g. UNIFAC, UNIQUAC, NRTL) also expect the temperature. Check the routine for units.

Some files are multicomponent, some are binary.
Binary files generally can evaluate a whole set of x values simultaneously for isothermal data.

See examples.m for example function calls or the 'caller' m-files in the UNIFAC folders.
