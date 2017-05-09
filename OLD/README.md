# PPM_Lagrange
Solves the Lagrange fluid equation using the Peicewise Parabolic Method (PPM).

To run the file just call the command 

python MakeFile.py

You may need open up MakeFile.py and change the first few variables to suit your system. There are two presets one called "PC" and another call "MAC". The PC presets assume you are running Cygwin. The "MAC" presests assume you are using a bash terminal.

If you would like to compile the Exact.f90 file on its own, you must do this from the basic directory, NOT the "RP-Exact" subdirectory. This is simply because of the architecture of MakeFile. This should not cause issues though.
