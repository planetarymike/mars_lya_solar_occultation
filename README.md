# Mars Lyman alpha Solar Occultation Routines

The included routines compute optical depths as a function of Doppler shift on a line of sight towards the Sun from a spacecraft embedded in the H corona of Mars. These optical depths can be used to extinct an assumed solar line shape and derive properties of the H corona. 

To use these routines, you will need a C++ compiler, a set of spacecraft altitudes and solar zenith angles, and the altitudes of the tangent point along the spacecraft line of sight towards the Sun. 

First, compile the code:

	make solar_LOS_profiles

To run the code, you'll need to specify an H exobase density and temperature, along with a filename (relative to the execution directory) pointing to a file of spacecraft positions. An example file including these positions is given in test_out.dat; this file was produced by the IDL routine sav2cppinput.pro, which you can inspect if you'd like to generate more such files. Then specify an output filename, again relative to the execution directory, where the output line of sight profiles should be stored. That's it!