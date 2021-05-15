A GNU Octave script to calculate the transmission, reflection and insertion phase delay 
of single or multilayer radomes.

The following files must be in the same folder/directory:

Multilayer_Radome_Analysis_MASTER_V1_1.m

sub_wall.m

cmult.m

pol2rect.m

rect2pol.m

complex_sqr_root.m

The code is written in GNU OCTAVE but should be compatible with Matlab.  If Matlab complains,
only minor changes are likely to be required.

See the PDF for instructions.

November 2020

Paul Klasmann


New in V1.2

Added section with an additional comparison of the same C-sandwich radome simulated with a homebrew 1D FDTD
simulation using the method described by EMPossible (Prof Rumph).  No changes to the analytical method.

May 2021

Paul Klasmann