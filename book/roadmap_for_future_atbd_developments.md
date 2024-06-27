# Roadmap for future ATBD development

The core algorithm of the Sea Ice Drift has been long used in other projects,
and there was as thorough investigation of the use of this algorithm for
swath-to-swath calculations as required by CIMR (paper ref). The future
development work is therefore largely focussed on the niche requirements
for CIMR test cards needed to thoroughly validate this algorithm.

## Expansion of core code to process all channels

The core cross-correlation algorithm used for the calculation of sea-ice
drift optimises multiple channels at once. This is currently able to optimise
up to 6 channels simultaneously, however, since we require 8 channels (Ku-band
and Ka-band, V and H-polarisations, and forward and back scans in all
permutations), this code needs expansion to ensure that it can run without
encountering memory issues. The current demonstration uses only Ka-band.

## Requirement for new test scenes

In order to calculate sea-ice drift, it is required to have two swaths
separated in time, between which features in the ice are cross-correlated
to determine the ice motion. In the initial validation, we were lacking any
such pair of test scenes, and therefore we made a simple test by shifting
the data from the Radiometric test scene by 3 pixels in x and 4 pixels in
y, as well as with a 24-hour offset.

...

## Investigation of optimum time offsets for swath combination

...