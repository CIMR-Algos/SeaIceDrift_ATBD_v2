# Background and justification of selected algorithm

## The Continuous Maximum Cross-Correlation method (CMCC)

As for many other motion tracking methods applied in geophysics
{cite:p}`nsidc:1998:drift-summary,notarstefano:2007:sst-motion`, sea ice
drift is tracked from a pair of images, with a block-based strategy.
Each block (*feature*, *sub-image*, \...) is composed of a limited
ensemble of pixels from the first image (the *reference* block) and
centred at a tracking location, for which the most similar block in the
second image (the *candidate* block) is looked for. The degree of
similarity is assessed by a metric, which often is the correlation
coefficient between the reference and candidate blocks. The maximum
correlation indicates the best match and the two-dimensional offset
between the centre points of the two blocks is the drift vector.

The description given above applies to the well known Maximum Cross
Correlation (MCC) technique which has been successfully applied by many
investigators {cite:p}`ninnis:1986:avhrr-drift,emery:1991:avhrr-drift,kwok:1998:ssmi85-drift,haarpaintner:2006:qscat-drift,notarstefano:2007:sst-motion,schmetz:1993:meteosat-wind`,
among others). In the MCC, wthe search for the best candidate block is
*discrete* and *exhaustive*. Discrete since the offsets between the
centre points are in whole number of pixels and exhaustive because all
candidate blocks are evaluated before the best can be chosen.

The same description applies to the Continuous Maximum Cross Correlation technique (CMCC)
method introduced by {cite:t}`lavergne:2010:cmcc-jgr`. CMCC is the strategy developed and
implemented first for the low resolution sea ice drift product of the EUMETSAT OSI SAF.
Conversely to the MCC, the search for the best candidate vector is
performed in a continuous manner over the two-dimensional plane and, as
a consequence, the search algorithm is not exhaustive.

```{figure} ./static_imgs/CMCC_vs_MCC.png
--- 
name: fig_cmcc_vs_mcc
---
Example se-ice displacements from AMSR-E (37 GHz H and V channels) over the Beaufort Sea and Canadian Basin from 29 to 31 January 2008. The product was processed with the (left) Continuous Maximum Cross Correlation (CMCC) and (right) Maximum Cross Correlation (MCC) method from the same satellite images. On the MCC product, zero-length vectors are depicted with a small square symbol while tiny arrows are used for the CMCC product. The MCC product clearly exhibits a quantization noise. Such an artifact is not present in the CMCC drift data set (reproduced from {cite:p}`lavergne:2010:cmcc-jgr`).
```

Although more complex to implement and, by nature, potentially less
robust than the MCC technique, the CMCC has the advantage of removing
the *quantization noise* which has hindered the retrieval of smooth
motion vector fields when the time span between the images is shortened
{cite:p}`ezraty:2007:amsr-pum,haarpaintner:2006:qscat-drift,kwok:1998:ssmi85-drift`.
**This is particularly important when the CIMR Level-2 product aims at sub-daily drift detection.**

## Choice of Ku and Ka as main microwave frequencies

As introduced above, sea-ice motion tracking from pairs of satellite images does not require specific microwave frequencies to work. As long as
the frequency allows a mostly unperturbed view of the sea-ice surface (window channels), the accuracy of the sea-ice drift retrieval does not
depend on the frequency per se, but rather on the spatial resolution achieved. Imagery with higher spatial resolution should result in higher
accuracy of the motion vectors.

This is the reason why the primary microwave frequencies for the Level-2 sea-ice drift product are the Ku (18.7 GHz) and Ka (36.5 GHz) channels
that will achieve spatial resolution of 5 km and better. Experience from the EUMETSAT OSI SAF near-real-time and climate processing confirms that
the Ka-band imagery of AMSR2 provides the best results. In CIMR, the spatial resolution of the Ku-band imagery will be only slightly worse than at
Ka-band, so that it should be included in the main processing. {cite:t}`kwok:2008:summer-drift` showed that the Ku-band imagery of AMSR2 provided
valuable sea-ice motion information during the summer melt season, but {cite:t}`lavergne:2021:s2s` did not find better performance of summer
sea-ice motion retrievals from Ku compared with Ka.

The lower frequency channels of CIMR (L, C, X) could possibly contribute to sea-ice motion tracking, but two challenges are their coarser spatial
resolution and the lack of intensity patterns to track from one image to the next (Ku and Ka imagery are sensitive to snow and sea-ice type, which has
large variability across the polar sea-ice and is stable over the tracking window of 0,5 to 2 days). This ATBD thus focus on Ku and Ka imagery but leaves
the door open for later inclusion of at least C and X.

## Exploitation of the forward and backward scan

In the past (e.g. {cite:t}`haarpaintner:2006:qscat-drift,kwok:1998:ssmi85-drift,girard-ardhuin:2012:drift`, among others), one sea-ice motion vector field would have been processed
for each microwave channel as input, i.e. one vector field from the Ku-H band, one from Ku-V, etc... the different vector fields would then be merged together a-posteriori.

{cite:t}`lavergne:2010:cmcc-jgr` introduced a different approach where a single sea-ice drift motion field is processed in one go from all the available imagery channels. This *implicit* merging
it implemented at the core of the motion tracking algorithm, by solving for the maximum value of the sum of the cross-correlation of several imagery channels, instead of just one imagery channel.

For CIMR, this can be further extended by considering the imagery from *forward* and *backward* scans as independent imaging channels. There are thus two images with Ku-H, two with Ku-V, etc... in
total eight imaging channels for each swath. When doing sea-ice motion tracking with CIMR, one can thus do an implicit merging with 16 pairs of images (fwd-fwd, fwd-bck, bck-fwd, and fwd-bck) for each of
Ku-V, Ku-H, Ka-V, and Ka-H. It is expected that using the forward and backward scans as part of the implicit merging will be benefitial to reduce the retrieval uncertainty and limit the number of
rogue vectors. This will have to be validated in the product development phase since CIMR is the first passive microwave mission offering full scans for sea-ice motion tracking.

## Swath-to-swath Motion Tracking (Level-2 strategy)

One of the key characteristics of the CIMR Level-2 sea-ice drift product is that it will be a "swath-to-swath" product, thus computed at the intersection of individual swaths. {numref}`fig_s2s` illustrates the concept.

```{figure} ./static_imgs/swath_to_swath.png
--- 
name: fig_s2s
width: 75%
---
(a) Daily average map of AMSR2 36.5 GHz V-pol TB on 1 December 2019 (greys) for the Arctic and two individual gridded swaths on the same day (blues: 01:16:34 UTC; reds: 19:24:55 UTC). The sea-ice region of overlap between the two swaths is highlighted in greens and is where S2S drift vectors can be computed. (b) Similar but for the Antarctic on 15 August 2019 (blues: 01:43:45 UTC; reds: 16:31:28 UTC). Reproduced from {cite:p}`lavergne:2021:s2s`.
```

The advantages of choosing a swath-to-swath approach (instead of a daily L3 product) are presented in {cite:t}`lavergne:2021:s2s`. We reproduce their conclusions below:

```{epigraph}
We investigate the feasibility and impact of adopting a swath-to-swath (S2S) vs. daily map (DM) framework for the processing of sea-ice motion from modern passive microwave satellite missions such as JAXA's AMSR2 in preparation for the future CIMR mission. We find that S2S sea-ice drift vectors obtained from AMSR2 imagery are more accurate than the corresponding DM vectors when compared to GPS trajectories from on-ice buoys in both the Arctic and Antarctic. An S2S configuration also results in many more drift vectors on a daily basis: the number varies with latitude and depends on the orbital and swath characteristics of the satellite mission. Since S2S drift vectors can be prepared for each new incoming swath, this configuration yields much better timeliness, which is beneficial for several operational applications such as support to navigation and short-term sea-ice forecasting. One potential limitation to the S2S configuration is that it is more sensitive to inaccurate geolocation, especially if the geolocation errors are systematic (e.g. a shift in the flight direction).

As far the CIMR mission is concerned, we recommend the adoption of an S2S configuration for the Level-2 sea-ice drift product in the operational ground segment. Considering the microwave frequency channels, target spatial resolution, swath width and geolocation accuracy specified for the CIMR imagery, its Level-2 sea-ice drift product will allow for unprecedented spatial resolution, coverage and accuracy for a microwave radiometer mission. Several other new characteristics of the CIMR mission (e.g. the relatively high spatial resolution at 6.9 and 10.8 GHz, the backward and forward scans) will also contribute to an enhanced sea-ice drift product, but this requires further research.
```



