# Abstract

Pushed by winds and ocean currents, sea ice is always on the move. Satellite remote sensing is an effective tool to monitor this motion, especially from microwave imagery. 

Monitoring of sea ice drift vectors by means of passive microwave is commonly based on motion tracking algorithms working from pairs of brightness temperature imagery in both horizontal and vertical polarisation.
Depending on the satellite missions, different microwave frequencies have been used, but one generally expects higher accuracy from the higher resolution channels.
Accordingly, this version of the Sea Ice Drift Level-2 product focuses on the Ku and Ka imagery.

Contrarily to other existing operational sea-ice drift product (such as those of IFREMER Cersat or EUMETSAT OSI SAF), the CIMR sea-ice drift product is a Level-2 product. It is thus computed from pairs of
overlapping swaths, not from pairs of daily averaged maps. This "swath-to-swath" setup was introduced in {cite:t}`lavergne:2021:s2s`.

Key assets of CIMR in terms of Sea Ice Drift monitoring are 1) its increased spatial resolution at Ku- and Ka-band (more accurate drift vectors), 2) larger swath (more drift vectors) and 3) the forward
and backward scans (better quality control). The main source of uncertainty for the CIMR Sea Ice Drift product will probably be the geolocation uncertainty, hence the importance of the pointing
accuracy (in orbit) and geolocation correction steps (in Level-1).

Sea-ice drift vectors are assessed against trajectories from in situ buoy deployed on the ice. A sufficient number of these buoys must be ensured throughout the lifetime of the mission, both in the northern and
southern hemisphere. Sea-ice drift vectors from Synthetic Aperture Radar (SAR) as derived in the Copernicus Marine Service Sea Ice Thematic Assembly Center from Sentinel-1 (in the future S1-NG and ROSE-L) can also
be used for validation.

