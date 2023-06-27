# Introduction, purpose and scope

The balance between air drag, ocean drag and lateral forces controls the motion of sea ice {cite:p}`lepparanta:2011:seaicedrift`.
At the local scale, sea-ice motion can both be a facilitator and impediment to ship navigation, opening and closing routes, opening leads, or forming pressure ridges.
At the larger regional to basin scales, sea-ice motion (a.k.a. sea-ice drift) exports sea ice to lower latitudes where it melts, contributing to the redistribution of fresh water {cite:p}`haine:2015:freshwater`.
Inside the Arctic Ocean, re-circulation of sea-ice (e.g. in the Beaufort Sea) leads to the ageing and thickening of the ice pack towards the northern coasts of the Canadian Arctic Archipelago and Greenland {cite:p}`timmermans:2020:arctic`.
Sea-ice drift also plays a role in sea-ice formation and ocean circulation via the formation of coastal latent heat polynyas {cite:p}`ohshima:2016:polynyas`, as well as in the transport of sediments and other tracers across ocean basins
{cite:p}`krumpen:2019:transpolar`. With climate change, the area and thickness of sea ice is reduced in the Arctic, which leads to a more mobile sea-ice cover and positive trends in sea-ice velocity {cite:p}`spreen:2011:drifttrend,kwok:2013:arctictrend`.
Trends in sea-ice motion, linked to trends in wind speed, are also observed in the Southern Hemisphere (SH) {cite:p}`holland:2012:antarctictrend,kwok:2017:antarctictrend`.

Satellite remote sensing has developed as an attractive option to monitor sea-ice drift consistently across the polar sea-ice cover at a daily to sub-daily frequency.
The initial work by {cite:t}`ninnis:1986:avhrr-drift` was followed by many investigators using a variety of satellite imaging sensor technologies as input, including visible and infrared radiometry {cite:p}`emery:1991:avhrr-drift`,
microwave radiometry and scatterometry {cite:p}`agnew:1997:drift,kwok:1998:ssmi85-drift,liu:1999:drift-wavelet,lavergne:2010:cmcc-jgr,girard-ardhuin:2012:drift`,
and synthetic aperture radar (SAR) {cite:p}`kwok:1990:rgps,komarov:2014:sar-drift,muckenhuber:2016:sar.drift`. The various imaging technologies however lead to sea-ice motion fields with different characteristics, e.g.
medium spatial resolution (∼ 20 km) and coverage limited by cloud cover for the visible and infrared radiometry, high spatial resolution (∼ 5–10 km) but coverage limited by acquisition repeat cycles for the SAR imagery,
and coarse spatial resolution (> 30 km) and daily complete coverage for the microwave radiometers and scatterometers. Despite the imaging technologies being very different from each other, the motion tracking algorithms
employed are quite similar and stem from the maximum cross-correlation (MCC) technique {cite:p}`emery:1986:sst-drift`.

```{note}
The characteristics of the CIMR mission, and especially the spatial resolution obtained at Ku (5 km)
and Ka (4 km) should allow global, year-round (irrespective of solar illumination and cloud cover), sub-daily moniroting at medium spatial resolution (25 km), which is unprecedented.
```

The Level-2 sea-ice drift product described in this ATBD has been specifically designed to fully exploit the CIMR mission, in particular:
1. the CIMR sea-ice drift product is a Level-2 product, adopting a swath-to-swath approach {cite:p}`lavergne:2021:s2s`. This is conversely to all existing sea-ice drift products from passive microwave missions (OSI SAF, IFREMER, NSIDC, etc...) that are Level-3 products (available once a day from daily averaged maps of brightness temperatures). SAR-based products (e.g. from the CMEMS Sea Ice Thematic Assembly Center) are generally swath-to-swath products (but with spatial coverage limited by the coverage of the SAR missions).
2. the CIMR sea-ice drift product aims at a 25 km grid spacing, taking full advantage of the Ku and Ka imagery.
3. sea-ice drift is derived both from the forward and backward scans separately to improve quality control (filtering of "rogue" vectors)
 

