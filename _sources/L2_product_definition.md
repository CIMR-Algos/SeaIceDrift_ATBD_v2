# Level-2 product definition

{numref}`l2_variables`, {numref}`l2_attributes`, {numref}`l2_global_attributes` define the content (TBC) of the Level-2 sea-ice drift product files. Their structure
is inspired by that of the EUMETSAT OSI SAF products.

```{table} NetCDF Group: Processed data (TBC)
:name: l2_variables
|  name                  | description | units | dimensions |
|  --------------------- | ----- | ---- | ---- |
|  crs           | Definition of the EASE2 product grid | n/a | n/a |
|  dX            | Component of the drift along the X-axis of the EASE2 grid | m | (nx,ny) |
|  dY            | Component of the drift along the Y-axis of the EASE2 grid | m | (nx,ny) |
|  lat_1         | Latitude of the end point of the drift vector | degrees North | (nx,ny) |
|  lon_1         | Longitude of the end point of the drift vector | degrees East | (nx,ny) |
|  sX            | Uncertainty (one standard deviation) of the dX component | m | (nx,ny) |
|  sY            | Uncertainty (one standard deviation) of the dY component | m | (nx,ny) |
|  cXY           | Correlation between the uncertainty in dX and in dY | m | (nx,ny) |
|  status_flag   | A flag indicating status of retrieval, e.g. "nominal", "over land", "close to ice edge", "over ocean", "corrected from rogue vector",\... | n/a | (nx,ny) |
```

```{note}
The {term}`CIMR` Level-2 {term}`Sea Ice Drift` product is an Earth-gridded product. The dimensions of each variable in the Level-2 file
refer to the (nx,ny) dimensions of the product grid in an {term}`EASE2` polar projection (see {numref}`grids`).
```

Each data variable in the `Processed data` group holds conventional attributes following CF-1.7 or above:

```{table} Standard variable attributes (TBC)
:name: l2_attributes
|  name                  | description |
|  --------------------- | ----- | 
| standard_name | A standard name that references a description of a variable\"s content                 |
| long_name     | A descriptive name that indicates a variable\"s content                                |
| \_FillValue   | value that will appear in the dataset when the data is either missing or undefined     |
| units         | unit of measure                                     |
```


```{table} Some global attributes in the Level-2 product files (TBC)
:name: l2_global_attributes
|  name                  | description |
|  --------------------- | ----- | 
|  title                 | CIMR L2 Sea Ice Drift product |
|  processing_level      | Level-2 |
|  time_coverage_start   | valid *start* time of the product |
|  time_coverage_end     | valid *end* time of the product |
|  area                  | (Northern or Southern) Hemisphere |
|  TBD                   | TBD |
```

