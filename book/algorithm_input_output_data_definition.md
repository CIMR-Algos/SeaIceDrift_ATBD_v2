# Algorithm Input and Output Data Definition (IODD)

Except for the input L1B TB data, the auxiliary and output data are all on {term}`EASE2` grids. Refer to {numref}`grids`.

## Input data

| Field | Description | Shape/Amount |
| ---   | ----------- | ------------ |
| L1B TB | L1B Brightness Temperature at K and KA-bands (both H and V polarization) | full swath or section of it (Nscans, Npos) |
| L1B NeΔT | Random radiometric uncertainty of the channels | full swath or section of it (Nscans, Npos) |

## Output data

| Field | Description | Shape/Amount |
| ---   | ----------- | ------------ |
| sea-ice drift vectors | dX and dY components of the motion vectors. | (nx,ny) on the `(n,s)h_ease2-250` grid. |
| sea-ice drift vectors uncertainty | dX and dY components of the uncertainty vector (1-sigma). | (nx,ny) on the `(n,s)h_ease2-250` grid. |
| status flags | indicate the reasons for missing, bad, or nominal vectors. | (nx,ny) on the `(n,s)h_ease2-250` grid. |

## Auxiliary data

| Field | Description | Shape/Amount |
| ---   | ----------- | ------------ |
| sea-ice and land mask | A recent sea-ice concentration field including land information | SIC and land-mask on the `(n,s)h_ease2-005` grid. |


