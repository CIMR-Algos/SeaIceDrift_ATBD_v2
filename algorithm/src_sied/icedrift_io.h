
#ifndef ICEDRIFT_IO_H
#define ICEDRIFT_IO_H

int load_nwp_field(char *filename,size_t *ne, size_t *nx, size_t *ny, float **DegHours);

int load_image_from_netcdf_new(char *filename, short nbWaveBands, char **WaveBands,
      size_t *ne, size_t *nx, size_t *ny, 
      float **lat, float **lon, short loadLatLon,
      float ***avgT,float ***obs, short ***flag,
      short **mask, float *fillvalue,
      float *map_Ax, float *map_Ay, float *map_Bx, float *map_By, char *map_projstr, fmsec1970 *time);

int load_image_from_netcdf(char *filename, short nbWaveBands, char **waveBands,
      size_t *ne, size_t *nx, size_t *ny, 
      float **lat, float **lon, short loadLatLon,
      float ***avgT,float ***obs,short ***flag,
      short **mask, float *fillvalue,
      float *map_Ax, float *map_Ay, float *map_Bx, float *map_By, char *map_projstr, fmsec1970 *time);

extern int loadOutPositions_ascii(int npoints, char posFile[],double *olat,double *olon);

#endif /* ICEDRIFT_IO_H */
