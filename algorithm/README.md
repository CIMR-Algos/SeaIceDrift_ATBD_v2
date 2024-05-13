## Compilation of the code

The python modules required in the notebook can be installed by using the
supplied .yml file, as follows:

conda env create -f reqs_cimr-devalgo-sid.yml
conda activate cimr-devalgo-sid

The core of the sea ice drift code is in C, which is called via cython.
This requires that some compilation is done in advance of running the
notebook. A list of compilation commands is found in algorithm/src_sied
or reproduced here:

gcc -g -I /usr/include fmerrmsg.c -c
gcc -g -I /usr/include fmextra.c -c
gcc -g -I /usr/include icedrift_instruments.c -c
gcc -g -I /usr/include icedrift_model.c -c
gcc -g -I /usr/include icedrift_common.c -c
gcc -g -I /usr/include icedrift_solve_common.c -c
gcc -g -I /usr/include icedrift_solve_filter.c -c
gcc -g -I /usr/include icedrift_prepost.c -c
gcc -g -I /usr/include -I src_simplex icedrift_solve_core.c -c
python icedrift_cffi_build.py
