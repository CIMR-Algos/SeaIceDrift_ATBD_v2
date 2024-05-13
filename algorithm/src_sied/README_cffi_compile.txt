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
