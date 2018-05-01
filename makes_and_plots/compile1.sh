rm a.out *.mod 2>/dev/null
cp -r /home/elena/Desktop/thesis/2d-functions/SA_code/SA_code/*.f95 ./
gfortran -Wall -Wextra -Wconversion -pedantic -g3 -fbacktrace -fbounds-check -O0 -Wunused-parameter main.f95 spsamod.f95 fdsamod.f95 othermod.f95
