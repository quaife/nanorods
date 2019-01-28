# Compilation Instructions

1. `cd src/iagm/iagm`
2. `make`
3. `cd ../..`
4. `make`
5. `cd ../laplace_fmm`
6. `make`
7. `cd ../stokes_fmm`
8. `make`
9. `cd ..`
10. `LD_PRELOAD="`*abolute path to src/iagm/iagm/lib/libiagm.so*`" matlab &`
11. run `createpaths.m`
12. run `extensional.m`
