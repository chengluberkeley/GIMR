# A faster algorithm solving a generalization of isotonic median regression and a class of fused lasso problems

This repo implements the algorithm in paper [**Hochbaum, D.S. and Lu, C.** A faster algorithm solving a generalization of isotonic median regression and a class of fused lasso problems. *SIAM J. Opt.*, 27(4): 2563-2596, 2017.](https://hochbaum.ieor.berkeley.edu/html/pub/Isotonic-HLuSIAM-Opt2017.pdf) Algorithms with and without [Dynamic Path](https://github.com/chengluberkeley/DynamicPath) are implemented.

## Run prebuilt executable
For a sample run of the implemented algorithms:
```
cd bin
./test_main
```
Results are output to console and written to `/tmp/` folder.

## Build from the source
This project is a `cmake` project. To build from the source:
```
mkdir build && cd build
cmake .. && make -j5
```

## Reference

Please cite the paper if you use the algorithms.

**Hochbaum, D.S. and Lu, C.** A faster algorithm solving a generalization of isotonic median regression and a class of fused lasso problems. *SIAM J. Opt.*, 27(4): 2563-2596, 2017.
