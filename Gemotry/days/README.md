# A Massive Fractal In Days, Not Years

This code complements the Journal of Computer Graphics Tools submission, *A Massive Fractal In Days, Not Years*, by Kim and Duff.

## Building the Code
### Mac OS

A `Makefile.osx` is provided for building on the Mac:

```
make -f Makefile.osx
```

You will probably have to point `make` to your CUDA installation by modifying the `CUDA_PATH` parameter on its first line.

### Linux

The default `Makefile` builds on Linux, so just call `make`. You will (still) probably need to point it to your CUDA installation by modifying `CUDA_PATH` on the first line of `Makefile`.

## Running the Code

Two scripts are provided for running the code.

If you execute `run`, it will generate a relatively small, 500^3 bunny.

If you try `run_massive`, it will try to compute the 11,500^3 bunny that previously took 17 days to generate. It should output checkpoint files, and if you kill the computation after a checkpoint, it should start from there the next time `run_massive` is called.

****

*Theodore Kim, 9/6/19*
