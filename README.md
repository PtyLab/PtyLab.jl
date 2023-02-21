# PtyLab.jl

| **Build Status**                          | **Code Coverage**               |
|:-----------------------------------------:|:-------------------------------:|
| [![][CI-img]][CI-url] | [![][codecov-img]][codecov-url] |

# License
Please read the [license](https://github.com/PtyLab/PtyLab.jl/blob/main/LICENSE.md) before usage!


# PtyLab.jl

Julia flavored implementation of PtyLab. PtyLab is a software for Ptychography reconstruction.
Please [cite this arXiv for the moment](https://arxiv.org/abs/2301.06595).

## Installation
Via the Julia-REPL:
```julia
julia> ] add https://github.com/PtyLab/PtyLab.jl/
```

## Feature Set
Right now, only CP reconstruction works. For more functionality, please see the [PtyLab.m](https://github.com/PtyLab/PtyLab.m) and the [PtyLab.py](https://github.com/PtyLab/PtyLab.py).

The following features are implemented:
* simple Ptychography reconstruction (loading from a .hdf5 file)
* ePIE
* CUDA support
* probe center of mass restriction
* regular randomized grid generation (for usage of a translation stages.

The structure should be flexible to add more solvers, etc.


## Tutorials
* [Real data reconstruction](examples/Introduction_real_data_reconstruction.ipynb)
* [Simple simulation and reconstruction](examples/simple_simulation_and_reconstruction.ipynb)
* [Another artificial simulation and reconstruction](examples/cuda_simulation_and_reconstruction.ipynb)


## JuliaCon 2022
* [Pluto notebook](examples/JuliaCon_2022.jl) for this video

<a  href="https://www.youtube.com/watch?v=pDp83OxBJ_I"><img src="docs/src/assets/juliacon.png"  width="300"></a>



[docs-dev-img]: https://img.shields.io/badge/docs-dev-pink.svg
[docs-dev-url]: https://ptylab.github.io/PtyLab.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-darkgreen.svg
[docs-stable-url]:  https://ptylab.github.io/PtyLab.jl/stable/

[CI-img]: https://github.com/ptylab/PtyLab.jl/actions/workflows/ci.yml/badge.svg
[CI-url]: https://github.com/ptylab/PtyLab.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/PtyLab/PtyLab.jl/branch/main/graph/badge.svg?token=OQ6BQCUFQB
[codecov-url]: https://codecov.io/gh/ptylab/PtyLab.jl
