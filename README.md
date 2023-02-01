# PtyLab.jl


| **Documentation**                       | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][CI-img]][CI-url] | [![][codecov-img]][codecov-url] |

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

The following features are implemented:
* simple Ptychography reconstruction (loading from a .hdf5 file)
* ePIE
* CUDA support
* probe center of mass restriction
* regular randomized grid generation

The structure should be flexible to add more solvers, etc.


## JuliaCon 2022
<a  href="https://www.youtube.com/watch?v=pDp83OxBJ_I"><img src="docs/src/assets/juliacon.png"  width="300"></a>



[docs-dev-img]: https://img.shields.io/badge/docs-dev-pink.svg
[docs-dev-url]: https://ptylab.github.io/PtyLab.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-darkgreen.svg
[docs-stable-url]:  https://ptylab.github.io/PtyLab.jl/stable/

[CI-img]: https://github.com/ptylab/PtyLab.jl/actions/workflows/ci.yml/badge.svg
[CI-url]: https://github.com/ptylab/PtyLab.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/PtyLab/PtyLab.jl/branch/main/graph/badge.svg?token=OQ6BQCUFQB
[codecov-url]: https://codecov.io/gh/ptylab/PtyLab.jl
