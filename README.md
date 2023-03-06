# PtyLab.jl

Conventional Ptychography is a lensless microscopy imaging technique which captures a sequence of light diffraction patterns to solve the optical phase problem. The resulting datasets are large and can typically not directly be solved. Instead, iterative reconstruction algorithms with low runtime memory footprint are employed. Here we present PtyLab.jl, a software for ptychographic data analysis and demonstrate how a functional programming style in Julia allows for performant iterative algorithms.


Please [see this arXiv article](https://arxiv.org/abs/2301.06595).


| **Build Status**                          | **Code Coverage**               |
|:-----------------------------------------:|:-------------------------------:|
| [![][CI-img]][CI-url] | [![][codecov-img]][codecov-url] |

# License
Please read the [license](https://github.com/PtyLab/PtyLab.jl/blob/main/LICENSE.md) before usage! 
One important part is the use of this software "for academic, non-commercial purposes" only.


## Installation
Via the [Julia-REPL](https://julialang.org/) add [PtyLab.jl](https://github.com/PtyLab/PtyLab.jl) with the following command:
```julia
julia> ] add https://github.com/PtyLab/PtyLab.jl/
```

To start the Jupyter notebooks, clone or download a [.zip](https://github.com/PtyLab/PtyLab.jl/archive/refs/heads/main.zip) of this repository. 
Extract it.
Then open a Julia REPL and type the following to install the Jupyter kernel:
```julia
julia> ] add IJulia

julia> using IJulia

julia> notebook()
```
A browser should open. Navigate to the location of the notebooks and open them.

## Feature Set
Right now, only CP reconstruction works. For more functionality, please see the [PtyLab.m](https://github.com/PtyLab/PtyLab.m) and the [PtyLab.py](https://github.com/PtyLab/PtyLab.py).

The following features are implemented:
* simple Ptychography reconstruction (loading from a .hdf5 file)
* ePIE
* CUDA support
* probe center of mass restriction
* regular randomized grid generation (for usage of translation stages)

The structure should be flexible to add more solvers, etc.


## Tutorials
* [Real data reconstruction](examples/Introduction_real_data_reconstruction.ipynb)
* [Simple simulation and reconstruction](examples/simple_simulation_and_reconstruction.ipynb)
* [Another artificial simulation and reconstruction](examples/cuda_simulation_and_reconstruction.ipynb)

## Questions?
Feel free to open an issue if anything does not work or is unclear!
You can also join [my conference room](https://epfl.zoom.us/j/9384091974). Give me a minute to join!

## JuliaCon 2022
* [Pluto notebook](examples/JuliaCon_2022.jl) for this video

<a  href="https://www.youtube.com/watch?v=pDp83OxBJ_I"><img src="docs/src/assets/juliacon.png"  width="300"></a>

## Citing
```
@misc{https://doi.org/10.48550/arxiv.2301.06595,
  doi = {10.48550/ARXIV.2301.06595},
  url = {https://arxiv.org/abs/2301.06595},
  author = {Loetgering, Lars and Du, Mengqi and Flaes, Dirk Boonzajer and Aidukas, Tomas and Wechsler, Felix and Molina, Daniel S. Penagos and Rose, Max and Pelekanidis, Antonios and Eschen, Wilhelm and Hess, Jürgen and Wilhein, Thomas and Heintzmann, Rainer and Rothhardt, Jan and Witte, Stefan},
  keywords = {Computational Physics (physics.comp-ph), Image and Video Processing (eess.IV), Optics (physics.optics), FOS: Physical sciences, FOS: Physical sciences, FOS: Electrical engineering, electronic engineering, information engineering, FOS: Electrical engineering, electronic engineering, information engineering},
  title = {PtyLab.m/py/jl: a cross-platform, open-source inverse modeling toolbox for conventional and Fourier ptychography},
  publisher = {arXiv},
  year = {2023},
  copyright = {Creative Commons Attribution 4.0 International}
}
```

[docs-dev-img]: https://img.shields.io/badge/docs-dev-pink.svg
[docs-dev-url]: https://ptylab.github.io/PtyLab.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-darkgreen.svg
[docs-stable-url]:  https://ptylab.github.io/PtyLab.jl/stable/

[CI-img]: https://github.com/ptylab/PtyLab.jl/actions/workflows/ci.yml/badge.svg
[CI-url]: https://github.com/ptylab/PtyLab.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/PtyLab/PtyLab.jl/branch/main/graph/badge.svg?token=OQ6BQCUFQB
[codecov-url]: https://codecov.io/gh/ptylab/PtyLab.jl
