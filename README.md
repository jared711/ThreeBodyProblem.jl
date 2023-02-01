![ThreeBodyProblem.jl](docs/src/assets/banner.png)

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://jared711.github.io/ThreeBodyProblem.jl/)
[![Build Status](https://travis-ci.com/jared711/ThreeBodyProblem.jl.svg?branch=master)](https://travis-ci.com/jared711/ThreeBodyProblem.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jared711/ThreeBodyProblem.jl?svg=true)](https://ci.appveyor.com/project/jared711/ThreeBodyProblem-jl)
[![Codecov](https://codecov.io/gh/jared711/ThreeBodyProblem.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jared711/ThreeBodyProblem.jl)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jared711/ThreeBodyProblem.jl/master)

An astrodynamics package for working in the three body problem.

## Requirements
Julia 1.6+ https://julialang.org/downloads/

## Installation
ThreebodyProblem.jl is not on the General Registry yet, so you must use Git to clone it to your computer.
```shell
git clone https://github.com/jared711/ThreeBodyProblem.jl.git
```
Then, you can add it to your Julia package environment by running the following from a Julia REPL.
```julia
using Pkg
Pkg.add(path="[local path]/ThreeBodyProblem.jl")
```
Note that this syntax may vary slightly on Windows.

## Examples
We've put together some working examples that demonstrate many of the functionalities of ThreeBodyProblem.jl. They can be found in the example folder [here](https://github.com/jared711/ThreeBodyProblem.jl/tree/master/example). Each example has a .jl file and a .ipynb file. The .jl file can be run line by line in your favorite IDE while the .ipynb file is a Jupyter notebook. The .ipynb 

### Jupyter Notebooks
To run Jupyter notebooks on your local machine, you need to have the IJulia package downloaded.
```julia
using Pkg
Pkg.add("IJulia")
```
Then run the following from your desired working directory
```julia
using IJulia
notebook(dir=".",detached=true)
```

You can also run the Jupyter notebooks on your browser through binder. Just click the binder icon [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jared711/ThreeBodyProblem.jl/master) and wait for a JupyterLab window to appear. It may take a few minutes to load. Eventually this window should appear.
![image](https://user-images.githubusercontent.com/25643720/216104189-4d60e01b-dc72-4946-b72f-0d774bd78187.png)
Navigate to the example folder on the left side of the screen, then pick the Jupyter notebook (filetype .ipynb) you want to run. Pressing CTRL+ENTER will run a block of code, while SHIFT+ENTER will run a block and move to the next one.

## Documentation
Documentation is found [here](https://jared711.github.io/ThreeBodyProblem.jl)
