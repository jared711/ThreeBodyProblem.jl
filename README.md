![ThreeBodyProblem.jl](docs/src/assets/banner.png)

<!--- [![Build Status](https://travis-ci.com/jared711/ThreeBodyProblem.jl.svg?branch=master)](https://travis-ci.com/jared711/ThreeBodyProblem.jl) -->
<!---[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jared711/ThreeBodyProblem.jl/master) -->
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jared711/ThreeBodyProblem.jl?svg=true)](https://ci.appveyor.com/project/jared711/ThreeBodyProblem-jl)
[![Codecov](https://codecov.io/gh/jared711/ThreeBodyProblem.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jared711/ThreeBodyProblem.jl)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://jared711.github.io/ThreeBodyProblem.jl/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jared711/ThreeBodyProblemExamples.jl/main)

An astrodynamics package for working in the three-body problem.

## Requirements
Julia 1.6+ https://julialang.org/downloads/

## Installation
ThreebodyProblem.jl is on the General Registry. That means you can download it and start using it quickly in the Julia repl
```julia
Pkg.add("ThreeBodyProblem")
using ThreeBodyProblem
```

## Examples
We've put together some working examples that demonstrate many of the functionalities of ThreeBodyProblem.jl. They can be found in the [ThreeBodyProblemExamples.jl](https://github.com/jared711/ThreeBodyProblemExamples.jl) repository. Each example has a .jl and .ipynb file. The .jl files can be run line by line in your favorite IDE. Just download the ThreeBodyProblemExamples.jl repo and start a Julia repl.
```bash
git clone https://github.com/jared711/ThreeBodyProblemExamples.jl.git
cd ThreeBodyProblemExamples.jl
julia
```
From the Julia repl run
```julia
Pkg.activate(".")
Pkg.instantiate()
```
to download all the necessary packages.

### Jupyter Notebooks

To run Jupyter notebooks on your local machine, you need to have the IJulia package downloaded, which it should be if you followed the steps above. Then run the notebook() command from your desired working directory
```julia
using IJulia
notebook(dir=".",detached=true)
```

VS Code has a [Jupyter notebook extension](https://code.visualstudio.com/docs/datascience/jupyter-notebooks), as do some other popular IDEs.

You can also run the Jupyter notebooks on your browser through binder. Just click the binder icon [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jared711/ThreeBodyProblem.jl/master) and wait for a JupyterLab window to appear. It may take a few minutes to load.  Pick the Jupyter notebook (filetype .ipynb) you want to run from the file tab on the left. Pressing CTRL+ENTER will run a block of code, while SHIFT+ENTER will run a block and move to the next one.
<!---Eventually this window should appear. ![image](https://user-images.githubusercontent.com/25643720/216104189-4d60e01b-dc72-4946-b72f-0d774bd78187.png)
Navigate to the example folder on the left side of the screen,-->

## Documentation
The documentation is found [here](https://jared711.github.io/ThreeBodyProblem.jl).
