[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jared711/ThreeBodyProblem.jl/master)
<!-- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jared711/ThreeBodyProblemExamples.jl/main) -->


# Examples
This section of the repository contains helpful examples, both in Julia scripts and Jupyter notebooks. 
If you are brand new to the package, go through the examples in order to get a feel for how everything works.
However, there are still elements of the package that are not included in the examples.

Each example has a .jl file and a .ipynb file. The .jl file can be run line by line in your favorite IDE while the .ipynb file is a Jupyter notebook.

## Jupyter Notebooks
To run Jupyter notebooks on your local machine, you need to have the IJulia package downloaded.
```julia
Pkg.add("IJulia")
using IJulia
notebook(dir=".",detached=true)
```

You can also run the Jupyter notebooks on your browser through binder. Just click the binder icon [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jared711/ThreeBodyProblem.jl/master) and wait for a JupyterLab window to appear. It may take a few minutes to load. Eventually this window should appear.
![image](https://user-images.githubusercontent.com/25643720/216104189-4d60e01b-dc72-4946-b72f-0d774bd78187.png)
Pick the Jupyter notebook (filetype .ipynb) you want to run. Pressing CTRL+ENTER will run a block of code, while SHIFT+ENTER will run a block and move to the next one.

Check out the documentation or the source code itself for more tips.
