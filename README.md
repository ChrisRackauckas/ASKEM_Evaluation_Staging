# ASKEM_Evaluation_Staging

See the results at https://chrisrackauckas.github.io/ASKEM_Evaluation_Staging/

## Getting Started with Running the Code

To run locally, start by downloading this Github repository. Then open Julia and run:

```julia
cd("to the ASKEM_Evaluation_Staging/docs directory")

# Setup packages
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Run the build script
include("make.jl")
```

## Dependencies Information

All information about the required dependencies is stored in the Project.toml file:

https://github.com/ChrisRackauckas/ASKEM_Evaluation_Staging/blob/main/docs/Project.toml

## Using the Pluto Notebooks

The Pluto notebooks are fully relocatable and download the dependencies and data on-demand, making them easy
install and use. To do this, first install Julia v1.8.5 from https://julialang.org/downloads/. Then
open up Julia and run:

```julia
using Pkg
Pkg.add("Pluto")
import Pluto
Pluto.run()
```

When this window opens in your browser, navigate to the Pluto notebook. Opening it will start the download
of all data and dependencies, along with the automatic installation and compilation of all dependencies.
This may take awhile due to the complexity of the downloads.

## Example of the Output Plots

![](https://user-images.githubusercontent.com/1814174/216075828-e7e7289b-1300-483f-bc3d-9bf73820fc33.png)
