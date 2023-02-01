# DARPA ASKEM 6 Month Evaluations

## Reproducing the Results

This system is a fully reproducer builder that generates the results from scratch. It downloads Julia,
downloads the packages, builds the results, and then hosts them at the website 
https://chrisrackauckas.github.io/ASKEM_Evaluation_Staging/dev/.

### Running the Builder

The simplest way to reproduce the results is to run the builder by opening a PR. The `Documentation` job is the
job generating the results, and if successful (green) then its artifacts are pushed to the repository on merge.

### Running Locally via the `make.jl` Script

To run locally, start by downloading the Github repository and set the current directory to 
`ASKEM_Evaluation_Staging/docs`. Then open Julia and run:

```julia
cd("to the ASKEM_Evaluation_Staging/docs directory")

# Setup packages
using Pkg
Pkg.activate(".") 
Pkg.instantiate()

# Run the build script
include("make.jl")
```

The outputs will be generated into the `ASKEM_Evaluation_Staging/docs/build` directory and can be opened in 
a web browser.

### Running Individual Files

The individual files can be built as well by running the code in the Markdown files. We recommend grabbing the
packages using the Project.toml via:

```julia
cd("to the ASKEM_Evaluation_Staging/docs directory")

# Setup packages
using Pkg
Pkg.activate(".") 
Pkg.instantiate()
```

Then open up the .md and run the code! If using the [Julia VS Code plugin](https://code.visualstudio.com/docs/languages/julia)
then inline evaluation will directly work.

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to SciML.

  - See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
  - There are a few community forums:
    
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Slack](https://julialang.org/slack/)
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
      + On the [Julia Discourse forums](https://discourse.julialang.org)
      + See also [SciML Community page](https://sciml.ai/community/)

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@raw html
You can also download the
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
       "/assets/Manifest.toml"
```

```@raw html
">manifest</a> file and the
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
       "/assets/Project.toml"
```

```@raw html
">project</a> file.
```
