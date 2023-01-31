#cd(@__DIR__)
#using Pkg
#Pkg.activate(".")
#Pkg.instantiate()

using Documenter, EasyModelAnalysis

# Make sure that plots don't throw a bunch of warnings / errors!
ENV["GKSwstype"] = "100"
using Plots

#using Pkg
#Pkg.add(url = "https://github.com/AlgebraicJulia/ASKEM-demos/", rev = "pas/hackathon",
#        subdir = "lib")

mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/require", "[tex]/mathtools"]),
                           :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                                        "packages" => [
                                            "base",
                                            "ams",
                                            "autoload",
                                            "mathtools",
                                            "require",
                                        ])))

makedocs(sitename = "DARPA-ASKEM Evalution",
         authors = "Chris Rackauckas",
         modules = Module[EasyModelAnalysis],
         clean = true, doctest = false,
         strict = [
             :doctest,
             :linkcheck,
             :parse_error,
             :example_block,
             # Other available options are
             # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
         ],
         format = Documenter.HTML(analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  mathengine = mathengine),
         pages = [
             "DARPA-ASKEM Evalution" => "index.md",
             #"Scenario1/Evaluation_Scenario_1.md",
             #"Scenario2/Evalution_Scenario_2.md",
             "Scenario3/Evalution_Scenario_3.md",
         ])

deploydocs(repo = "github.com/ChrisRackauckas/ASKEM_Evaluation_Staging.jl")
