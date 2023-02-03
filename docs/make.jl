#cd(@__DIR__)
#using Pkg
#Pkg.activate(".")
#Pkg.instantiate()

using Documenter, EasyModelAnalysis

# Make sure that plots don't throw a bunch of warnings / errors!
ENV["GKSwstype"] = "100"
using Plots

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/require", "[tex]/mathtools"]),
                           :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                                        "packages" => [
                                            "base",
                                            "ams",
                                            "autoload",
                                            "mathtools",
                                            "require",
                                        ])))

makedocs(sitename = "DARPA-ASKEM Evaluation",
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
             "DARPA-ASKEM Evaluation" => "index.md",
             "Scenario1/Evaluation_Scenario_1.md",
             "Scenario2/Evaluation_Scenario_2.md",
             "Scenario3_prep/scenario3.md",
             "Scenario3/Evaluation_Scenario_3.md",
         ])

deploydocs(repo = "github.com/ChrisRackauckas/ASKEM_Evaluation_Staging")
