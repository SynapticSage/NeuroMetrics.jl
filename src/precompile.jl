# Uses https://github.com/JuliaLang/PackageCompiler.jl to speed up use of my library
# Spends less time in precompile when a function/module is first called

using PackageCompiler
function precompile
end
