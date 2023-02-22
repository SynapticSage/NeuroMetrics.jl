using GLM, Lasso, Distributions, ThreadSafeDicts
@time include(scriptsdir("isolated","load_isolated.jl"))
# @time include("./load_isolated.jl")
fn = path_iso(opt; append="_cyclewise")
@load "/home/ryoung/Projects/goal-code/data/isolated/iso_animal=RY16_day=36_tet=ca1ref.jld2"
@load "/home/ryoung/Projects/goal-code/data/isolated/iso_animal=RY16_day=36_tet=ca1ref_cyclewise.jld2"

# Clean data frame
col_all_zero = map(v->all(v.==0), eachcol(df))
df = df[!, Not(names(df)[col_all_zero])]


function get_dx_dy(df, relcyc)
    dx = @subset(df, :relcycs .== relcyc)
    dy = @subset(df, :relcycs .== 0)
    dx = groupby(dx, [:cyc_batch, :cyc_match])
    dy = groupby(dy, [:cyc_batch, :cyc_match])
    kx, ky = keys(dx.keymap), keys(dy.keymap)
    k = intersect(kx,ky)
    dx = sort(vcat([dx[kk] for kk in k]...), [:cyc_batch, :cyc_match])
    dy = sort(vcat([dy[kk] for kk in k]...), [:cyc_batch, :cyc_match])
    dx, dy
end
@time dx_dy = Dict(relcyc => get_dx_dy(df, relcyc) for relcyc in -8:8)

