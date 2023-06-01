# This script just explores whether my underling data ought to be predictive
# if its changing enough over the 2 hours of recording that prediction is
# ruined

include("../imports_isolated.jl")
include("./pfc_cyclewise.jl")

# Make some helper functions!
    """
    string representations of columsn denoting iso spiking cells
    """
    is(df) = replace.(uicellcols(df), ["_i"=>""])
    """
    string representations of columsn denoting spiking cells
    """
    us(df) = ucellcols(df)
    """
    int rep of iso spiking cell columns
    """
    ic(df) = parse.(Int,is(df))
    """
    int rep of spiking cell columns
    """
    uc(df) = parse.(Int,us(df))

    """
    subset of cells dataframe by relevant columns
    """
    cu(cells::DataFrame, df::DataFrame)       = @subset(cells, :unit .∈ (uc(df),))
    icu(cells::DataFrame, df::DataFrame)      = @subset(cells, :unit .∈ (ic(df),))
    cu(cells::Vector{String}, df::DataFrame)  = cu(parse.(Int,cells), df)
    icu(cells::Vector{String}, df::DataFrame) = icu(parse.(Int,cells), df)
    cu(cells::Vector{Int}, df::DataFrame)     = cells .∈ (ic(df),)
    icu(cells::Vector{Int}, df::DataFrame)    = cells .∈ (uc(df),)
    function Base.:*(x::Int64, s::String) 
        string(x) * s
    end

iset = is(df)
uset = us(df)
icells = ic(df)
ucells = uc(df)
CU  = cu(cells, df)
ICU = icu(cells, df)

@assert CU.unit == parse.(Int,ucellcols(df)) "Cells not in same order as cell table"
@assert ICU.unit == parse.(Int,iset) "Cells not in same order as cell table"


# Should we be able to predict between the cells?
area_sorted_cells = string.(sort(CU, :area).unit)
pfc_cells = string.(@subset(CU, :area .== "PFC").unit)
ca1_cells = string.(@subset(ICU, :area .== "CA1").unit)
pfc = area_sorted_cells .∈ (pfc_cells,)
ca1 = area_sorted_cells .∈ (ca1_cells,)

function blank_lowerdiagonal!(C::Matrix)
    for i in 1:size(C,1), j in 1:i
        C[i,j] = NaN
    end
    nothing
end

function get_correlation_matrices(df)
    C = cor(Matrix(df[!,area_sorted_cells]))
    blank_lowerdiagonal!(C)
    Cca1    = C[ca1, ca1]
    Cpfc    = C[pfc, pfc]
    Cca1pfc = C[ca1, pfc]
    (;C, Cca1, Cpfc, Cca1pfc)
end

function plot_correlation_matrices(;Cca1, Cpfc, Cca1pfc, kws...)

    g= grid(2,2)
    ca1pct = sum(ca1)/length(ca1) * Plots.pct
    pfcpct = sum(pfc)/length(pfc) * Plots.pct
    g.widths  = [ca1pct, pfcpct]
    g.heights = [ca1pct, pfcpct]
    m = 0.2
    # Get yticks for Cca1 (cell labels)
    ca1ticks = string.(sort(CU[ca1,:], :area).unit)
    # Set yticks to be equal (locations, labels, size)
    ca1ticks = (collect(1:size(Cca1,1)), ca1ticks)
    ffont = Plots.font(3, "Courier")
    pfcticks = string.(sort(CU[pfc,:], :area).unit)
    pfcticks = (collect(1:size(Cpfc,2)), pfcticks)
    hCA1=heatmap(collect(1:size(Cca1,1)), collect(1:size(Cca1,2)), Cca1; 
                xticks=ca1ticks, yticks=ca1ticks, 
                    xtickfont=ffont, ytickfont=ffont,
                    colorbar_title="CA1",     
                    c=:PiYG_5,    clim=(-m, m),  colorbar=:left)
    hCA1PFC=heatmap(1:size(Cca1pfc,2), 1:size(Cca1pfc,1), Cca1pfc,
    colorbar_tickfontsize=3, colorbar_title="CA1-PFC", c=:vik,       clim=(-m,
        m), xticks=pfcticks, yticks=ca1ticks, xtickfont=ffont, ytickfont=ffont,
    colorbar=:right)

    hPFC = heatmap(collect(1:size(Cpfc,1)), collect(1:size(Cpfc,2)), Cpfc,
    colorbar_tickfontsize=3, colorbar_title="PFC",     c=:RdYlGn_11, clim=(-m,
        m))
    plot(hCA1,
        hCA1PFC,
        Plot.create_blank_plot(),
        hPFC, layout=g,
            xticks=pfcticks, yticks=pfcticks, xtickfont=ffont, ytickfont=ffont,
    )
end

overall = get_correlation_matrices(df)
plot_correlation_matrices(;overall...)

function checkbins(var::Symbol)
    println("Unique $var bins: ", length(unique(getproperty(Main,var).bins)))
end

# Prepare for a 10 minute bucketed approach
dT = diff(beh.time)
bins = Int(maximum(floor.(cumsum(dT)/60/10))) #get number of 10 minute buckets
beh.bins = DIutils.binning.digitize(beh.time, bins);
checkbins(:beh)
DIutils.filtreg.register(beh, cycles, on="time",
    transfer=["epoch","bins"]);
checkbins(:cycles);
df.cycle = df.cyc_central;
DIutils.filtreg.register(cycles, df, on="cycle",
    transfer=["epoch","time","bins"]);
checkbins(:df);
DIutils.filtreg.register(cycles, Rdf, on="cycle", transfer=["epoch","bins"]);
checkbins(:Rdf);
DIutils.filtreg.register(cycles, Rdf, on="time", transfer=["epoch","bins"]);
checkbins(:Rdf)

using ColorSchemes
function plot_bins(df, x=nothing; cmap=:Reds)
    y = :bins
    vars = x === nothing ? [y] : [x,y] 
    dfc = dropmissing(df, vars) 
    x = x === nothing ? (1:size(dfc,1)) : dfc[!,x]
    c=DIutils.plotutils.getplotcolor(dfc[!,y], cmap)
    y = dfc[!,y]
    scatter(x, y; c )
end
plot(
    plot_bins(beh, :time),
    plot_bins(df, :time)
)

# CORRELATION STRUCTURE OF 10 MINUTE BINS
using LinearAlgebra
function test_correlation_structure(prop=:bins; desc="", ploton=true)
    SIM=[]
    DF = groupby(df, prop)
    prog =Progress(length(DF);desc="CC")
    CC = Vector{Array}(undef, length(DF))
    for (i,d) in enumerate(DF)
        CC[i] = begin
            Q = get_correlation_matrices(d).C
            Q[isnan.(Q)] .= 0
            next!(prog)
            Q[:]
        end
    end

    a2d = DIutils.arr.atleast2d
    println("Computing simlarity")
    similarity = Matrix{Float64}(undef, length(CC), length(CC))
    prog = Progress(length(CC)^2)
    Threads.@threads for (i,(q1,q2)) in 
        collect(enumerate(Iterators.product(CC,CC)))
        similarity[i]= (a2d(q1)' * a2d(q2))[1]/(norm(q1)*norm(q2))
        next!(prog)
    end
    push!(SIM,similarity)
    println("PLotting CC")
    hCC = ploton ? heatmap(similarity, 
        title="Correlation structure over time bins") : nothing
    ploton ? savefig(plotsdir("isolated","glm",
        "correlation structure of $desc $prop for cyclemean df.pdf")) : nothing

    println("Processing Rdf")
    RDF = groupby(Rdf, prop)
    prog =Progress(length(RDF);desc="RCC")
    println("Computing simlarity")
    RCC = Vector{Array}(undef, length(RDF))
    Threads.@threads for (i,d) in collect(enumerate(RDF))
        RCC[i] = begin
            Q = get_correlation_matrices(begin
                unstack(d, :time, :unit, :value, combine=mean) 
                end).C
            Q[isnan.(Q)] .= 0
            next!(prog)
            Q[:]
        end
    end

    println("Computing RCC similarity")
    similarity = Matrix{Float64}(undef, length(RCC), length(RCC))
    prog = Progress(length(CC)^2)
    Threads.@threads for (i,(q1,q2)) in 
        collect(enumerate(Iterators.product(RCC,RCC)))
        similarity[i]= (a2d(q1)' * a2d(q2))[1] / (norm(q1) * norm(q2))
        next!(prog)
    end
    push!(SIM,similarity)
    println("PLotting RCC")
    hRCC = ploton ? heatmap(similarity, 
        title="Correlation structure over time bins") : nothing
    ploton ? savefig(plotsdir("isolated","glm",
        "correlation structure of $desc $prop for cyclemean Rdf.pdf")) : nothing

    println("PLotting CC vs RCC")
    similarity = Matrix{Float64}(undef, length(CC), length(RCC))
    prog = Progress(length(CC)^2)
    Threads.@threads for (i,(q1,q2)) in 
        collect(enumerate(Iterators.product(CC,RCC)))
        similarity[i]= (a2d(q1)' * a2d(q2))[1]/(norm(q1)*norm(q2))
        next!(prog)
    end
    push!(SIM,similarity)
    hCCRCC = ploton ? heatmap(similarity, xlabel="cycle df", 
        ylabel="firing rate df",
        title="Correlation structure over time bins") : nothing
    return (;hCC,hRCC, hCCRCC, SIM)
    # (;SIM)
end
test_correlation_structure(:cycle; desc="", ploton=false)

outs = test_correlation_structure(:bins; desc="10 minute", ploton=true)
 # heatmap(cor(Matrix(df[!,uicellcols(df)])), c=:delta, clim=(0,0.1))


    cycles[!,:hasmis] = vec(any(Matrix(ismissing.(cycles[!,matchprops])),
        dims=2))
    cycles.hasocc = (!).(ismissing.(occ.datainds))
    @assert size(occ.datainds,1) == size(cycles.hasocc,1)

    q1= @df dropmissing(cycles,[matchprops..., :bins]) plot(:time,:bins, 
        label="cycles have behavior and bins")
    @df dropmissing(df,[:bins]) scatter!(:time,:bins,label="df have bins")
    x = cycles.time[abs.([0;diff((!).(cycles.hasmis))]) .> 0]
    vspan!(x, alpha=0.40, labels="cycles has non miss", c=:blue)

    x=cycles.time[abs.([0;diff(cycles.hasocc)]) .> 0]
    q2=vspan!(x, alpha=0.40, labels="cycles has occ", c=:red)
    plot(q1,q2)

    Z=ismissing.([
        Matrix(@subset(cycles, :cycle .∈ (iso_cycles,))[:,[:x,:stopWell,:speed]]) occ.datainds[iso_cycles,:]])
    heatmap(Z,xticks=([1,2,3,4],["x","stopwell","speed","datainds"]))

using .Threads
Base.throwto.(collect(1:16), [InterruptException()])
