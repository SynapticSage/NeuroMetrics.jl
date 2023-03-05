include("../imports_isolated.jl")

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
    cu(cells::DataFrame, df::DataFrame)  = @subset(cells, :unit .∈ (uc(df),))
    icu(cells::DataFrame, df::DataFrame) = @subset(cells, :unit .∈ (ic(df),))
    cu(cells::Vector{String}, df::DataFrame)  = cu(parse.(Int,cells), df)
    icu(cells::Vector{String}, df::DataFrame) = icu(parse.(Int,cells), df)
    cu(cells::Vector{Int}, df::DataFrame)  = cells .∈ (ic(df),)
    icu(cells::Vector{Int}, df::DataFrame)  = cells .∈ (uc(df),)
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
    ca1pct = sum(ca1)/length(ca1) * Plots.pct
    pfcpct = sum(pfc)/length(pfc) * Plots.pct
    (;C, Cca1, Cpfc, Cca1pfc)
end

function plot_correlation_matrices(;Cca1, Cpfc, Cca1pfc, kws...)
    g= grid(2,2)
    g.widths  = [ca1pct, pfcpct]
    g.heights = [ca1pct, pfcpct]
    m = 0.2
    plot(
        heatmap(Cca1,             colorbar_tickfontsize=3, colorbar_title="CA1",     c=:PiYG_5,    clim=(-m, m),  colorbar=:left),
        heatmap(Cca1pfc,          colorbar_tickfontsize=3, colorbar_title="CA1-PFC", c=:vik,       clim=(-m, m)),
        Plot.create_blank_plot(),
        heatmap(Cpfc,             colorbar_tickfontsize=3, colorbar_title="PFC",     c=:RdYlGn_11, clim=(-m, m)),
        layout=g,
    )
end


overall = get_correlation_matrices(df)
plot_correlation_matrices(;overall...)


# Prepare for a 10 minute bucketed approach
dT = diff(beh.time)
bins = Int(maximum(floor.(cumsum(dT)/60/10))) #get number of 10 minute buckets
beh.bins = DIutils.binning.digitize(beh.time, bins)
DIutils.filtreg.register(beh, cycles, on="time", transfer=["epoch","bins"])
df.cycle = df.cycs
DIutils.filtreg.register(cycles, df, on="cycle", transfer=["epoch","time","bins"])

@df df scatter(:cycle)

 # heatmap(cor(Matrix(df[!,uicellcols(df)])), c=:delta, clim=(0,0.1))
