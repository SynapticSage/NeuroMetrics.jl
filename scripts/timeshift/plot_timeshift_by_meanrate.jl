
cells = raw.load_cells(animal, day, "*")

using LazyGrids: ndgrid

x = replace(cells.origunit, missing=>NaN)
y = replace(cells[!,colorby], missing=>NaN)
z = replace(cells.meanrate, missing=>NaN)

Plots.scatter(x,y,z,
              markercolor=get.([ColorSchemes.vik], z),
                               xlabel="addition-order", ylabel="shift", zlabel="firing rate")

xx, zz = ndgrid(nanminimum(x):0.1:nanmaximum(x), nanminimum(z):0.1:nanmaximum(z))
xx = Array(xx)
yy = Array(yy)
yy = ones(size(xx)).*0.00
Plots.surface!(vec(xx),vec(yy),vec(zz))

Plots.savefig(plotsdir("timeshift", "addition-order_shiftOptimal_rate_(shift in minutes).png"))
Plots.savefig(plotsdir("timeshift", "addition-order_shiftOptimal_rate_(shift in minutes).pdf"))


x = cells.origunit
y = cells[!,colorby]
z = cells.spikewidth

Plots.scatter(x,y,z,
              markercolor=get.([ColorSchemes.vik], cells.meanrate),
                               xlabel="addition-order", ylabel="shift", zlabel="firing rate")

xx, zz = ndgrid(nanminimum(x):0.1:nanmaximum(x), nanminimum(z):0.1:nanmaximum(z))
yy = ones(size(xx)).*0.00
Plots.surface!(xx,yy,zz)
