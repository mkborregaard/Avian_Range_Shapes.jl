using ArchGDAL, Rasters, DataFrames, Plots, JLD2
using AvianRangeShapes

#loading geographical domain
dom_master = .!ismissing.(Raster("data/sf1_mainland.tif") )

#loading topographical raster
top = resample(Raster("data/top_q_proj.tif"); to = dom_master)
top = Float32.(replace!(top, missing => NaN))

#loading species elevational range limits
elv=load_object("data/elevational range limits.jld2")
ele_range = Dict(r.Species => (min = r.minimum_elevation, max = r.Maximu_elevation) for r in eachrow(elv))

#loading the geographic geographic ranges of six example species
dis_elv=load_object("data/bird ranges.jld2")
species=["Acropternis orthonyx","Aglaeactis castelnaudii","Coeligena lutetiae","Diglossa mystacalis","Phlogophilus harterti","Scytalopus griseicollis"]
geo_range = Dict(zip(species, dis_elv))

#loading standardized range sizes
formated_rs=load_object("data/standardized_range_sizes.jld2")
stand_range = Dict(r.nam => r.rank_range for r in eachrow(formated_rs))

# initialize objects for analysis
bg = Background(dom_master, top[Band = 1], top[Band = 2])
sp = SpeciesInfo(species, ele_range, geo_range, stand_range)
nm = NullModeller(dom_master)

############################ run the null model 

example_species="Phlogophilus harterti" # name of one of the example species from the nam object
rs_std=false # should the null model use the standardized range size (true) or the empirical range size (false)
nrep=10 # nuber of repetitions

#### empirical range
map_emp = decode_range(geo_range[example_species], dom_master)
plot(map_emp)

# run all four null models

results = [null_model!(nm, example_species, sp, bg; split_groups = split, cut_domain = cut) 
                for cut in (false, true), split in (false, true)]

# and plot the results
plot([plot(res[400:480, 320:400], title = "Null model $i", ticks = false) for (i, res) in enumerate(results)]...)

