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
# load into a usable object
bg = Background(dom_master, top)
sp = SpeciesInfo(species, ele_range, geo_range, stand_range)
nm = NullModeller(dom_master)

############################ run the null model 

example_species="Phlogophilus harterti" # name of one of the example species from the nam object
rs_std=false # should the null model use the standardized range size (true) or the empirical range size (false)
nrep=10 # nuber of repetitions

#### empirical range

function decode_range(presences, domain)
    ret = falses(dims(domain))
    ret[presences] .= true
    ret .&= domain
    ret
end

function prep_map(res_nm,dom;trim_map=true,crop_to_ext=nothing)
    map_nm = decode_range(res_nm, dom)
    if crop_to_ext === nothing
        if trim_map
            map_nm=Rasters.trim(map_nm,pad=10)
        end
    else
        map_nm = Rasters.crop(map_nm, to=crop_to_ext)
    end
    plot(map_nm)
    map_nm
end

map_emp = decode_range(geo_range[example_species], dom_master)
plot(map_emp)

cut(x) = x[400:480, 320:400]

results = []
for cut in (false, true), split in (false, true)
    res = null_model!(nm, example_species, sp, bg; split_groups = split, cut_domain = cut)
    push!(results, copy(res))
end

plot([plot(cut(res), title = "Null model $i", ticks = false) for (i, res) in enumerate(results)]...)

