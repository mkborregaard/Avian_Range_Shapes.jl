using SpreadingDye, NearestNeighbors, StatsBase, Rasters, ImageMorphology

# some convenience structs to hold the data
struct Background
    domain::Raster{Bool}
    elevation::Raster
end
Base.show(io::IO, x::Background) = print("A Background object with `domain` and `elevation` rasters of size $(size(domain, 1)) x $(size(domain, 2))")

struct SpeciesInfo
    names::Vector{String}
    ele_range::Dict{String, @NamedTuple{min::Int64, max::Int64}}
    geo_range::Dict{String, Vector{Int}}
    stand_range::Dict{String, Int}
end
Base.show(io::IO, x::SpeciesInfo) = print("A SpeciesInfo object with names and elevational ranges of species, along with their geographic ranges and standardized range sizes")


struct NullModeller
    domain::Raster{Bool}
    emp::Raster{Bool}
    final_sim::Raster{Bool}
    patch_sim::Raster{Bool}
    patches::Raster{Int}
end
Base.show(io::IO, x::NullModeller) = print("A NullModeller object with intermediate rasters for the empirical and simulated ranges of a species")


NullModeller(r::Raster) = NullModeller(copy(r), falses(dims(r)), falses(dims(r)), falses(dims(r)), zeros(Int, dims(r)))

function decode_range!(ret, presences, domain)
    fill!(ret, false)
    ret[presences] .= true
    ret .&= domain
    ret
end
decode_range(presences, domain) = decode_range!(falses(dims(domain)), presences, domain)

# Filters the geographic domain `dom` by the species elevational range limits
cut_elevation!(dom, top, min, max) = (dom .&= top[Band = 1] .< max .&& top[Band = 2] .> min)

"""
    find_groups(emp2, max_dist=5, min_prop=0.1)

Identify groups of isolated range patches.

Arguments:
    - max_dist: maximum distance between range patches
    - min_prop: minimum patch size relative to the species’ largest patch 
"""
function find_groups(emp2, max_dist=5, min_prop=0.1) 
    #Identify groups of isolated range patches 
   groups = collect_groups(emp2)
   groups=groups[findall(length.(groups).>0)]

   #### algorithm treating tiny range patches as extensions of the adjacent larger coherent range rather than independent range patches
   groups=join_neighbours(groups; max_dist, min_prop) 
   groups
end

#standardizing range sizes of groups
#if the standardized range size is smaller than the empirical, randomly subtract grid cells from the groups in stepwise fashion, weighted by the patch size
#i.e. large patches have greater chance of beeing modified by the standardization
function update_group_size!(new_range,total_rangesize,group_size)   
    dif = new_range - total_rangesize
    group_size .+= sign(dif) .* sample(length(group_size), Weights(group_size), abs(dif))
end

spreading_dye_patches!(nm::NullModeller) = spreading_dye_patches!(nm.final_sim, nm.patch_sim, nm.patches, nm.domain)
function spreading_dye_patches!(final_sim::Raster{Bool}, patch_sim::Raster{Bool}, patches::Raster{Int}, dom::Raster{Bool})
    final_sim .= false
    finalrange = count(!=(0), patches)
    for i in 1:maximum(patches) # 0 is outside
        patch_sim .= patches .== i
        patchsize = sum(patch_sim)
        spreading_dye!(patch_sim, patchsize, dom, random_point_on_domain(patch_sim))
        final_sim .|= patch_sim
    end
    sum(final_sim) < finalrange && SpreadingDye.expand_spreading!(final_sim, finalrange - sum(final_sim), dom)
    final_sim
end

### constructing neighborhood matrix for the range patches
function reldists(point)
    a = fill(Inf, length(point), length(point))
    sort!(point, by = length, rev = true)
    for i in 1:length(point)-1
        m1 = Float64.(stack(point[i]))
        kdtree = KDTree(m1; leafsize = 8)
        for j in i+1:length(point)
            m2 = Float64.(stack(point[j]))
            ind, dist = knn(kdtree, m2, 1)
            a[i, j] = only(minimum(dist)); a[j, i] = a[i,j]
        end
    end
    a
end

#parameters: maximum distance between patches and minimum percentage size of patches
function join_neighbours(groups;max_dist::Int64=5,min_prop::Float64=0.1)
    #organizing patch groups
    groups=sort(groups, by = length, rev = true)
    group_size=length.(groups)
    group_size_prop=group_size./maximum(group_size)

    ## convert raster ID's to matrix coordinates
    point=Any[]
    for t in 1:length(groups)
        zz=copy(Float64.(dom));zz.=NaN;zz[groups[t]].=1;push!(point,Tuple.(collect(CartesianIndices(zz))[isfinite.(zz)]))
    end

    #construct neigbborhood matrix
    dm=reldists(point)
    p = Set.(1:length(groups))
    ll = group_size
    id=copy(dm);id[:].=1:length(id[:])

    # identify small patches
    small=reverse(findall(group_size_prop.<min_prop)) # reverse to start merging from the smallest patch to the largest

    while length(small)>0
        po=small
        mer=Tuple.(argmin(dm[po,:])) # picks first element (will prioritize merging with the largest patch)
        mer2=Tuple.(findall(id.==id[po,:][mer...])[1])
        mi=minimum(mer2);ma=maximum(mer2)

        if dm[mer2...]<max_dist
            p[mi] = union(p[mi], p[ma])
            empty!(p[ma])
            ll[mi]=ll[ma]+ll[mi]; ll[ma]=0
            
            comp=[reshape(dm[mi,:], 1, :); reshape(dm[ma,:], 1, :)]
            mins=minimum.(eachcol(comp))
            infs=.!isfinite.(dm)
            dm[ma,:]=mins
            dm[:,mi]=mins
            dm[infs].=Inf
        else ll[ma]=0 end

        dm[mi,ma]=Inf
        dm[ma,mi]=Inf
        group_size_prop=ll./maximum(ll)
        small= reverse(findall((group_size_prop.<min_prop) .& (group_size_prop.>0)))# reverse to start merging from the smallest patch to the largest
    end

    out=Any[]
    ss=findall(length.(p).>0)
    for a in 1:length(ss)
        push!(out,reduce(vcat, groups[collect(p[ss[a]])]))
    end

    out
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

"""
    null_model!(nm::NullModeller, species::String, sp::SpeciesInfo, bg::Background, nrep::Int64=1; 
    cut_domain::Bool, split_groups::Bool, rs_std::Bool)

The outward_facing function to run each null model

Arguments:
    - nm:               an object containing the null model intermediates
    - species:          the name of the species
    - sp:               an object containing all data on species ranges
    - bg:               an object of regional data rasters, domain, elevation etc
    - nrep:             how many times to redo the simulation
    - cut_domain:       should the domain be adapted to the elevational range of the species?
    - split_groups:     should non-cohesive ranges remain split?
    - rs_std::Bool:     should the null model use the standardized range size (true) or the empirical range size (false)
"""
function null_model!(nm::NullModeller, species::String, sp::SpeciesInfo, bg::Background, nrep::Int64=1; 
    cut_domain::Bool=false, split_groups::Bool=false, rs_std::Bool=false)

    # initialise the domain and possibly filter it by the species elevational range limits     
    nm.domain .= bg.domain
    cut_domain && cut_elevation!(nm.domain, bg.elevation, sp.ele_range[species]...)

    #constructing a raster of the species empirical range
    decode_range!(nm.emp, sp.geo_range[species], nm.domain)

    # define groups
    nm.patches .= split_groups ? label_components!(nm.patches, nm.emp, strel_box((3, 3))) : nm.emp # queen style neighbourhood

   # standardizing the range size frequency distribution
    if rs_std
        new_range = stand_range[species]
        update_group_size!(new_range,total_rangesize,group_size)    
    end

    spreading_dye_patches!(nm)
    nm.final_sim
end

"""
    null_model!(final_sim::Raster{Bool}, patch_sim::Raster{Bool}, patches::Raster{Int}, 
    dom::Raster{Bool}, domain_master::Raster{Bool}, species::String, cut_domain::Bool, split_groups::Bool, 
    geo_range::Dict, rs_std::Bool, top::Raster, ele_range::Dict, stand_range::Dict, nrep::Int64)

The outward_facing function to run each null model

Arguments:
    - final_sim:        a target raster for the simulated range
    - patch_sim:        an intermediate convenience raster for simulating multiple patches
    - patches:          a raster to identify separate patches
    - dom:              an intermediate raster to hold the smaller domain for some species
    - domain_master     the biogeographical domain
    - species::String:  the name of the species
    - cut_domain:       should the domain be adapted to the elevational range of the species?
    - split_groups:     should non-cohesive ranges remain split?
    - geo_range::Dict:  grid cell ids comprising the species empirical range
    - rs_std::Bool:     should the null model use the standardized range size (true) or the empirical range size (false)
    - top:              topographical raster
    - ele_range:        dict with the species' elevational range limits
    - stand_range:      dict with standardized range sizes (only used if rs_std=true)
    - nrep:             number of repetitions
"""
function null_model!(final_sim::Raster{Bool}, patch_sim::Raster{Bool}, patches::Raster{Int}, 
    dom::Raster{Bool}, domain_master::Raster{Bool}, species::String, cut_domain::Bool, split_groups::Bool, 
    geo_range::Dict, rs_std::Bool, top::Raster, ele_range::Dict, stand_range::Dict, nrep::Int64)

    # initialise the domain and possibly filter it by the species elevational range limits     
    dom .= domain_master
    cut_domain && cut_elevation!(dom, top, ele_range[species]...)

    #construct a raster of the species empirical range
    emp = decode_range(geo_range[species], dom)

    # define patches
    patches .= split_groups ? label_components!(patches, emp, strel_box((3, 3))) : emp # queen style neighbourhood

   # standardizing the range size frequency distribution
    if rs_std
        new_range = stand_range[species]
        update_group_size!(new_range,total_rangesize,group_size)    
    end

    spreading_dye_patches!(final_sim, patch_sim, patches, dom)
end