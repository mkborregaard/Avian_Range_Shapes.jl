using SpreadingDye, NearestNeighbors, SkipNan, StatsBase, Rasters, ImageMorphology

# some convenience structs to hold the data
struct Background
    domain::Raster{Bool}
    elevation::Raster
end

struct SpeciesInfo
    names::Vector{String}
    ele_range::Dict{String, @NamedTuple{min::Int64, max::Int64}}
    geo_range::Dict{String, Vector{Int}}
    stand_range::Dict{String, Int}
end

struct NullModeller
    domain::Raster{Bool}
    final_sim::Raster{Bool}
    patch_sim::Raster{Bool}
    patches::Raster{Int}
end

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

#stadardizing range sizes of groups
function update_group_size!(new_range,total_rangesize,group_size)   
    dif = new_range - total_rangesize
    group_size .+= sign(dif) .* sample(length(group_size), Weights(group_size), abs(dif))
end

function spreading_dye_patches!(final_sim::Raster{Bool}, patch_sim::Raster{Bool}, patches::Raster{Int}, dom::Raster{Bool})
    final_sim .= false
    finalrange = count(!=(0), patches)
    for i in 1:maximum(patches) # 0 is outside
        patch_sim .= patches .== i
        patchsize = sum(patches)
        spreading_dye!(patch_sim, patchsize, dom, random_point_on_domain(patch_sim))
        final_sim .|= patch_sim
    end
    sum(final_sim) < finalrange && expand_spreading!(final_sim, total_rangesize - sum(final_sim), dom)
end

"""
    null_models(species::String, geo_range::Raster{Bool}, rs_std::Bool, dom::Raster{Bool}, top::Raster, elv::Dict,
    formated_rs::Dict, nrep::Int64)

The outward_facing function to run each null model

Arguments:
    - species::String:  state name of example species
    - geo_range::Dict:  grid cell ids comprising the species empirical range
    - rs_std::Bool:     should the null model use the standardized range size (true) or the empirical range size (false)
    - dom:              biogeographical domain
    - top:              topographical raster
    - elv:              dict with the species' elevational range limits
    - formated_rs:      dict with standardized range sizes (only used if rs_std=true)
    - nrep:             number of repetitions
"""
function null_model!(final_sim::Raster{Bool}, patch_sim::Raster{Bool}, patches::Raster{Int}, 
    dom::Raster{Bool}, domain_master::Raster{Bool}, species::String, cut_domain::Bool, split_groups::Bool, 
    rs_std::Bool, top::Raster, ele_range::Dict, stand_range::Dict, nrep::Int64)

    # initialise the domain and possibly filter it by the species elevational range limits     
    dom .= domain_master
    cut_domain && cut_elevation!(dom, top, ele_range[species]...)

    #constructing a raster of the species empirical range
    emp = decode_range(geo_range[species], dom)

    # define groups
    if split_groups 
        label_components!(patches, emp, strel_box((3, 3))) # queen style neighbourhood
    else
        patches .= emp
    end

    ######### ## standardizing the range size frequency distribution
    if rs_std
        #if the standardized range size is smaller than the empirical, randomly subtract grid cells from the groups in stepwise fashion, weighted by the patch size
        #i.e. large patches have greater chance of beeing modified by the standardization
        new_range=formated_rs[formated_rs.nam.==species,:rank_range][1]
        update_group_size!(new_range,total_rangesize,group_size)    
        total_rangesize =new_range
        ppo=findall(group_size.>0)
        group_size=group_size[ppo]
        groups=groups[ppo]
    end

    spreading_dye_patches!(final_sim, patch_sim, patches, dom)
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


