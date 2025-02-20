module AvianRangeShapes

using ArchGDAL,Rasters, DataFrames,Plots, SpreadingDye, NearestNeighbors, Images, JLD2, SkipNan


include("functions.jl")

export null_models, prep_map
end