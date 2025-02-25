module AvianRangeShapes

using ArchGDAL,Rasters, DataFrames,Plots, SpreadingDye, NearestNeighbors, Images, JLD2, SkipNan


include("functions.jl")

export null_model!, prep_map, decode_range, Background, SpeciesInfo, NullModeller
end