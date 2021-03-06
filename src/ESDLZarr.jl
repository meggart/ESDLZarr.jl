module ESDLZarr
using ESDL
import ZarrNative: ZGroup, zopen
import ESDL.CubeAPI: getNv, gettoffsnt, getvartype, Cube
import ESDL.Cubes: cubechunks, iscompressed
using Dates

struct ZarrConfig
  start_time::Date
  end_time::Date
  grid_x0::Int
  grid_y0::Int
  grid_height::Int
  grid_width::Int
  spatial_res::Float64
  temporal_res::Int
  chunk_sizes::NTuple{3,Int}
  compressed::Bool
end
function ZarrConfig(g::ZGroup,firstvar)
  timeax = g["time"]
  timevals = timeax[:]
  refdate = Date(split(timeax.attrs["units"]," ")[3])
  start_time = refdate+Dates.Day(first(timevals))-Day(4)
  end_time = refdate+Dates.Day(last(timevals))-Day(4)
  lonax = g["lon"]
  latax = g["lat"]
  sres = abs(latax[2][1]-latax[1][1])
  tres = timevals[2][1]-timevals[1][1]
  cs = g[firstvar].metadata.chunks
  compressed = !isa(g[firstvar].metadata.compressor,ZarrNative.NoCompressor)
  ZarrConfig(start_time,end_time,0,0,length(latax),2*length(latax),sres,tres,cs,compressed)
end

struct ZarrCube <: ESDL.CubeAPI.UCube
  path::String
  group::ZGroup
  dataset_files::Vector{String}
  var_name_to_var_index::Dict{String,Int}
  config::ZarrConfig
end
function ZarrCube(p::String)
  Cube(zopen(p))
end
function Cube(zg::ZGroup)
  spattempars = filter(v->ndims(v[2])==3,zg.arrays)
  varlist = collect(keys(spattempars))
  vni = Dict(i[2]=>i[1] for i in enumerate(varlist))
  ZarrCube(zg.storage.folder,zg,varlist,vni,ZarrConfig(zg,varlist[1]))
end
getvartype(c::ZarrCube, n::String)=eltype(c.group[n])
cubechunks(c::ZarrCube)=c.config.chunk_sizes
iscompressed(c::ZarrCube)=c.config.compressed
function ESDL.Cubes._read(s::ESDL.CubeAPI.SubCube{<:Any,ZarrCube},t::AbstractArray,r::CartesianIndices)
  grid_y1,grid_y2,grid_x1,grid_x2 = s.sub_grid
  y1,i1,y2,i2,ntime,NpY           = s.sub_times
  toffs,nt= gettoffsnt(s,r)
  toffs = (y1 - year(s.cube.config.start_time))*NpY + i1 + toffs
  rcor = CartesianIndices((r.indices[1].+(grid_x1-1), r.indices[2].+(grid_y1-1), (toffs:toffs.+nt-1)))
  #voffs,nv = getNv(r)
  singvar_zarr(t,s.cube,s.variable,rcor)
end
function ESDL.Cubes._read(s::ESDL.CubeAPI.SubCubeV{<:Any,ZarrCube},t::AbstractArray,r::CartesianIndices)
  grid_y1,grid_y2,grid_x1,grid_x2 = s.sub_grid
  y1,i1,y2,i2,ntime,NpY           = s.sub_times

  toffs,nt= gettoffsnt(s,r)
  toffs = (y1 - year(s.cube.config.start_time))*NpY + i1 + toffs
  voffs,nv = getNv(r)
  #@show (r.indices[1].+(grid_x1-1), r.indices[2].+(grid_y1-1), toffs:toffs+nt-1)
  rcor = CartesianIndices((r.indices[1].+(grid_x1-1), r.indices[2].+(grid_y1-1), toffs:toffs+nt-1))
  #@show rcor.indices
  for (iiv,iv)=enumerate((voffs+1):(voffs+nv))
    singvar_zarr(view(t,:,:,:,iiv),s.cube,s.variable[iv],rcor)
  end
end
function ESDL.Cubes._read(s::ESDL.CubeAPI.SubCubeStatic{<:Any,ZarrCube},t::AbstractArray,r::CartesianIndices)
  grid_y1,grid_y2,grid_x1,grid_x2 = s.sub_grid
  y1,i1,y2,i2,ntime,NpY           = s.sub_times
  toffs = (y1 - year(s.cube.config.start_time))*NpY + i1
  rcor = CartesianIndices((r.indices[1].+(grid_x1-1), r.indices[2].+(grid_y1-1),(toffs):(toffs)))
  singvar_zarr(t,s.cube,s.variable,rcor)
end

import ZarrNative
function singvar_zarr(outar,cube,variable,r)
  ZarrNative.readblock!(outar, cube.group[variable],r)
  mv  = cube.group[variable].metadata.fill_value
  if !isa(mv,Nothing)
    replace!(outar,mv=>missing)
  end
end

end # module
