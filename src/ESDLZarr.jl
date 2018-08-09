module ESDLZarr
using ESDL
import Zarr: ZarrGroup, zopen
import ESDL.CubeAPI: getNv, gettoffsnt, getvartype, Cube
import Base.Dates: Day, year

struct ZarrConfig
  start_time::Date
  end_time::Date
  grid_x0::Int
  grid_y0::Int
  grid_height::Int
  grid_width::Int
  spatial_res::Float64
  temporal_res::Int
end
function ZarrConfig(g::ZarrGroup)
  timeax = g["time"]
  timevals = timeax[:]
  refdate = Date(split(timeax.atts["units"]," ")[3])
  start_time = refdate+Base.Dates.Day(first(timevals))-Day(4)
  end_time = refdate+Base.Dates.Day(last(timevals))-Day(4)
  lonax = g["lon"]
  latax = g["lat"]
  sres = abs(latax[2]-latax[1])
  tres = timevals[2]-timevals[1]
  ZarrConfig(start_time,end_time,0,0,length(latax),2*length(latax),sres,tres)
end

struct ZarrCube <: ESDL.CubeAPI.UCube
  path::String
  group::ZarrGroup
  dataset_files::Vector{String}
  var_name_to_var_index::Dict{String,Int}
  config::ZarrConfig
end
function ZarrCube(p::String)
  Cube(zopen("/home/fgans/zarr/esdc-8d-0.25deg-1x720x1440-1.0.1_1_zarr","r"))
end
function Cube(zg::ZarrGroup)
  vni = Dict(i[2]=>i[1] for i in enumerate(zg.arrays))
  ZarrCube(zg.g[:_store][:path],zg,zg.arrays,vni,ZarrConfig(zg))
end
getvartype(c::ZarrCube, n::String)=eltype(c.group[n])
function ESDL.Cubes._read(s::ESDL.CubeAPI.SubCube{<:Any,ZarrCube},t::Tuple,r::CartesianRange)
  outar,mask=t
  grid_y1,grid_y2,grid_x1,grid_x2 = s.sub_grid
  y1,i1,y2,i2,ntime,NpY           = s.sub_times

  grid_x1 = grid_x1 + r.start.I[1] - 1
  nx      = r.stop.I[1] - r.start.I[1] + 1
  grid_y1 = grid_y1 + r.start.I[2] - 1
  ny      = r.stop.I[2] - r.start.I[2] +1
  toffs,nt= gettoffsnt(s,r)
  toffs = (y1 - year(s.cube.config.start_time))*NpY + i1 + toffs
  @show (toffs,nt)
  #voffs,nv = getNv(r)
  singvar_zarr(outar,s.cube,s.variable,grid_x1,nx,grid_y1,ny,toffs,nt)
  lsmask = s.cube.group["water_mask"][grid_x1:grid_x1+nx-1, grid_y1:grid_y1+ny-1,1:1]
  broadcast!((v,m)->((UInt8(m)-0x01)*0x05) | isnan(v),mask,outar,lsmask)
end
function ESDL.Cubes._read(s::ESDL.CubeAPI.SubCubeV{<:Any,ZarrCube},t::Tuple,r::CartesianRange)
  outar,mask=t
  grid_y1,grid_y2,grid_x1,grid_x2 = s.sub_grid
  y1,i1,y2,i2,ntime,NpY           = s.sub_times

  grid_x1 = grid_x1 + r.start.I[1] - 1
  nx      = r.stop.I[1] - r.start.I[1] + 1
  grid_y1 = grid_y1 + r.start.I[2] - 1
  ny      = r.stop.I[2] - r.start.I[2] +1
  toffs,nt= gettoffsnt(s,r)
  toffs = y1 - year(s.cube.config.start_time)*NpY + i1 + toffs
  voffs,nv = getNv(r)
  lsmask = reshape(s.cube.group["water_mask"][grid_x1:grid_x1+nx-1, grid_y1:grid_y1+ny-1,1],nx,ny,1,1)
  for iv=(voffs+1):(voffs+nv)
    singvar_zarr(view(outar,:,:,:,iv),s.cube,s.variable[iv],grid_x1,nx,grid_y1,ny,toffs,nt)
  end
  broadcast!((v,m)->((UInt8(m)-0x01)*0x05) | isnan(v),mask,outar,lsmask)
end

function singvar_zarr(outar,cube,variable,grid_x1,nx,grid_y1,ny,toffs,nt)
  outar[:,:,:] = cube.group[variable][grid_x1:grid_x1+nx-1, grid_y1:grid_y1+ny-1,toffs:toffs+nt-1]
end

end # module
