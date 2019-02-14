module ESDLZarr
using ESDL
import ZarrNative: ZGroup, zopen, ZArray
import ESDL.CubeAPI: getNv, gettoffsnt, getvartype, Cube
import ESDL.Cubes: cubechunks, iscompressed, AbstractCubeData, getCubeDes,
  caxes,chunkoffset, gethandle, subsetCube, axVal2Index, findAxis
import ESDL.Cubes.Axes: axname
import Dates: Day,Hour,Minute,Second,Month,Year, Date
const spand = Dict("days"=>Day,"months"=>Month,"years"=>Year,"seconds"=>Second,"minutes"=>Minute)

struct ZArrayCube{T,N,A<:ZArray{T,N},S} <: AbstractCubeData{T,N}
  a::A
  axes::Vector{CubeAxis}
  subset::S
end
import ZarrNative: readblock!
getCubeDes(::ZArrayCube)="ZArray Cube"
caxes(z::ZArrayCube)=z.axes
iscompressed(z::ZArrayCube)=!isa(z.a.metadata.compressor,ZarrNative.NoCompressor)
cubechunks(z::ZArrayCube)=z.a.metadata.chunks
chunkoffset(z::ZArrayCube{<:Any,N}) where N=ntuple(i->0,N)
Base.size(z::ZArrayCube) = length.(z.subset)
Base.size(z::ZArrayCube{<:Any,<:Any,<:ZArray,Nothing}) = size(z.a)
#ESDL.Cubes.gethandle(z::ZArrayCube) = z.a
function ESDL.Cubes._read(z::ZArrayCube{<:Any,N,<:Any,<:Nothing},thedata::AbstractArray{<:Any,N},r::CartesianIndices{N}) where N
  readblock!(thedata,z.a,r)
end

#Helper functions for subsetting indices
_getinds(s1,s,i) = s1[firstarg(i...)],Base.tail(i)
_getinds(s1::Int,s,i) = s1,i
function getsubinds(subset,inds)
    el,rest = _getinds(firstarg(subset...),subset,inds)
    (el,getsubinds(Base.tail(subset),rest)...)
end
getsubinds(subset::Tuple{},inds) = ()
firstarg(x,s...) = x


function ESDL.Cubes._read(z::ZArrayCube{<:Any,N,<:Any},thedata::AbstractArray{<:Any,N},r::CartesianIndices{N}) where N
  allinds = CartesianIndices(map(Base.OneTo,size(z.a)))
  subinds = map(getindex,allinds.indices,z.subset)
  @show subinds, r.indices
  r2 = getsubinds(subinds,r.indices)
  readblock!(thedata,z.a,CartesianIndices(r2))
end

function _write(y::ZArrayCube{<:Any,N,<:Any,<:Nothing},thedata::AbstractArray,r::CartesianIndices{N}) where N
  readblock!(thedata,z.a,r,readmode=false)
end

function infervarlist(g::ZGroup)
  dimsdict = Dict{Tuple,Vector{String}}()
  foreach(g.arrays) do ar
    k,v = ar
    vardims = reverse((v.attrs["_ARRAY_DIMENSIONS"]...,))
    haskey(dimsdict,vardims) ? push!(dimsdict[vardims],k) : dimsdict[vardims] = [k]
  end
  filter!(p->!in("bnds",p[1]),dimsdict)
  llist = Dict(p[1]=>length(p[2]) for p in dimsdict)
  _,dims = findmax(llist)
  varlist = dimsdict[dims]
end

function parsetimeunits(unitstr)
    unitstr = "days since 1982-01-05"
    re = r"(\w+) since (\d\d\d\d)-(\d\d)-(\d\d)"

    m = match(re,unitstr)

    refdate = Date(map(i->parse(Int,m[i]),2:4)...)
    refdate,spand[m[1]]
end
function toaxis(dimname,g)
    axname = dimname in ("lon","lat","time") ? uppercasefirst(dimname) : dimname
    ar = g[dimname]
    if axname=="Time" && haskey(ar.attrs,"units")
        refdate,span = parsetimeunits(ar.attrs["units"])
        tsteps = refdate.+span.(ar[:])
        TimeAxis(tsteps)
    else
      axdata = testrange(ar[:])
      RangeAxis(axname,axdata)
    end
end

"Test if data in x can be approximated by a step range"
function testrange(x)
  r = range(first(x),last(x),length=length(x))
  all(i->isapprox(i...),zip(x,r)) ? r : x
end
import DataStructures: counter

function Cube(g::ZGroup;varlist=nothing,joinname="Variable")

  if varlist===nothing
    varlist = infervarlist(g)
  end
  v1 = g[varlist[1]]
  s = size(v1)
  vardims = reverse((v1.attrs["_ARRAY_DIMENSIONS"]...,))
  inneraxes = toaxis.(vardims,Ref(g))
  iax = collect(CubeAxis,inneraxes)
  s == length.(inneraxes) || throw(DimensionMismatch("Array dimensions do not fit"))
  allcubes = map(varlist) do iv
    v = g[iv]
    size(v) == s || throw(DimensionMismatch("All variables must have the same shape. $iv does not match $(varlist[1])"))
    ZArrayCube(v,iax,nothing)
  end
  # Filter out minority element types
  c = counter(eltype(i) for i in allcubes)
  _,et = findmax(c)
  indtake = findall(i->eltype(i)==et,allcubes)
  allcubes = allcubes[indtake]
  varlist  = varlist[indtake]
  if length(allcubes)==1
    return allcubes[1]
  else
    return concatenateCubes(allcubes,CategoricalAxis(joinname,varlist))
  end
end

sorted(x,y) = x<y ? (x,y) : (y,x)

interpretsubset(subexpr::Union{CartesianIndices{1},LinearIndices{1}},ax) = subexpr.indices[1]
interpretsubset(subexpr::CartesianIndex{1},ax)   = subexpr.I[1]
interpretsubset(subexpr,ax)                      = axVal2Index(ax,subexpr)
interpretsubset(subexpr::NTuple{2,Any},ax)       = Colon()(sorted(axVal2Index(ax,subexpr[1]),axVal2Index(ax,subexpr[2]))...)
#interpretsubset(subexpr::AbstractVector,ax)      = axVal2Index.(Ref(ax),subexpr)

import ESDL.Cubes: subsetCube
using InteractiveUtils

axcopy(ax::RangeAxis,vals) = RangeAxis(axname(ax),vals)
axcopy(ax::CategoricalAxis,vals) = CategoricalAxis(axname(ax),vals)

function subsetCube(z::ZArrayCube;kwargs...)
  subs = isnothing(z.subset) ? collect(Any,map(Base.OneTo,size(z))) : collect(z.subset)
  newaxes = deepcopy(caxes(z))
  foreach(kwargs) do kw
    axdes,subexpr = kw
    axdes = string(axdes)
    iax = findAxis(axdes,caxes(z))
    if isnothing(iax)
      throw(ArgumentError("Axis $axdes not found in cube"))
    else
      oldax = newaxes[iax]
      subinds = interpretsubset(subexpr,oldax)
      subs2 = subs[iax][subinds]
      subs[iax] = subs2
      newaxes[iax] = axcopy(oldax,oldax.values[subinds])
    end
  end
  newaxes = filter(ax->length(ax)>1,newaxes) |> collect
  ZArrayCube(z.a,newaxes,ntuple(i->subs[i],length(subs)))
end
end # module
