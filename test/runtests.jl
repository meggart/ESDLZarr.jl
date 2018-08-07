using ESDLZarr
using ESDL
using Zarr
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

c1 = Cube("/BGI/scratch/DataCube/v1.0.0/low-res/")
d1 = getCubeData(c1,variable="air_temperature_2m",longitude=(30,31),latitude=(50,51),time=(Date("2002-01-01"),Date("2008-12-31")))


c2=Cube(zopen("/home/fgans/zarr/esdc-8d-0.25deg-1x720x1440-1.0.1_1_zarr"))
d2 = getCubeData(c2,variable="air_temperature_2m",longitude=(30,31),latitude=(50,51),time=(Date("2002-01-01"),Date("2008-12-31")))


a1 = readCubeData(d1)

a2 = readCubeData(d2)
