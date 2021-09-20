# geodesic

[![GoDoc](https://godoc.org/github.com/tidwall/geodesic?status.svg)](https://godoc.org/github.com/tidwall/geodesic)

`geodesic` is a Go package providing operations for performing accurate measurements of Earth. Includes a Go port of the geodesic routines from [GeographicLib](https://geographiclib.sourceforge.io) along with general-purpose spherical algorithms.

## Features

- Precise distance calculations (5-15 nanometer)
- Calculate the area and perimeter of a polygon
- Get the length of a polyline
- Direct and Inverse solutions
- Optional spherical formulas for efficency, such as Haversine.

## Installing

To start using `geodesic`, install Go and run `go get`:

```sh
$ go get -u github.com/tidwall/geodesic
```

This will retrieve the library.

## Using 

All operations are performed on an Ellipsoid object.

To create a custom ellipsoid:
- `geodesic.NewEllipsoid`: For highly accurate geodesic routines.
- `geodesic.NewSpherical`: For faster spherical routines such as Haversine.

```go
// The first argument is the radius and the second is the flattening.
myEllipsoid := geodesic.NewEllipsoid(6378137.0, 1.0/298.257223563)
```

Or, you can simply use the built-in `geodesic.WGS84` non-spherical ellipsoid, which conforms to the [WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System) standard.

```go
// The argument is the radius.
mySphere := geodesic.NewSpherical(6378137.0)
```

Or, you can use the `geodesic.Globe` spherical ellipsoid, which represents Earth as a terrestrial globe.

## Examples

Calculate distance between two points.

```go
var dist float64
geodesic.WGS84.Inverse(33.4911, -112.4223, 32.1189, -113.1123, &dist, nil, nil)
fmt.Printf("%f meters\n", dist)
// output: 
// 165330.214571 meters
```

Calculate initial azimuth from point A to point B.

```go
var azi float64
geodesic.WGS84.Inverse(33.4911, -112.4223, 32.1189, -113.1123, nil, &azi, nil)
fmt.Printf("%f degrees\n", azi)
// output:
// -156.803310 degrees
```

Calculate final azimuth from point A to point B.

```go
var azi float64
geodesic.WGS84.Inverse(33.4911, -112.4223, 32.1189, -113.1123, nil, nil, &azi)
fmt.Printf("%f degrees\n", azi)
// output:
// -157.177169 degrees
```

Calculate distance and azimuths at the same time.

```go
var dist, azi1, azi2 float64
geodesic.WGS84.Inverse(33.4911, -112.4223, 32.1189, -113.1123, &dist, &azi1, &azi2)
fmt.Printf("%f meters, %f degrees, %f degrees\n", dist,azi1,azi2)
// output:
// 165330.214571 meters, -156.803310 degrees, -157.177169 degrees
```

Calculate destination using an initial point, azimuth, and distance.

```go
var lat, lon float64
geodesic.WGS84.Direct(33.4911, -112.4223, -156.803310, 165330.214571, &lat, &lon, nil)
fmt.Printf("%f, %f\n", lat, lon)
// output: 
// 32.118900, -113.112300
```

Calculate the area and perimeter of a polygon.

```go
// Arizona
poly := WGS84.PolygonInit(false)
poly.AddPoint(36.99377, -109.050292)
poly.AddPoint(36.96744, -114.049072)
poly.AddPoint(36.26199, -114.016113)
poly.AddPoint(36.08462, -114.279785)
poly.AddPoint(36.11125, -114.730224)
poly.AddPoint(34.86790, -114.631347)
poly.AddPoint(34.47939, -114.367675)
poly.AddPoint(34.29806, -114.104003)
poly.AddPoint(33.89777, -114.532470)
poly.AddPoint(33.58716, -114.543457)
poly.AddPoint(33.35806, -114.708251)
poly.AddPoint(33.09154, -114.708251)
poly.AddPoint(32.87036, -114.444580)
poly.AddPoint(32.74108, -114.719238)
poly.AddPoint(32.50049, -114.818115)
poly.AddPoint(31.33487, -111.093749)
poly.AddPoint(31.35363, -109.050292)
poly.AddPoint(36.99377, -109.050292)
var area, perimeter float64
poly.Compute(false, false, &area, &perimeter)
fmt.Printf("%f area (km²), %f perimeter (meters)\n", area/1000000, perimeter)
// output:
// 294838.722804 area (km²), 2254910.021767 perimeter (meters)
```
