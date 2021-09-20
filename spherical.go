// API for the shperical routines in Go
//
// Copyright (c) Joshua Baker (2021) and licensed under the MIT License.
//
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Latitude/longitude spherical geodesy tools   (c) Chris Veness 2002-2019 */
/*                                                             MIT Licence */
/* www.movable-type.co.uk/scripts/latlong.html                             */
/* www.movable-type.co.uk/scripts/geodesy-library.html#latlon-spherical    */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

package geodesic

import "math"

func sphericalInverse(
	radius float64,
	lat1 float64, lon1 float64,
	lat2 float64, lon2 float64,
	s12 *float64, azi1 *float64, azi2 *float64,
) {
	if s12 != nil {
		*s12 = distance(radius, lat1, lon1, lat2, lon2)
	}
	if azi1 != nil {
		*azi1 = bearing(lat1, lon1, lat2, lon2)
	}
	if azi2 != nil {
		*azi2 = wrap180(bearing(lat2, lon2, lat1, lon1) + 180)
	}
}

func sphericalDirect(
	radius float64,
	lat1 float64, lon1 float64, azi1 float64,
	s12 float64,
	lat2 *float64, lon2 *float64, azi2 *float64,
) {
	la2, lo2 := destination(radius, lat1, lon1, s12, azi1)
	if lat2 != nil {
		*lat2 = la2
	}
	if lon2 != nil {
		*lon2 = lo2
	}
	if azi2 != nil {
		*azi2 = wrap180(bearing(la2, lo2, lat1, lon1) + 180)
	}
}

const radians = math.Pi / 180
const degrees = 180 / math.Pi

func destination(radius float64, lat1, lon1, meters, bearingDegrees float64) (lat2, lon2 float64) {
	// sinφ2 = sinφ1⋅cosδ + cosφ1⋅sinδ⋅cosθ
	// tanΔλ = sinθ⋅sinδ⋅cosφ1 / cosδ−sinφ1⋅sinφ2
	// see mathforum.org/library/drmath/view/52049.html for derivation
	δ := meters / radius
	θ := bearingDegrees * radians
	φ1 := lat1 * radians
	λ1 := lon1 * radians
	φ2 := math.Asin(math.Sin(φ1)*math.Cos(δ) +
		math.Cos(φ1)*math.Sin(δ)*math.Cos(θ))
	λ2 := λ1 + math.Atan2(math.Sin(θ)*math.Sin(δ)*math.Cos(φ1),
		math.Cos(δ)-math.Sin(φ1)*math.Sin(φ2))
	λ2 = math.Mod(λ2+3*math.Pi, 2*math.Pi) - math.Pi // normalise to -180..+180°
	return φ2 * degrees, λ2 * degrees
}

func distance(radius float64, lat1, lon1, lat2, lon2 float64) float64 {
	// haversine formula
	φ1 := lat1 * radians
	λ1 := lon1 * radians
	φ2 := lat2 * radians
	λ2 := lon2 * radians
	Δφ := φ2 - φ1
	Δλ := λ2 - λ1
	sΔφ2 := math.Sin(Δφ / 2)
	sΔλ2 := math.Sin(Δλ / 2)
	haver := sΔφ2*sΔφ2 + math.Cos(φ1)*math.Cos(φ2)*sΔλ2*sΔλ2
	return radius * 2 * math.Asin(math.Sqrt(haver))
}

func bearing(lat1, lon1, lat2, lon2 float64) float64 {
	// tanθ = sinΔλ⋅cosφ2 / cosφ1⋅sinφ2 − sinφ1⋅cosφ2⋅cosΔλ
	// see mathforum.org/library/drmath/view/55417.html for derivation
	φ1 := lat1 * radians
	φ2 := lat2 * radians
	Δλ := (lon2 - lon1) * radians
	y := math.Sin(Δλ) * math.Cos(φ2)
	x := math.Cos(φ1)*math.Sin(φ2) - math.Sin(φ1)*math.Cos(φ2)*math.Cos(Δλ)
	θ := math.Atan2(y, x)
	return wrap180(θ * degrees)
}

func wrap180(degs float64) float64 {
	if degs < -180 || degs > 180 {
		degs = math.Mod(degs, 360)
		if degs < -180 {
			degs += 360
		} else if degs > 180 {
			degs -= 360
		}
	}
	return degs
}
