// Copyright 2021 Joshua J Baker. All rights reserved.
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.

package geodesic

import "math"

const (
	pi                         = float64(math.Pi)
	degree                     = pi / 180
	epsilon                    = float64(7.)/3 - float64(4.)/3 - float64(1.)
	digits                     = 53
	geographicLibGeodesicOrder = 6
	nA1                        = geographicLibGeodesicOrder
	nC1                        = geographicLibGeodesicOrder
	nC1p                       = geographicLibGeodesicOrder
	nA2                        = geographicLibGeodesicOrder
	nC2                        = geographicLibGeodesicOrder
	nA3                        = geographicLibGeodesicOrder
	nA3x                       = nA3
	nC3                        = geographicLibGeodesicOrder
	nC3x                       = ((nC3 * (nC3 - 1)) / 2)
	nC4                        = geographicLibGeodesicOrder
	nC4x                       = ((nC4 * (nC4 + 1)) / 2)
	nC                         = (geographicLibGeodesicOrder + 1)
	tol0                       = epsilon
	realmin                    = math.SmallestNonzeroFloat64
)

const (
	geodNone          = 0                   /**< Calculate nothing */
	geodLatitude      = 1 << 7              /**< Calculate latitude */
	geodLongitude     = 1<<8 | 1<<3         /**< Calculate longitude */
	geodAzimuth       = 1 << 9              /**< Calculate azimuth */
	geosDistance      = 1<<10 | 1<<0        /**< Calculate distance */
	geodDistanceIn    = 1<<11 | 1<<0 | 1<<1 /**< Allow distance as input  */
	geodReducedLength = 1<<12 | 1<<0 | 1<<2 /**< Calculate reduced length */
	geodGeodesicScale = 1<<13 | 1<<0 | 1<<2 /**< Calculate geodesic scale */
	geodArea          = 1<<14 | 1<<4        /**< Calculate reduced length */
	geodAll           = 0x7F80 | 0x1F       /**< Calculate everything */
)

/**
 * flag values for the \e flags argument to geod_gendirect() and
 * geod_genposition()
 **********************************************************************/
const (
	geodNoFlags    = 0       /**< No flags */
	geodArcMode    = 1 << 0  /**< Position given in terms of arc distance */
	geodLongUnroll = 1 << 15 /**< Unroll the longitude */
)

const (
	capNone = 0
	capC1   = 1 << 0
	capC1p  = 1 << 1
	capC2   = 1 << 2
	capC3   = 1 << 3
	capC4   = 1 << 4
	capAll  = 0x1F
	outAll  = 0x7F80
)

var (
	tiny         = sqrt(realmin)
	tol1         = 200 * tol0
	tol2         = sqrt(tol0)
	xthresh      = 1000 * tol2
	maxit1  uint = 20
	maxit2  uint = maxit1 + digits + 10
	tolb         = tol0 * tol2
)

type geodGeodesic struct {
	a float64 /**< the equatorial radius */
	f float64 /**< the flattening */
	/**< @cond SKIP */
	f1, e2, ep2, n, b, c2, etol2 float64
	A3x                          [6]float64
	C3x                          [15]float64
	C4x                          [21]float64
	/**< @endcond */
}

func geodInverse(g *geodGeodesic,
	lat1 float64, lon1 float64,
	lat2 float64, lon2 float64,
	ps12 *float64, pazi1 *float64, pazi2 *float64) {
	geodGenInverse(g, lat1, lon1, lat2, lon2, ps12, pazi1, pazi2,
		nil, nil, nil, nil)
}

func geodGenInverse(g *geodGeodesic,
	lat1 float64, lon1 float64, lat2 float64, lon2 float64,
	ps12 *float64, pazi1 *float64, pazi2 *float64,
	pm12 *float64, pM12 *float64, pM21 *float64, pS12 *float64,
) float64 {
	var salp1, calp1, salp2, calp2 float64
	a12 := geodGenInverseInt(g, lat1, lon1, lat2, lon2, ps12,
		&salp1, &calp1, &salp2, &calp2,
		pm12, pM12, pM21, pS12)
	if pazi1 != nil {
		*pazi1 = atan2dx(salp1, calp1)
	}
	if pazi2 != nil {
		*pazi2 = atan2dx(salp2, calp2)
	}
	return a12
}

func atan2dx(y float64, x float64) float64 {
	/* In order to minimize round-off errors, this function rearranges the
	 * arguments so that result of atan2 is in the range [-pi/4, pi/4] before
	 * converting it to degrees and mapping the result to the correct
	 * quadrant. */
	var q = 0
	var ang float64
	if fabs(y) > fabs(x) {
		x, y = y, x
		q = 2
	}
	if x < 0 {
		x = -x
		q++
	}
	/* here x >= 0 and x >= abs(y), so angle is in [-pi/4, pi/4] */
	ang = atan2(y, x) / degree
	switch q {
	/* Note that atan2d(-0.0, 1.0) will return -0.  However, we expect that
	 * atan2d will not be called with y = -0.  If need be, include
	 *
	 *   case 0: ang = 0 + ang; break;
	 */
	case 1:
		var v float64
		if y >= 0 {
			v = 180
		} else {
			v = -180
		}
		ang = v - ang
	case 2:
		ang = 90 - ang
	case 3:
		ang = -90 + ang
	}
	return ang
}

func geodGenInverseInt(g *geodGeodesic,
	lat1 float64, lon1 float64, lat2 float64, lon2 float64,
	ps12 *float64,
	psalp1 *float64, pcalp1 *float64,
	psalp2 *float64, pcalp2 *float64,
	pm12 *float64, pM12 *float64, pM21 *float64,
	pS12 *float64,
) float64 {
	var s12, m12, M12, M21, S12 float64
	var lon12, lon12s float64
	var latsign, lonsign, swapp int
	var sbet1, cbet1, sbet2, cbet2, s12x, m12x float64
	var dn1, dn2, lam12, slam12, clam12 float64
	var a12, sig12, calp1, salp1, calp2, salp2 float64
	var Ca [nC]float64
	var meridian bool
	/* somg12 > 1 marks that it needs to be calculated */
	var omg12, somg12, comg12 float64
	somg12 = 2

	var outmask uint
	if ps12 != nil {
		outmask |= geosDistance
	}
	if pm12 != nil {
		outmask |= geodReducedLength
	}
	if pM12 != nil || pM21 != nil {
		outmask |= geodGeodesicScale
	}
	if pS12 != nil {
		outmask |= geodArea
	}
	outmask &= outAll
	/* Compute longitude difference (AngDiff does this carefully).  Result is
	* in [-180, 180] but -180 is only for west-going geodesics.  180 is for
	* east-going and meridional geodesics. */
	lon12 = angDiff(lon1, lon2, &lon12s)
	/* Make longitude difference positive. */
	if lon12 >= 0 {
		lonsign = 1
	} else {
		lonsign = -1
	}
	/* If very close to being on the same half-meridian, then make it so. */
	lon12 = float64(lonsign) * angRound(lon12)
	lon12s = angRound((180 - lon12) - float64(lonsign)*lon12s)
	lam12 = lon12 * degree
	if lon12 > 90 {
		sincosdx(lon12s, &slam12, &clam12)
		clam12 = -clam12
	} else {
		sincosdx(lon12, &slam12, &clam12)
	}

	/* If really close to the equator, treat as on equator. */
	lat1 = angRound(latFix(lat1))
	lat2 = angRound(latFix(lat2))
	/* Swap points so that point with higher (abs) latitude is point 1
	* If one latitude is a nan, then it becomes lat1. */
	if fabs(lat1) < fabs(lat2) {
		swapp = -1
	} else {
		swapp = 1
	}
	if swapp < 0 {
		lonsign *= -1
		lat1, lat2 = lat2, lat1
	}
	/* Make lat1 <= 0 */
	if lat1 < 0 {
		latsign = 1
	} else {
		latsign = -1
	}
	lat1 *= float64(latsign)
	lat2 *= float64(latsign)

	/* Now we have
	*
	*     0 <= lon12 <= 180
	*     -90 <= lat1 <= 0
	*     lat1 <= lat2 <= -lat1
	*
	* longsign, swapp, latsign register the transformation to bring the
	* coordinates to this canonical form.  In all cases, 1 means no change was
	* made.  We make these transformations so that there are few cases to
	* check, e.g., on verifying quadrants in atan2.  In addition, this
	* enforces some symmetries in the results returned. */

	sincosdx(lat1, &sbet1, &cbet1)
	sbet1 *= g.f1
	/* Ensure cbet1 = +epsilon at poles */
	norm2(&sbet1, &cbet1)
	cbet1 = maxx(tiny, cbet1)

	sincosdx(lat2, &sbet2, &cbet2)
	sbet2 *= g.f1
	/* Ensure cbet2 = +epsilon at poles */

	norm2(&sbet2, &cbet2)
	cbet2 = maxx(tiny, cbet2)

	/* If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
	* |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
	* a better measure.  This logic is used in assigning calp2 in Lambda12.
	* Sometimes these quantities vanish and in that case we force bet2 = +/-
	* bet1 exactly.  An example where is is necessary is the inverse problem
	* 48.522876735459 0 -48.52287673545898293 179.599720456223079643
	* which failed with Visual Studio 10 (Release and Debug) */

	if cbet1 < -sbet1 {
		if cbet2 == cbet1 {
			if sbet2 < 0 {
				sbet2 = sbet1
			} else {
				sbet2 = -sbet1
			}
		}
	} else {
		if fabs(sbet2) == -sbet1 {
			cbet2 = cbet1
		}
	}

	dn1 = sqrt(1 + g.ep2*sq(sbet1))
	dn2 = sqrt(1 + g.ep2*sq(sbet2))

	meridian = lat1 == -90 || slam12 == 0

	if meridian {

		/* Endpoints are on a single full meridian, so the geodesic might lie on
		* a meridian. */

		var ssig1, csig1, ssig2, csig2 float64
		calp1 = clam12
		salp1 = slam12 /* Head to the target longitude */
		calp2 = 1
		salp2 = 0 /* At the target we're heading north */

		/* tan(bet) = tan(sig) * cos(alp) */
		ssig1 = sbet1
		csig1 = calp1 * cbet1
		ssig2 = sbet2
		csig2 = calp2 * cbet2

		/* sig12 = sig2 - sig1 */
		sig12 = atan2(maxx(0, csig1*ssig2-ssig1*csig2),
			csig1*csig2+ssig1*ssig2)
		var a, b *float64
		if outmask&geodGeodesicScale != 0 {
			a = &M12
			b = &M21
		}
		lengths(g, g.n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
			cbet1, cbet2, &s12x, &m12x, nil, a, b, Ca[:])
		/* Add the check for sig12 since zero length geodesics might yield m12 <
		* 0.  Test case was
		*
		*    echo 20.001 0 20.001 0 | GeodSolve -i
		*
		* In fact, we will have sig12 > pi/2 for meridional geodesic which is
		* not a shortest path. */
		if sig12 < 1 || m12x >= 0 {
			/* Need at least 2, to handle 90 0 90 180 */
			if sig12 < 3*tiny ||
				// Prevent negative s12 or m12 for short lines
				(sig12 < tol0 && (s12x < 0 || m12x < 0)) {
				sig12, m12x, s12x = 0, 0, 0
			}
			m12x *= g.b
			s12x *= g.b
			a12 = sig12 / degree
		} else {
			/* m12 < 0, i.e., prolate and too close to anti-podal */
			meridian = false
		}
	}

	if !meridian &&
		sbet1 == 0 && /* and sbet2 == 0 */
		/* Mimic the way Lambda12 works with calp1 = 0 */
		(g.f <= 0 || lon12s >= g.f*180) {
		/* Geodesic runs along equator */
		calp1, calp2 = 0, 0
		salp1, salp2 = 1, 1
		s12x = g.a * lam12
		omg12 = lam12 / g.f1
		sig12 = omg12
		m12x = g.b * sin(sig12)
		if outmask&geodGeodesicScale != 0 {
			M21 = cos(sig12)
			M12 = M21
		}
		a12 = lon12 / g.f1
	} else if !meridian {

		/* Now point1 and point2 belong within a hemisphere bounded by a
		* meridian and geodesic is neither meridional or equatorial. */

		/* Figure a starting point for Newton's method */
		var dnm float64

		sig12 = inverseStart(g, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
			lam12, slam12, clam12,
			&salp1, &calp1, &salp2, &calp2, &dnm,
			Ca[:])

		if sig12 >= 0 {
			/* Short lines (InverseStart sets salp2, calp2, dnm) */
			s12x = sig12 * g.b * dnm
			m12x = sq(dnm) * g.b * sin(sig12/dnm)
			if outmask&geodGeodesicScale != 0 {
				M21 = cos(sig12 / dnm)
				M12 = M21
			}
			a12 = sig12 / degree
			omg12 = lam12 / (g.f1 * dnm)
		} else {
			/* Newton's method.  This is a straightforward solution of f(alp1) =
			* lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
			* root in the interval (0, pi) and its derivative is positive at the
			* root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
			* alp1.  During the course of the iteration, a range (alp1a, alp1b) is
			* maintained which brackets the root and with each evaluation of
			* f(alp) the range is shrunk, if possible.  Newton's method is
			* restarted whenever the derivative of f is negative (because the new
			* value of alp1 is then further from the solution) or if the new
			* estimate of alp1 lies outside (0,pi); in this case, the new starting
			* guess is taken to be (alp1a + alp1b) / 2. */
			var ssig1, csig1, ssig2, csig2, eps, domg12 float64
			var numit uint
			/* Bracketing range */
			var salp1a, calp1a, salp1b, calp1b = tiny, 1.0, tiny, -1.0
			var tripn bool
			var tripb bool

			for ; numit < maxit2; numit++ {
				/* the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
				* WGS84 and random input: mean = 2.85, sd = 0.60 */
				var dv float64
				v := lambda12(g, sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
					slam12, clam12,
					&salp2, &calp2, &sig12, &ssig1, &csig1, &ssig2, &csig2,
					&eps, &domg12, numit < maxit1, &dv, Ca[:])

				/* Reversed test to allow escape with NaNs */
				var vv float64

				if tripn {
					vv = 8
				} else {
					vv = 1
				}
				if tripb || !(fabs(v) >= vv*tol0) {
					break
				}
				/* Update bracketing values */
				if v > 0 && (numit > maxit1 || calp1/salp1 > calp1b/salp1b) {
					salp1b = salp1
					calp1b = calp1
				} else if v < 0 && (numit > maxit1 || calp1/salp1 < calp1a/salp1a) {
					salp1a = salp1
					calp1a = calp1
				}
				if numit < maxit1 && dv > 0 {
					dalp1 := -v / dv
					sdalp1, cdalp1 := sincos(dalp1)
					nsalp1 := salp1*cdalp1 + calp1*sdalp1
					if nsalp1 > 0 && fabs(dalp1) < pi {
						calp1 = calp1*cdalp1 - salp1*sdalp1
						salp1 = nsalp1
						norm2(&salp1, &calp1)
						/* In some regimes we don't get quadratic convergence because
						* slope -> 0.  So use convergence conditions based on epsilon
						* instead of sqrt(epsilon). */
						tripn = fabs(v) <= 16*tol0
						continue
					}
				}
				/* Either dv was not positive or updated value was outside legal
				* range.  Use the midpoint of the bracket as the next estimate.
				* This mechanism is not needed for the WGS84 ellipsoid, but it does
				* catch problems with more eccentric ellipsoids.  Its efficacy is
				* such for the WGS84 test set with the starting guess set to alp1 =
				* 90deg:
				* the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
				* WGS84 and random input: mean = 4.74, sd = 0.99 */
				salp1 = (salp1a + salp1b) / 2
				calp1 = (calp1a + calp1b) / 2
				norm2(&salp1, &calp1)
				tripn = false
				tripb = (fabs(salp1a-salp1)+(calp1a-calp1) < tolb ||
					fabs(salp1-salp1b)+(calp1-calp1b) < tolb)
			}

			var v1 *float64
			var v2 *float64
			if outmask&geodGeodesicScale != 0 {
				v1 = &M12
				v2 = &M21
			}

			lengths(g, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
				cbet1, cbet2, &s12x, &m12x, nil,
				v1,
				v2, Ca[:])
			m12x *= g.b
			s12x *= g.b
			a12 = sig12 / degree

			if outmask&geodArea != 0 {
				/* omg12 = lam12 - domg12 */
				sdomg12, cdomg12 := sincos(domg12)
				somg12 = slam12*cdomg12 - clam12*sdomg12
				comg12 = clam12*cdomg12 + slam12*sdomg12
			}
		}
	}

	if outmask&geosDistance != 0 {
		s12 = 0 + s12x /* Convert -0 to 0 */
	}

	if outmask&geodReducedLength != 0 {
		m12 = 0 + m12x /* Convert -0 to 0 */
	}

	if outmask&geodArea != 0 {
		/* From Lambda12: sin(alp1) * cos(bet1) = sin(alp0) */
		salp0 := salp1 * cbet1
		calp0 := hypot(calp1, salp1*sbet1) /* calp0 > 0 */
		var alp12 float64
		if calp0 != 0 && salp0 != 0 {
			/* From Lambda12: tan(bet) = tan(sig) * cos(alp) */
			ssig1 := sbet1
			csig1 := calp1 * cbet1
			ssig2 := sbet2
			csig2 := calp2 * cbet2
			k2 := sq(calp0) * g.ep2
			eps := k2 / (2*(1+sqrt(1+k2)) + k2)
			/* Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0). */
			A4 := sq(g.a) * calp0 * salp0 * g.e2
			var B41, B42 float64
			norm2(&ssig1, &csig1)
			norm2(&ssig2, &csig2)
			c4f(g, eps, Ca[:])
			B41 = sinCosSeries(iFalse, ssig1, csig1, Ca[:], nC4)
			B42 = sinCosSeries(iFalse, ssig2, csig2, Ca[:], nC4)
			S12 = A4 * (B42 - B41)
		} else {
			/* Avoid problems with indeterminate sig1, sig2 on equator */
			S12 = 0
		}
		if !meridian && somg12 > 1 {
			somg12, comg12 = sincos(omg12)
		}

		if !meridian &&
			/* omg12 < 3/4 * pi */
			comg12 > -(real)(0.7071) && /* Long difference not too big */
			sbet2-sbet1 < (real)(1.75) { /* Lat difference not too big */
			/* Use tan(Gamma/2) = tan(omg12/2)
			* * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
			* with tan(x/2) = sin(x)/(1+cos(x)) */
			domg12 := 1 + comg12
			dbet1 := 1 + cbet1
			dbet2 := 1 + cbet2
			alp12 = 2 * atan2(somg12*(sbet1*dbet2+sbet2*dbet1),
				domg12*(sbet1*sbet2+dbet1*dbet2))
		} else {
			/* alp12 = alp2 - alp1, used in atan2 so no need to normalize */
			salp12 := salp2*calp1 - calp2*salp1
			calp12 := calp2*calp1 + salp2*salp1
			/* The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
			* salp12 = -0 and alp12 = -180.  However this depends on the sign
			* being attached to 0 correctly.  The following ensures the correct
			* behavior. */
			if salp12 == 0 && calp12 < 0 {
				salp12 = tiny * calp1
				calp12 = -1
			}
			alp12 = atan2(salp12, calp12)
		}
		S12 += g.c2 * alp12
		S12 *= float64(swapp * lonsign * latsign)
		/* Convert -0 to 0 */
		S12 += 0
	}

	/* Convert calp, salp to azimuth accounting for lonsign, swapp, latsign. */
	if swapp < 0 {
		salp1, salp2 = salp2, salp1
		calp1, calp2 = calp2, calp1
		if outmask&geodGeodesicScale != 0 {
			M12, M21 = M21, M12
		}
	}

	salp1 *= float64(swapp * lonsign)
	calp1 *= float64(swapp * latsign)
	salp2 *= float64(swapp * lonsign)
	calp2 *= float64(swapp * latsign)

	if psalp1 != nil {
		*psalp1 = salp1
	}
	if pcalp1 != nil {
		*pcalp1 = calp1
	}
	if psalp2 != nil {
		*psalp2 = salp2
	}
	if pcalp2 != nil {
		*pcalp2 = calp2
	}

	if outmask&geosDistance != 0 {
		*ps12 = s12
	}
	if outmask&geodReducedLength != 0 {
		*pm12 = m12
	}
	if outmask&geodGeodesicScale != 0 {
		if pM12 != nil {
			*pM12 = M12
		}
		if pM21 != nil {
			*pM21 = M21
		}
	}
	if outmask&geodArea != 0 {
		*pS12 = S12
	}

	/* Returned value in [0, 180] */
	return a12
}

func angDiff(x float64, y float64, e *float64) float64 {
	var t float64
	var d = angNormalize(sumx(angNormalize(-x), angNormalize(y), &t))
	/* Here y - x = d + t (mod 360), exactly, where d is in (-180,180] and
	 * abs(t) <= eps (eps = 2^-45 for doubles).  The only case where the
	 * addition of t takes the result outside the range (-180,180] is d = 180
	 * and t > 0.  The case, d = -180 + eps, t = -eps, can't happen, since
	 * sum would have returned the exact result in such a case (i.e., given t
	 * = 0). */
	var v float64
	if d == 180 && t > 0 {
		v = -180
	} else {
		v = d
	}
	return sumx(v, t, e)
}

func sumx(u float64, v float64, t *float64) float64 {
	var s = u + v
	var up = s - v
	var vpp = s - up
	up -= u
	vpp -= v
	if t != nil {
		*t = -(up + vpp)
	}
	/* error-free sum:
	 * u + v =       s      + t
	 *       = round(u + v) + t */
	return s
}

func angNormalize(x float64) float64 {
	x = remainder(x, 360)
	if x != -180 {
		return x
	}
	return 180
}

func angRound(x float64) float64 {
	const z = 1.0 / 16.0
	var y float64
	if x == 0 {
		return 0
	}
	y = fabs(x)
	/* The compiler mustn't "simplify" z - (z - y) to y */
	if y < z {
		y = z - (z - y)
	}
	if x < 0 {
		return -y
	}
	return y
}

func remquo(x float64, y float64, q *float64) float64 {
	*q = x / y
	return mod(x, y)
}

func sincosdx(x float64, sinx *float64, cosx *float64) {
	/* In order to minimize round-off errors, this function exactly reduces
	 * the argument to the range [-45, 45] before converting it to radians. */
	var r, s, c float64
	var q float64
	r = remquo(x, 90, &q)
	/* now abs(r) <= 45 */
	r *= degree
	/* Possibly could call the gnu extension sincos */
	s, c = sincos(r)
	switch uint64(int64(q)) & 3 {
	case 0:
		*sinx = s
		*cosx = c
	case 1:
		*sinx = c
		*cosx = -s
	case 2:
		*sinx = -s
		*cosx = -c
	default:
		*sinx = -c
		*cosx = s
		/* case 3U */
	}
	if x != 0 {
		*sinx += 0
		*cosx += 0
	}
}

func latFix(x float64) float64 {
	if fabs(x) > 90 {
		return math.NaN()
	}
	return x
}

func norm2(sinx *float64, cosx *float64) {
	r := hypot(*sinx, *cosx)
	*sinx /= r
	*cosx /= r
}
func sq(x float64) float64 { return x * x }

func lengths(g *geodGeodesic,
	eps float64, sig12 float64,
	ssig1 float64, csig1 float64, dn1 float64,
	ssig2 float64, csig2 float64, dn2 float64,
	cbet1 float64, cbet2 float64,
	ps12b *float64, pm12b *float64, pm0 *float64,
	pM12 *float64, pM21 *float64,
	/* Scratch area of the right size */
	Ca []float64,
) {
	var m0, J12, A1, A2 float64
	var Cb [nC]float64

	/* Return m12b = (reduced length)/b; also calculate s12b = distance/b,
	* and m0 = coefficient of secular term in expression for reduced length. */
	redlp := pm12b != nil || pm0 != nil || pM12 != nil || pM21 != nil
	if ps12b != nil || redlp {
		A1 = a1m1f(eps)
		c1f(eps, Ca)
		if redlp {
			A2 = a2m1f(eps)
			c2f(eps, Cb[:])
			m0 = A1 - A2
			A2 = 1 + A2
		}
		A1 = 1 + A1
	}
	if ps12b != nil {
		B1 := sinCosSeries(iTrue, ssig2, csig2, Ca, nC1) -
			sinCosSeries(iTrue, ssig1, csig1, Ca, nC1)
		/* Missing a factor of b */
		*ps12b = A1 * (sig12 + B1)
		if redlp {
			B2 := sinCosSeries(iTrue, ssig2, csig2, Cb[:], nC2) -
				sinCosSeries(iTrue, ssig1, csig1, Cb[:], nC2)
			J12 = m0*sig12 + (A1*B1 - A2*B2)
		}
	} else if redlp {
		/* Assume here that nC1 >= nC2 */
		var l int
		for l = 1; l <= nC2; l++ {
			Cb[l] = A1*Ca[l] - A2*Cb[l]
		}
		J12 = m0*sig12 + (sinCosSeries(iTrue, ssig2, csig2, Cb[:], nC2) -
			sinCosSeries(iTrue, ssig1, csig1, Cb[:], nC2))
	}
	if pm0 != nil {
		*pm0 = m0
	}
	if pm12b != nil {
		/* Missing a factor of b.
		* Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
		* accurate cancellation in the case of coincident points. */
		*pm12b = dn2*(csig1*ssig2) - dn1*(ssig1*csig2) -
			csig1*csig2*J12
	}
	if pM12 != nil || pM21 != nil {
		csig12 := csig1*csig2 + ssig1*ssig2
		t := g.ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
		if pM12 != nil {
			*pM12 = csig12 + (t*ssig2-csig2*J12)*ssig1/dn1
		}
		if pM21 != nil {
			*pM21 = csig12 - (t*ssig1-csig1*J12)*ssig2/dn2
		}
	}
}

/* The scale factor A1-1 = mean value of (d/dsigma)I1 - 1 */
func a1m1f(eps float64) float64 {
	var coeff = [...]float64{
		/* (1-eps)*A1-1, polynomial in eps2 of order 3 */
		1, 4, 64, 0, 256,
	}
	m := nA1 / 2
	t := polyval(m, coeff[:], sq(eps)) / coeff[m+1]
	return (t + eps) / (1 - eps)
}

func polyval(N int, p []float64, x float64) (y float64) {
	for i := 0; i <= N; i++ {
		y = y*x + p[i]
	}
	return y
}

var coeffC1f = [...]float64{
	/* C1[1]/eps^1, polynomial in eps2 of order 2 */
	-1, 6, -16, 32,
	/* C1[2]/eps^2, polynomial in eps2 of order 2 */
	-9, 64, -128, 2048,
	/* C1[3]/eps^3, polynomial in eps2 of order 1 */
	9, -16, 768,
	/* C1[4]/eps^4, polynomial in eps2 of order 1 */
	3, -5, 512,
	/* C1[5]/eps^5, polynomial in eps2 of order 0 */
	-7, 1280,
	/* C1[6]/eps^6, polynomial in eps2 of order 0 */
	-7, 2048,
}

/* The coefficients C1[l] in the Fourier expansion of B1 */
func c1f(eps float64, c []float64) {
	var eps2 = sq(eps)
	var d = eps
	var o int
	for l := 1; l <= nC1; l++ { /* l is index of C1p[l] */
		m := (nC1 - l) / 2 /* order of polynomial in eps^2 */
		c[l] = d * polyval(m, coeffC1f[o:], eps2) / coeffC1f[o+m+1]
		o += m + 2
		d *= eps
	}
}

/* The scale factor A2-1 = mean value of (d/dsigma)I2 - 1 */
func a2m1f(eps float64) float64 {
	var coeff = [...]float64{
		/* (eps+1)*A2-1, polynomial in eps2 of order 3 */
		-11, -28, -192, 0, 256,
	}
	m := nA2 / 2
	t := polyval(m, coeff[:], sq(eps)) / coeff[m+1]
	return (t - eps) / (1 + eps)
}

var coeffC2f = [...]float64{
	/* C2[1]/eps^1, polynomial in eps2 of order 2 */
	1, 2, 16, 32,
	/* C2[2]/eps^2, polynomial in eps2 of order 2 */
	35, 64, 384, 2048,
	/* C2[3]/eps^3, polynomial in eps2 of order 1 */
	15, 80, 768,
	/* C2[4]/eps^4, polynomial in eps2 of order 1 */
	7, 35, 512,
	/* C2[5]/eps^5, polynomial in eps2 of order 0 */
	63, 1280,
	/* C2[6]/eps^6, polynomial in eps2 of order 0 */
	77, 2048,
}

/* The coefficients C2[l] in the Fourier expansion of B2 */
func c2f(eps float64, c []float64) {
	var eps2 = eps * eps //sq(eps)
	var d = eps
	var o int
	for l := 1; l <= nC2; l++ { /* l is index of C2[l] */
		m := (nC2 - l) / 2 /* order of polynomial in eps^2 */
		c[l] = d * polyval(m, coeffC2f[o:], eps2) / coeffC2f[o+m+1]
		o += m + 2
		d *= eps
	}
}

type ibool int

const (
	iFalse ibool = 0
	iTrue  ibool = 1
)

func sinCosSeries(sinpb ibool, sinx, cosx float64, c []float64, n int) float64 {
	/* Evaluate
	 * y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
	 *            sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
	 * using Clenshaw summation.  N.B. c[0] is unused for sin series
	 * Approx operation count = (n + 5) mult and (2 * n + 2) add */
	ci := n + int(sinpb)                    /* Point to one beyond last element */
	ar := 2 * (cosx - sinx) * (cosx + sinx) /* 2 * cos(2 * x) */
	y0 := 0.0
	if n&1 != 0 {
		ci--
		y0 = c[0]
	}
	y1 := 0.0 /* accumulators for sum */
	/* Now n is even */
	n /= 2
	for n != 0 {
		n--
		/* Unroll loop x 2, so accumulators return to their original role */
		ci--
		y1 = ar*y0 - y1 + c[ci]
		ci--
		y0 = ar*y1 - y0 + c[ci]

	}
	if sinpb == iTrue {
		return 2 * sinx * cosx * y0 /* sin(2 * x) * y0 */
	}
	return cosx * (y0 - y1) /* cos(x) * (y0 - y1) */
}

func inverseStart(g *geodGeodesic,
	sbet1 float64, cbet1 float64, dn1 float64,
	sbet2 float64, cbet2 float64, dn2 float64,
	lam12 float64, slam12 float64, clam12 float64,
	psalp1 *float64, pcalp1 *float64,
	/* Only updated if return val >= 0 */
	psalp2 *float64, pcalp2 *float64,
	/* Only updated for short lines */
	pdnm *float64,
	/* Scratch area of the right size */
	Ca []float64,
) float64 {
	var salp1, calp1, salp2, calp2, dnm float64

	/* Return a starting point for Newton's method in salp1 and calp1 (function
	* value is -1).  If Newton's method doesn't need to be used, return also
	* salp2 and calp2 and function value is sig12. */
	var sig12 = -1.0 /* Return value */
	/* bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0] */
	var sbet12 = sbet2*cbet1 - cbet2*sbet1
	var cbet12 = cbet2*cbet1 + sbet2*sbet1
	var sbet12a float64
	var shortline = cbet12 >= 0 && sbet12 < 0.5 &&
		cbet2*lam12 < 0.5
	var somg12, comg12, ssig12, csig12 float64
	sbet12a = sbet2*cbet1 + cbet2*sbet1
	if shortline {
		var sbetm2 = sq(sbet1 + sbet2)
		var omg12 float64
		/* sin((bet1+bet2)/2)^2
		* =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2) */
		sbetm2 /= sbetm2 + sq(cbet1+cbet2)
		dnm = sqrt(1 + g.ep2*sbetm2)
		omg12 = lam12 / (g.f1 * dnm)
		somg12, comg12 = sincos(omg12)
	} else {
		somg12 = slam12
		comg12 = clam12
	}

	salp1 = cbet2 * somg12
	if comg12 >= 0 {
		calp1 = sbet12 + cbet2*sbet1*sq(somg12)/(1+comg12)
	} else {
		calp1 = sbet12a - cbet2*sbet1*sq(somg12)/(1-comg12)
	}

	ssig12 = hypot(salp1, calp1)
	csig12 = sbet1*sbet2 + cbet1*cbet2*comg12

	if shortline && ssig12 < g.etol2 {
		// /* really short lines */
		salp2 = cbet1 * somg12
		var v float64
		if comg12 >= 0 {
			v = sq(somg12) / (1 + comg12)
		} else {
			v = 1 - comg12
		}
		calp2 = sbet12 - cbet1*sbet2*v
		norm2(&salp2, &calp2)
		/* Set return value */
		sig12 = atan2(ssig12, csig12)
	} else if fabs(g.n) > 0.1 || /* No astroid calc if too eccentric */
		csig12 >= 0 ||
		ssig12 >= 6*fabs(g.n)*pi*sq(cbet1) {
		/* Nothing to do, zeroth order spherical approximation is OK */
	} else {
		/* Scale lam12 and bet2 to x, y coordinate system where antipodal point
		* is at origin and singular point is at y = 0, x = -1. */
		var y, lamscale, betscale float64
		/* Volatile declaration needed to fix inverse case
		* 56.320923501171 0 -56.320923501171 179.664747671772880215
		* which otherwise fails with g++ 4.4.4 x86 -O3 */
		var x float64
		lam12x := atan2(-slam12, -clam12) /* lam12 - pi */
		if g.f >= 0 {                     /* In fact f == 0 does not get here */
			/* x = dlong, y = dlat */
			{

				k2 := sq(sbet1) * g.ep2
				eps := k2 / (2*(1+sqrt(1+k2)) + k2)
				lamscale = g.f * cbet1 * a3f(g, eps) * pi
			}
			betscale = lamscale * cbet1

			x = lam12x / lamscale
			y = sbet12a / betscale
		} else { /* f < 0 */
			/* x = dlat, y = dlong */

			cbet12a := cbet2*cbet1 - sbet2*sbet1
			bet12a := atan2(sbet12a, cbet12a)
			var m12b, m0 float64
			/* In the case of lon12 = 180, this repeats a calculation made in
			* Inverse. */
			lengths(g, g.n, pi+bet12a,
				sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
				cbet1, cbet2, nil, &m12b, &m0, nil, nil, Ca)
			x = -1 + m12b/(cbet1*cbet2*m0*pi)
			if x < -0.01 {
				betscale = sbet12a / x
			} else {
				betscale = -g.f * sq(cbet1) * pi
			}
			lamscale = betscale / cbet1
			y = lam12x / lamscale
		}

		if y > -tol1 && x > -1-xthresh {
			/* strip near cut */
			if g.f >= 0 {
				salp1 = minx(1, -x)
				calp1 = -sqrt(1 - sq(salp1))
			} else {
				var v float64
				if x > -tol1 {
					v = 0
				} else {
					v = -1
				}
				calp1 = maxx(v, x)
				salp1 = sqrt(1 - sq(calp1))
			}
		} else {
			/* Estimate alp1, by solving the astroid problem.
			*
			* Could estimate alpha1 = theta + pi/2, directly, i.e.,
			*   calp1 = y/k; salp1 = -x/(1+k);  for f >= 0
			*   calp1 = x/(1+k); salp1 = -y/k;  for f < 0 (need to check)
			*
			* However, it's better to estimate omg12 from astroid and use
			* spherical formula to compute alp1.  This reduces the mean number of
			* Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
			* (min 0 max 5).  The changes in the number of iterations are as
			* follows:
			*
			* change percent
			*    1       5
			*    0      78
			*   -1      16
			*   -2       0.6
			*   -3       0.04
			*   -4       0.002
			*
			* The histogram of iterations is (m = number of iterations estimating
			* alp1 directly, n = number of iterations estimating via omg12, total
			* number of trials = 148605):
			*
			*  iter    m      n
			*    0   148    186
			*    1 13046  13845
			*    2 93315 102225
			*    3 36189  32341
			*    4  5396      7
			*    5   455      1
			*    6    56      0
			*
			* Because omg12 is near pi, estimate work with omg12a = pi - omg12 */
			k := astroid(x, y)
			var v float64
			if g.f >= 0 {
				v = -x * k / (1 + k)
			} else {
				v = -y * (1 + k) / k
			}
			omg12a := lamscale * v
			somg12, comg12 = sincos(omg12a)
			comg12 = -comg12
			/* Update spherical estimate of alp1 using omg12 instead of lam12 */
			salp1 = cbet2 * somg12
			calp1 = sbet12a - cbet2*sbet1*sq(somg12)/(1-comg12)
		}
	}
	/* Sanity check on starting guess.  Backwards check allows NaN through. */
	if !(salp1 <= 0) {
		norm2(&salp1, &calp1)
	} else {
		salp1 = 1
		calp1 = 0
	}

	*psalp1 = salp1
	*pcalp1 = calp1
	if shortline {
		*pdnm = dnm
	}
	if sig12 >= 0 {
		*psalp2 = salp2
		*pcalp2 = calp2
	}
	return sig12
}

func geodInit(g *geodGeodesic, a float64, f float64) {
	g.a = a
	g.f = f
	g.f1 = 1 - g.f
	g.e2 = g.f * (2 - g.f)
	g.ep2 = g.e2 / sq(g.f1) /* e2 / (1 - e2) */
	g.n = g.f / (2 - g.f)
	g.b = g.a * g.f1
	var v float64
	if g.e2 == 0 {
		v = 1
	} else if g.e2 > 0 {
		v = atanh(sqrt(g.e2))
	} else {
		v = atan(sqrt(-g.e2))
	}
	g.c2 = (sq(g.a) + sq(g.b)*(v/sqrt(fabs(g.e2)))) / 2 /* authalic radius squared */

	/* The sig12 threshold for "really short".  Using the auxiliary sphere
	 * solution with dnm computed at (bet1 + bet2) / 2, the relative error in the
	 * azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.  (Error
	 * measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a given f and
	 * sig12, the max error occurs for lines near the pole.  If the old rule for
	 * computing dnm = (dn1 + dn2)/2 is used, then the error increases by a
	 * factor of 2.)  Setting this equal to epsilon gives sig12 = etol2.  Here
	 * 0.1 is a safety factor (error decreased by 100) and max(0.001, abs(f))
	 * stops etol2 getting too large in the nearly spherical case. */
	g.etol2 = 0.1 * tol2 /
		sqrt(maxx(0.001, fabs(g.f))*minx(1, 1-g.f/2)/2)

	a3coeff(g)
	c3coeff(g)
	c4coeff(g)
}

/* The scale factor A3 = mean value of (d/dsigma)I3 */
func a3coeff(g *geodGeodesic) {
	var coeff = [...]float64{
		/* A3, coeff of eps^5, polynomial in n of order 0 */
		-3, 128,
		/* A3, coeff of eps^4, polynomial in n of order 1 */
		-2, -3, 64,
		/* A3, coeff of eps^3, polynomial in n of order 2 */
		-1, -3, -1, 16,
		/* A3, coeff of eps^2, polynomial in n of order 2 */
		3, -1, -2, 8,
		/* A3, coeff of eps^1, polynomial in n of order 1 */
		1, -1, 2,
		/* A3, coeff of eps^0, polynomial in n of order 0 */
		1, 1,
	}
	var o, k, j int
	for j = nA3 - 1; j >= 0; j-- { /* coeff of eps^j */
		var m int
		if nA3-j-1 < j {
			m = nA3 - j - 1
		} else {
			m = j
		} /* order of polynomial in n */
		g.A3x[k] = polyval(m, coeff[o:], g.n) / coeff[o+m+1]
		k++
		o += m + 2
	}
}

/* The coefficients C3[l] in the Fourier expansion of B3 */
func c3coeff(g *geodGeodesic) {
	var coeff = [...]float64{
		/* C3[1], coeff of eps^5, polynomial in n of order 0 */
		3, 128,
		/* C3[1], coeff of eps^4, polynomial in n of order 1 */
		2, 5, 128,
		/* C3[1], coeff of eps^3, polynomial in n of order 2 */
		-1, 3, 3, 64,
		/* C3[1], coeff of eps^2, polynomial in n of order 2 */
		-1, 0, 1, 8,
		/* C3[1], coeff of eps^1, polynomial in n of order 1 */
		-1, 1, 4,
		/* C3[2], coeff of eps^5, polynomial in n of order 0 */
		5, 256,
		/* C3[2], coeff of eps^4, polynomial in n of order 1 */
		1, 3, 128,
		/* C3[2], coeff of eps^3, polynomial in n of order 2 */
		-3, -2, 3, 64,
		/* C3[2], coeff of eps^2, polynomial in n of order 2 */
		1, -3, 2, 32,
		/* C3[3], coeff of eps^5, polynomial in n of order 0 */
		7, 512,
		/* C3[3], coeff of eps^4, polynomial in n of order 1 */
		-10, 9, 384,
		/* C3[3], coeff of eps^3, polynomial in n of order 2 */
		5, -9, 5, 192,
		/* C3[4], coeff of eps^5, polynomial in n of order 0 */
		7, 512,
		/* C3[4], coeff of eps^4, polynomial in n of order 1 */
		-14, 7, 512,
		/* C3[5], coeff of eps^5, polynomial in n of order 0 */
		21, 2560,
	}
	var o, k, l, j int
	for l = 1; l < nC3; l++ { /* l is index of C3[l] */
		for j = nC3 - 1; j >= l; j-- { /* coeff of eps^j */
			var m int
			if nC3-j-1 < j {
				m = nC3 - j - 1
			} else {
				m = j
			} /* order of polynomial in n */
			g.C3x[k] = polyval(m, coeff[o:], g.n) / coeff[o+m+1]
			k++
			o += m + 2
		}
	}
}

/* The coefficients C4[l] in the Fourier expansion of I4 */
func c4coeff(g *geodGeodesic) {
	var coeff = [...]float64{
		/* C4[0], coeff of eps^5, polynomial in n of order 0 */
		97, 15015,
		/* C4[0], coeff of eps^4, polynomial in n of order 1 */
		1088, 156, 45045,
		/* C4[0], coeff of eps^3, polynomial in n of order 2 */
		-224, -4784, 1573, 45045,
		/* C4[0], coeff of eps^2, polynomial in n of order 3 */
		-10656, 14144, -4576, -858, 45045,
		/* C4[0], coeff of eps^1, polynomial in n of order 4 */
		64, 624, -4576, 6864, -3003, 15015,
		/* C4[0], coeff of eps^0, polynomial in n of order 5 */
		100, 208, 572, 3432, -12012, 30030, 45045,
		/* C4[1], coeff of eps^5, polynomial in n of order 0 */
		1, 9009,
		/* C4[1], coeff of eps^4, polynomial in n of order 1 */
		-2944, 468, 135135,
		/* C4[1], coeff of eps^3, polynomial in n of order 2 */
		5792, 1040, -1287, 135135,
		/* C4[1], coeff of eps^2, polynomial in n of order 3 */
		5952, -11648, 9152, -2574, 135135,
		/* C4[1], coeff of eps^1, polynomial in n of order 4 */
		-64, -624, 4576, -6864, 3003, 135135,
		/* C4[2], coeff of eps^5, polynomial in n of order 0 */
		8, 10725,
		/* C4[2], coeff of eps^4, polynomial in n of order 1 */
		1856, -936, 225225,
		/* C4[2], coeff of eps^3, polynomial in n of order 2 */
		-8448, 4992, -1144, 225225,
		/* C4[2], coeff of eps^2, polynomial in n of order 3 */
		-1440, 4160, -4576, 1716, 225225,
		/* C4[3], coeff of eps^5, polynomial in n of order 0 */
		-136, 63063,
		/* C4[3], coeff of eps^4, polynomial in n of order 1 */
		1024, -208, 105105,
		/* C4[3], coeff of eps^3, polynomial in n of order 2 */
		3584, -3328, 1144, 315315,
		/* C4[4], coeff of eps^5, polynomial in n of order 0 */
		-128, 135135,
		/* C4[4], coeff of eps^4, polynomial in n of order 1 */
		-2560, 832, 405405,
		/* C4[5], coeff of eps^5, polynomial in n of order 0 */
		128, 99099,
	}
	var o, k, l, j int
	for l = 0; l < nC4; l++ { /* l is index of C4[l] */
		for j = nC4 - 1; j >= l; j-- { /* coeff of eps^j */
			m := nC4 - j - 1 /* order of polynomial in n */
			g.C4x[k] = polyval(m, coeff[o:], g.n) / coeff[o+m+1]
			k++
			o += m + 2
		}
	}
}

func astroid(x, y float64) float64 {
	/* Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
	 * This solution is adapted from Geocentric::Reverse. */
	var k float64
	p := sq(x)
	q := sq(y)
	r := (p + q - 1) / 6
	if !(q == 0 && r <= 0) {
		// 	/* Avoid possible division by zero when r = 0 by multiplying equations
		// 	 * for s and t by r^3 and r, resp. */
		S := p * q / 4 /* S = r^3 * s */
		r2 := sq(r)
		r3 := r * r2
		/* The discriminant of the quadratic equation for T3.  This is zero on
		 * the evolute curve p^(1/3)+q^(1/3) = 1 */
		disc := S * (S + 2*r3)
		u := r
		var v, uv, w float64
		if disc >= 0 {
			T3, T := S+r3, 0.0
			/* Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
			 * of precision due to cancellation.  The result is unchanged because
			 * of the way the T is used in definition of u. */
			if T3 < 0 {
				T3 += -sqrt(disc)
			} else {
				T3 += sqrt(disc)
			} /* T3 = (r * t)^3 */
			/* N.B. cbrt always returns the real root.  cbrt(-8) = -2. */
			T = cbrt(T3) /* T = r * t */
			/* T can be zero; but then r2 / T -> 0. */
			var v float64
			if T != 0 {
				v = r2 / T
			} else {
				v = 0
			}
			u += T + v
		} else {
			/* T is complex, but the way u is defined the result is real. */
			ang := atan2(sqrt(-disc), -(S + r3))
			/* There are three possible cube roots.  We choose the root which
			 * avoids cancellation.  Note that disc < 0 implies that r < 0. */
			u += 2 * r * cos(ang/3)
		}
		v = sqrt(sq(u) + q) /* guaranteed positive */
		/* Avoid loss of accuracy when u < 0. */
		if u < 0 {
			uv = q / (v - u)
		} else {
			uv = u + v
		} /* u+v, guaranteed positive */
		w = (uv - q) / (2 * v) /* positive? */
		/* Rearrange expression for k to avoid loss of accuracy due to
		 * subtraction.  Division by 0 not possible because uv > 0, w >= 0. */
		k = uv / (sqrt(uv+sq(w)) + w) /* guaranteed positive */
	} else { /* q == 0 && r <= 0 */
		/* y = 0 with |x| <= 1.  Handle this case directly.
		 * for y small, positive root is k = abs(y)/sqrt(1-x^2) */
		k = 0
	}
	return k
}

func a3f(g *geodGeodesic, eps float64) float64 {
	/* Evaluate A3 */
	return polyval(nA3-1, g.A3x[:], eps)
}

func lambda12(g *geodGeodesic,
	sbet1 float64, cbet1 float64, dn1 float64,
	sbet2 float64, cbet2 float64, dn2 float64,
	salp1 float64, calp1 float64,
	slam120 float64, clam120 float64,
	psalp2 *float64, pcalp2 *float64,
	psig12 *float64,
	pssig1 *float64, pcsig1 *float64,
	pssig2 *float64, pcsig2 *float64,
	peps *float64,
	pdomg12 *float64,
	diffp bool, pdlam12 *float64,
	/* Scratch area of the right size */
	Ca []float64,
) float64 {
	var salp2, calp2, sig12,
		ssig1, csig1, ssig2, csig2, eps,
		domg12, dlam12 float64
	var salp0, calp0 float64
	var somg1, comg1, somg2, comg2, somg12, comg12, lam12 float64
	var B312, eta, k2 float64

	if sbet1 == 0 && calp1 == 0 {
		/* Break degeneracy of equatorial line.  This case has already been
		* handled. */
		calp1 = -tiny
	}

	/* sin(alp1) * cos(bet1) = sin(alp0) */
	salp0 = salp1 * cbet1
	calp0 = hypot(calp1, salp1*sbet1) /* calp0 > 0 */

	/* tan(bet1) = tan(sig1) * cos(alp1)
	* tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1) */
	ssig1 = sbet1
	somg1 = salp0 * sbet1
	comg1 = calp1 * cbet1
	csig1 = comg1
	norm2(&ssig1, &csig1)
	/* norm2(&somg1, &comg1); -- don't need to normalize! */

	/* Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
	* about this case, since this can yield singularities in the Newton
	* iteration.
	* sin(alp2) * cos(bet2) = sin(alp0) */
	if cbet2 != cbet1 {
		salp2 = salp0 / cbet2
	} else {
		salp2 = salp1
	}

	/* calp2 = sqrt(1 - sq(salp2))
	*       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
	* and subst for calp0 and rearrange to give (choose positive sqrt
	* to give alp2 in [0, pi/2]). */
	if cbet2 != cbet1 || fabs(sbet2) != -sbet1 {
		var v float64
		if cbet1 < -sbet1 {
			v = (cbet2 - cbet1) * (cbet1 + cbet2)
		} else {
			v = (sbet1 - sbet2) * (sbet1 + sbet2)
		}
		calp2 = sqrt(sq(calp1*cbet1)+v) / cbet2
	} else {
		calp2 = fabs(calp1)
	}

	/* tan(bet2) = tan(sig2) * cos(alp2)
	* tan(omg2) = sin(alp0) * tan(sig2). */
	ssig2 = sbet2
	somg2 = salp0 * sbet2
	comg2 = calp2 * cbet2
	csig2 = comg2
	norm2(&ssig2, &csig2)
	/* norm2(&somg2, &comg2); -- don't need to normalize! */

	/* sig12 = sig2 - sig1, limit to [0, pi] */
	aa1 := maxx(0, csig1*ssig2-ssig1*csig2)
	aa2 := csig1*csig2 + ssig1*ssig2
	sig12 = atan2(aa1, aa2)

	/* omg12 = omg2 - omg1, limit to [0, pi] */
	somg12 = maxx((real)(0), comg1*somg2-somg1*comg2)
	comg12 = comg1*comg2 + somg1*somg2
	/* eta = omg12 - lam120 */
	eta = atan2(somg12*clam120-comg12*slam120,
		comg12*clam120+somg12*slam120)
	k2 = sq(calp0) * g.ep2
	eps = k2 / (2*(1+sqrt(1+k2)) + k2)

	c3f(g, eps, Ca)
	B312 = (sinCosSeries(iTrue, ssig2, csig2, Ca, nC3-1) -
		sinCosSeries(iTrue, ssig1, csig1, Ca, nC3-1))

	domg12 = -g.f * a3f(g, eps) * salp0 * (sig12 + B312)
	lam12 = eta + domg12

	if diffp {
		if calp2 == 0 {
			dlam12 = -2 * g.f1 * dn1 / sbet1
		} else {
			lengths(g, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
				cbet1, cbet2, nil, &dlam12, nil, nil, nil, Ca)
			dlam12 *= g.f1 / (calp2 * cbet2)
		}
		*pdlam12 = dlam12
	}

	*psalp2 = salp2
	*pcalp2 = calp2
	*psig12 = sig12
	*pssig1 = ssig1
	*pcsig1 = csig1
	*pssig2 = ssig2
	*pcsig2 = csig2
	*peps = eps
	*pdomg12 = domg12

	return lam12
}

func c3f(g *geodGeodesic, eps float64, c []float64) {
	/* Evaluate C3 coeffs
	 * Elements c[1] through c[nC3 - 1] are set */
	mult := 1.0
	var o, l int
	for l = 1; l < nC3; l++ { /* l is index of C3[l] */
		m := nC3 - l - 1 /* order of polynomial in eps */
		mult *= eps
		c[l] = mult * polyval(m, g.C3x[o:], eps)
		o += m + 1
	}
}

func c4f(g *geodGeodesic, eps float64, c []float64) {
	/* Evaluate C4 coeffs
	 * Elements c[0] through c[nC4 - 1] are set */
	mult := 1.0
	var o, l int
	for l = 0; l < nC4; l++ { /* l is index of C4[l] */
		m := nC4 - l - 1 /* order of polynomial in eps */
		c[l] = mult * polyval(m, g.C4x[o:], eps)
		o += m + 1
		mult *= eps
	}
}

func copysign(x, y float64) float64 {
	return math.Copysign(x, y)
}

func mod(x, y float64) float64 {
	return math.Mod(x, y)
}

func atan(x float64) float64 {
	return math.Atan(x)
}

func atanh(x float64) float64 {
	return math.Atanh(x)
}

func remainder(x, y float64) float64 {
	return math.Remainder(x, y)
}

func cbrt(x float64) float64 {
	return math.Cbrt(x)
}

func atan2(y, x float64) float64 {
	return math.Atan2(y, x)
}

func sqrt(a float64) float64 {
	return math.Sqrt(a)
}

func hypot(p, q float64) float64 {
	// intentionally not using math.Hypot()
	return sqrt(p*p + q*q)
}

func fabs(x float64) float64 {
	return math.Abs(x)
}

func maxx(x, y float64) float64 {
	if x > y {
		return x
	}
	return y
}

func minx(x, y float64) float64 {
	if x < y {
		return x
	}
	return y
}

func sin(x float64) float64 {
	return math.Sin(x)
}

func cos(x float64) float64 {
	return math.Cos(x)
}

func sincos(x float64) (float64, float64) {
	return math.Sincos(x)
}

func geodDirect(g *geodGeodesic,
	lat1 float64, lon1 float64, azi1 float64,
	s12 float64,
	plat2 *float64, plon2 *float64, pazi2 *float64) {
	geodGenDirect(g, lat1, lon1, azi1, geodNoFlags, s12, plat2, plon2, pazi2,
		nil, nil, nil, nil, nil)
}

type geodGeodesicLine struct {
	lat1  float64 /**< the starting latitude */
	lon1  float64 /**< the starting longitude */
	azi1  float64 /**< the starting azimuth */
	a     float64 /**< the equatorial radius */
	f     float64 /**< the flattening */
	salp1 float64 /**< sine of \e azi1 */
	calp1 float64 /**< cosine of \e azi1 */
	a13   float64 /**< arc length to reference point */
	s13   float64 /**< distance to reference point */
	/**< @cond SKIP */
	b, c2, f1, salp0, calp0, k2,
	ssig1, csig1, dn1, stau1, ctau1, somg1, comg1,
	A1m1, A2m1, A3c, B11, B21, B31, A4, B41 float64

	C1a  [6 + 1]float64
	C1pa [6 + 1]float64
	C2a  [6 + 1]float64
	C3a  [6]float64
	C4a  [6]float64
	/**< @endcond */
	caps uint /**< the capabilities */
}

func geodGenDirect(
	g *geodGeodesic,
	lat1 float64, lon1 float64, azi1 float64,
	flags uint, s12A12 float64,
	plat2 *float64, plon2 *float64, pazi2 *float64,
	ps12 *float64, pm12 *float64, pM12 *float64, pM21 *float64,
	pS12 *float64,
) float64 {
	var l geodGeodesicLine
	var outmask uint
	if plat2 != nil {
		outmask |= geodLatitude
	}
	if plon2 != nil {
		outmask |= geodLongitude
	}
	if pazi2 != nil {
		outmask |= geodAzimuth
	}
	if ps12 != nil {
		outmask |= geosDistance
	}
	if pm12 != nil {
		outmask |= geodReducedLength
	}
	if pM12 != nil || pM21 != nil {
		outmask |= geodGeodesicScale
	}
	if pS12 != nil {
		outmask |= geodArea
	}
	var moremask uint
	if flags&geodArcMode == 0 {
		moremask = geodDistanceIn
	}
	geodLineInit(&l, g, lat1, lon1, azi1,
		/* Automatically supply GEOD_DISTANCE_IN if necessary */
		outmask|moremask)

	return geodGenPosition(&l, flags, s12A12,
		plat2, plon2, pazi2, ps12, pm12, pM12, pM21, pS12)
}

func geodLineInit(
	l *geodGeodesicLine,
	g *geodGeodesic,
	lat1 float64, lon1 float64, azi1 float64, caps uint) {
	var salp1, calp1 float64
	azi1 = angNormalize(azi1)
	/* Guard against underflow in salp0 */

	sincosdx(angRound(azi1), &salp1, &calp1)
	geodLineInitInt(l, g, lat1, lon1, azi1, salp1, calp1, caps)
}

func geodLineInitInt(
	l *geodGeodesicLine,
	g *geodGeodesic,
	lat1 float64, lon1 float64,
	azi1 float64, salp1 float64, calp1 float64,
	caps uint) {
	var cbet1, sbet1, eps float64
	l.a = g.a
	l.f = g.f
	l.b = g.b
	l.c2 = g.c2
	l.f1 = g.f1
	/* If caps is 0 assume the standard direct calculation */
	if caps != 0 {
		l.caps = geodDistanceIn | geodLongitude
	}
	/* always allow latitude and azimuth and unrolling of longitude */
	l.caps |= geodLatitude | geodAzimuth | geodLongUnroll

	l.lat1 = latFix(lat1)
	l.lon1 = lon1
	l.azi1 = azi1
	l.salp1 = salp1
	l.calp1 = calp1

	sincosdx(angRound(l.lat1), &sbet1, &cbet1)
	sbet1 *= l.f1
	/* Ensure cbet1 = +epsilon at poles */
	norm2(&sbet1, &cbet1)
	cbet1 = maxx(tiny, cbet1)
	l.dn1 = sqrt(1 + g.ep2*sq(sbet1))

	/* Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0), */
	l.salp0 = l.salp1 * cbet1 /* alp0 in [0, pi/2 - |bet1|] */
	/* Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
	* is slightly better (consider the case salp1 = 0). */
	l.calp0 = hypot(l.calp1, l.salp1*sbet1)
	/* Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
	* sig = 0 is nearest northward crossing of equator.
	* With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
	* With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
	* With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
	* Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
	* With alp0 in (0, pi/2], quadrants for sig and omg coincide.
	* No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
	* With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi. */
	l.ssig1 = sbet1
	l.somg1 = l.salp0 * sbet1
	if sbet1 != 0 || l.calp1 != 0 {
		l.comg1 = cbet1 * l.calp1
	} else {
		l.comg1 = 1
	}
	l.csig1 = l.comg1

	norm2(&l.ssig1, &l.csig1) /* sig1 in (-pi, pi] */
	/* norm2(somg1, comg1); -- don't need to normalize! */

	l.k2 = sq(l.calp0) * g.ep2
	eps = l.k2 / (2*(1+sqrt(1+l.k2)) + l.k2)

	if l.caps&capC1 != 0 {
		var s, c float64
		l.A1m1 = a1m1f(eps)
		c1f(eps, l.C1a[:])
		l.B11 = sinCosSeries(iTrue, l.ssig1, l.csig1, l.C1a[:], nC1)
		s, c = sincos(l.B11)
		/* tau1 = sig1 + B11 */
		l.stau1 = l.ssig1*c + l.csig1*s
		l.ctau1 = l.csig1*c - l.ssig1*s
		/* Not necessary because C1pa reverts C1a
		*    B11 = -SinCosSeries(TRUE, stau1, ctau1, C1pa, nC1p); */
	}

	if l.caps&capC1p != 0 {
		c1pf(eps, l.C1pa[:])
	}

	if l.caps&capC2 != 0 {
		l.A2m1 = a2m1f(eps)
		c2f(eps, l.C2a[:])
		l.B21 = sinCosSeries(iTrue, l.ssig1, l.csig1, l.C2a[:], nC2)
	}

	if l.caps&capC3 != 0 {
		c3f(g, eps, l.C3a[:])
		l.A3c = -l.f * l.salp0 * a3f(g, eps)
		l.B31 = sinCosSeries(iTrue, l.ssig1, l.csig1, l.C3a[:], nC3-1)
	}

	if l.caps&capC4 != 0 {
		c4f(g, eps, l.C4a[:])
		/* Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0) */
		l.A4 = sq(l.a) * l.calp0 * l.salp0 * g.e2
		l.B41 = sinCosSeries(iFalse, l.ssig1, l.csig1, l.C4a[:], nC4)
	}
	l.s13 = math.NaN() //NaN
	l.a13 = l.s13

}

/* The coefficients C1p[l] in the Fourier expansion of B1p */
func c1pf(eps float64, c []float64) {
	coeff := [...]float64{
		/* C1p[1]/eps^1, polynomial in eps2 of order 2 */
		205, -432, 768, 1536,
		/* C1p[2]/eps^2, polynomial in eps2 of order 2 */
		4005, -4736, 3840, 12288,
		/* C1p[3]/eps^3, polynomial in eps2 of order 1 */
		-225, 116, 384,
		/* C1p[4]/eps^4, polynomial in eps2 of order 1 */
		-7173, 2695, 7680,
		/* C1p[5]/eps^5, polynomial in eps2 of order 0 */
		3467, 7680,
		/* C1p[6]/eps^6, polynomial in eps2 of order 0 */
		38081, 61440,
	}
	eps2 := sq(eps)
	d := eps
	var o, l int
	for l = 1; l <= nC1p; l++ { /* l is index of C1p[l] */
		m := (nC1p - l) / 2 /* order of polynomial in eps^2 */
		c[l] = d * polyval(m, coeff[o:], eps2) / coeff[o+m+1]
		o += m + 2
		d *= eps
	}
}

func geodGenPosition(
	l *geodGeodesicLine,
	flags uint, s12A12 float64,
	plat2 *float64, plon2 *float64, pazi2 *float64,
	ps12 *float64, pm12 *float64,
	pM12 *float64, pM21 *float64,
	pS12 *float64,
) float64 {

	var lat2, lon2, azi2, s12, m12, M12, M21, S12 float64
	/* Avoid warning about uninitialized B12. */
	var sig12, ssig12, csig12, B12, AB1 float64
	var omg12, lam12, lon12 float64
	var ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2, dn2 float64
	var outmask uint
	if plat2 != nil {
		outmask |= geodLatitude
	}
	if plon2 != nil {
		outmask |= geodLongitude
	}
	if pazi2 != nil {
		outmask |= geodAzimuth
	}
	if ps12 != nil {
		outmask |= geosDistance
	}
	if pm12 != nil {
		outmask |= geodReducedLength
	}
	if pM12 != nil || pM21 != nil {
		outmask |= geodGeodesicScale
	}
	if pS12 != nil {
		outmask |= geodArea
	}
	outmask &= l.caps & outAll

	if !(flags&geodArcMode != 0 || (l.caps&(geodDistanceIn&outAll)) != 0) {
		// /* Impossible distance calculation requested */
		return math.NaN()
	}

	if flags&geodArcMode != 0 {
		/* Interpret s12_a12 as spherical arc length */
		sig12 = s12A12 * degree
		sincosdx(s12A12, &ssig12, &csig12)
	} else {
		/* Interpret s12_a12 as distance */
		tau12 := s12A12 / (l.b * (1 + l.A1m1))
		s, c := sincos(tau12)

		/* tau2 = tau1 + tau12 */
		B12 = -sinCosSeries(iTrue,
			l.stau1*c+l.ctau1*s,
			l.ctau1*c-l.stau1*s,
			l.C1pa[:], nC1p)

		sig12 = tau12 - (B12 - l.B11)
		ssig12, csig12 = sincos(sig12)
		if fabs(l.f) > 0.01 {
			/* Reverted distance series is inaccurate for |f| > 1/100, so correct
			* sig12 with 1 Newton iteration.  The following table shows the
			* approximate maximum error for a = WGS_a() and various f relative to
			* GeodesicExact.
			*     erri = the error in the inverse solution (nm)
			*     errd = the error in the direct solution (series only) (nm)
			*     errda = the error in the direct solution (series + 1 Newton) (nm)
			*
			*       f     erri  errd errda
			*     -1/5    12e6 1.2e9  69e6
			*     -1/10  123e3  12e6 765e3
			*     -1/20   1110 108e3  7155
			*     -1/50  18.63 200.9 27.12
			*     -1/100 18.63 23.78 23.37
			*     -1/150 18.63 21.05 20.26
			*      1/150 22.35 24.73 25.83
			*      1/100 22.35 25.03 25.31
			*      1/50  29.80 231.9 30.44
			*      1/20   5376 146e3  10e3
			*      1/10  829e3  22e6 1.5e6
			*      1/5   157e6 3.8e9 280e6 */
			ssig2 := l.ssig1*csig12 + l.csig1*ssig12
			csig2 := l.csig1*csig12 - l.ssig1*ssig12
			B12 = sinCosSeries(iTrue, ssig2, csig2, l.C1a[:], nC1)
			serr := (1+l.A1m1)*(sig12+(B12-l.B11)) - s12A12/l.b
			sig12 = sig12 - serr/sqrt(1+l.k2*sq(ssig2))
			ssig12, csig12 = sincos(sig12)
			/* Update B12 below */
		}
	}

	/* sig2 = sig1 + sig12 */
	ssig2 = l.ssig1*csig12 + l.csig1*ssig12
	csig2 = l.csig1*csig12 - l.ssig1*ssig12
	dn2 = sqrt(1 + l.k2*sq(ssig2))
	if outmask&(geosDistance|geodReducedLength|geodGeodesicScale) != 0 {
		if flags&geodArcMode != 0 || fabs(l.f) > 0.01 {
			B12 = sinCosSeries(iTrue, ssig2, csig2, l.C1a[:], nC1)
		}
		AB1 = (1 + l.A1m1) * (B12 - l.B11)
	}
	/* sin(bet2) = cos(alp0) * sin(sig2) */
	sbet2 = l.calp0 * ssig2
	/* Alt: cbet2 = hypot(csig2, salp0 * ssig2); */
	cbet2 = hypot(l.salp0, l.calp0*csig2)
	if cbet2 == 0 {
		// /* I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case */
		csig2 = tiny
		cbet2 = csig2
	}
	/* tan(alp0) = cos(sig2)*tan(alp2) */
	salp2 = l.salp0
	calp2 = l.calp0 * csig2 /* No need to normalize */

	if (outmask & geosDistance) != 0 {
		if flags&geodArcMode != 0 {
			s12 = l.b * ((1+l.A1m1)*sig12 + AB1)
		} else {
			s12 = s12A12
		}
	}

	if outmask&geodLongitude != 0 {
		E := copysign(1, l.salp0) /* east or west going? */
		/* tan(omg2) = sin(alp0) * tan(sig2) */
		somg2 = l.salp0 * ssig2
		comg2 = csig2 /* No need to normalize */
		/* omg12 = omg2 - omg1 */

		if (flags & geodLongUnroll) != 0 {
			omg12 = E * (sig12 - (atan2(ssig2, csig2) - atan2(l.ssig1, l.csig1)) + (atan2(E*somg2, comg2) - atan2(E*l.somg1, l.comg1)))
		} else {
			omg12 = atan2(somg2*l.comg1-comg2*l.somg1, comg2*l.comg1+somg2*l.somg1)
		}
		lam12 = omg12 + l.A3c*
			(sig12+(sinCosSeries(iTrue, ssig2, csig2, l.C3a[:], nC3-1)-l.B31))
		lon12 = lam12 / degree

		if (flags & geodLongUnroll) != 0 {
			lon2 = l.lon1 + lon12
		} else {
			lon2 = angNormalize(angNormalize(l.lon1) + angNormalize(lon12))
		}
	}

	if outmask&geodLatitude != 0 {
		lat2 = atan2dx(sbet2, l.f1*cbet2)
	}

	if (outmask & geodAzimuth) != 0 {
		azi2 = atan2dx(salp2, calp2)
	}

	if outmask&(geodReducedLength|geodGeodesicScale) != 0 {
		B22 := sinCosSeries(iTrue, ssig2, csig2, l.C2a[:], nC2)
		AB2 := (1 + l.A2m1) * (B22 - l.B21)
		J12 := (l.A1m1-l.A2m1)*sig12 + (AB1 - AB2)
		if outmask&geodReducedLength != 0 {
			/* Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
			* accurate cancellation in the case of coincident points. */
			m12 = l.b * ((dn2*(l.csig1*ssig2) - l.dn1*(l.ssig1*csig2)) - l.csig1*csig2*J12)
		}
		if outmask&geodGeodesicScale != 0 {
			t := l.k2 * (ssig2 - l.ssig1) * (ssig2 + l.ssig1) /
				(l.dn1 + dn2)
			M12 = csig12 + (t*ssig2-csig2*J12)*l.ssig1/l.dn1
			M21 = csig12 - (t*l.ssig1-l.csig1*J12)*ssig2/dn2
		}
	}

	if outmask&geodArea != 0 {
		B42 := sinCosSeries(iFalse, ssig2, csig2, l.C4a[:], nC4)
		var salp12, calp12 float64
		if l.calp0 == 0 || l.salp0 == 0 {
			/* alp12 = alp2 - alp1, used in atan2 so no need to normalize */
			salp12 = salp2*l.calp1 - calp2*l.salp1
			calp12 = calp2*l.calp1 + salp2*l.salp1
		} else {
			/* tan(alp) = tan(alp0) * sec(sig)
			* tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
			* = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
			* If csig12 > 0, write
			*   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
			* else
			*   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
			* No need to normalize */
			var v float64
			if csig12 <= 0 {
				v = l.csig1*(1-csig12) + ssig12*l.ssig1
			} else {
				v = ssig12 * (l.csig1*ssig12/(1+csig12) + l.ssig1)
			}
			salp12 = l.calp0 * l.salp0 * v
			calp12 = sq(l.salp0) + sq(l.calp0)*l.csig1*csig2
		}
		S12 = l.c2*atan2(salp12, calp12) + l.A4*(B42-l.B41)
	}

	/* In the pattern
	*
	*   if ((outmask & GEOD_XX) && pYY)
	*     *pYY = YY;
	*
	* the second check "&& pYY" is redundant.  It's there to make the CLang
	* static analyzer happy.
	 */
	if (outmask&geodLatitude) != 0 && plat2 != nil {
		*plat2 = lat2
	}
	if (outmask&geodLongitude) != 0 && plon2 != nil {
		*plon2 = lon2
	}
	if (outmask&geodAzimuth) != 0 && pazi2 != nil {
		*pazi2 = azi2
	}
	if (outmask&geosDistance) != 0 && ps12 != nil {
		*ps12 = s12
	}
	if (outmask&geodReducedLength) != 0 && pm12 != nil {
		*pm12 = m12
	}
	if (outmask & geodGeodesicScale) != 0 {
		if pM12 != nil {
			*pM12 = M12
		}
		if pM21 != nil {
			*pM21 = M21
		}
	}
	if (outmask&geodArea) != 0 && pS12 != nil {
		*pS12 = S12
	}
	if (flags & geodArcMode) != 0 {
		return s12A12
	}
	return sig12 / degree
}

/**
 * The struct for accumulating information about a geodesic polygon.  This is
 * used for computing the perimeter and area of a polygon.  This must be
 * initialized by geod_polygon_init() before use.
 **********************************************************************/
type geodPolygon struct {
	lat float64 /**< the current latitude */
	lon float64 /**< the current longitude */
	/**< @cond SKIP */
	lat0      float64
	lon0      float64
	A         [2]float64
	P         [2]float64
	polyline  bool
	crossings int
	/**< @endcond */
	num uint /**< the number of points so far */
}

func geodPolygonInit(p *geodPolygon, polylinep bool) {
	p.polyline = polylinep
	geodPolygonClear(p)
}

func geodPolygonClear(p *geodPolygon) {
	p.lat0, p.lon0, p.lat, p.lon =
		math.NaN(), math.NaN(), math.NaN(), math.NaN()
	accini(&p.P)
	accini(&p.A)
	p.num, p.crossings = 0, 0
}

func accini(s *[2]float64) {
	/* Initialize an accumulator; this is an array with two elements. */
	s[0], s[1] = 0, 0
}

func accadd(s *[2]float64, y float64) {
	/* Add y to an accumulator. */
	var u float64
	z := sumx(y, s[1], &u)
	s[0] = sumx(z, s[0], &s[1])
	if s[0] == 0 {
		s[0] = u
	} else {
		s[1] = s[1] + u
	}
}

func transit(lon1, lon2 float64) int {
	var lon12 float64
	/* Return 1 or -1 if crossing prime meridian in east or west direction.
	 * Otherwise return zero. */
	/* Compute lon12 the same way as Geodesic::Inverse. */
	lon1 = angNormalize(lon1)
	lon2 = angNormalize(lon2)
	lon12 = angDiff(lon1, lon2, nil)
	if lon1 <= 0 && lon2 > 0 && lon12 > 0 {
		return 1
	} else if lon2 <= 0 && lon1 > 0 && lon12 < 0 {
		return -1
	}
	return 0
}

func geodPolygonAddPoint(g *geodGeodesic, p *geodPolygon, lat, lon float64) {
	lon = angNormalize(lon)
	if p.num == 0 {
		p.lat0, p.lat = lat, lat
		p.lon0, p.lon = lon, lon
	} else {
		var s12, S12 float64
		if p.polyline {
			geodGenInverse(g, p.lat, p.lon, lat, lon,
				&s12, nil, nil, nil, nil, nil, nil)
		} else {
			geodGenInverse(g, p.lat, p.lon, lat, lon,
				&s12, nil, nil, nil, nil, nil, &S12)
		}
		accadd(&p.P, s12)
		if !p.polyline {
			accadd(&p.A, S12)
			p.crossings += transit(p.lon, lon)
		}
		p.lat = lat
		p.lon = lon
	}
	p.num++
}

func acccopy(s [2]float64, t *[2]float64) {
	/* Copy an accumulator; t = s. */
	t[0] = s[0]
	t[1] = s[1]
}

func accsum(s [2]float64, y float64) float64 {
	/* Return accumulator + y (but don't add to accumulator). */
	var t [2]float64
	acccopy(s, &t)
	accadd(&t, y)
	return t[0]
}

func accrem(s *[2]float64, y float64) {
	/* Reduce to [-y/2, y/2]. */
	s[0] = remainder(s[0], y)
	accadd(s, (real)(0))
}

func accneg(s *[2]float64) {
	/* Negate an accumulator. */
	s[0] = -s[0]
	s[1] = -s[1]
}

func areareduceA(area *[2]float64, area0 float64,
	crossings int, reverse, sign bool,
) float64 {
	accrem(area, area0)
	if crossings != 0 {
		var x float64
		if area[0] < 0 {
			x = 1
		} else {
			x = -1
		}
		accadd(area, x*area0/2)
	}
	/* area is with the clockwise sense.  If !reverse convert to
	 * counter-clockwise convention. */
	if !reverse {
		accneg(area)
	}
	/* If sign put area in (-area0/2, area0/2], else put area in [0, area0) */
	if sign {
		if area[0] > area0/2 {
			accadd(area, -area0)
		} else if area[0] <= -area0/2 {
			accadd(area, +area0)
		}
	} else {
		if area[0] >= area0 {
			accadd(area, -area0)
		} else if area[0] < 0 {
			accadd(area, +area0)
		}
	}
	return 0 + area[0]
}
func geodPolygonCompute(
	g *geodGeodesic, p *geodPolygon, reverse, sign bool, pA, pP *float64,
) uint {
	if p.num < 2 {
		if pP != nil {
			*pP = 0
		}
		if !p.polyline && pA != nil {
			*pA = 0
		}
		return p.num
	}
	if p.polyline {
		if pP != nil {
			*pP = p.P[0]
		}
		return p.num
	}
	var s12, S12 float64
	var t [2]float64
	geodGenInverse(g, p.lat, p.lon, p.lat0, p.lon0,
		&s12, nil, nil, nil, nil, nil, &S12)
	if pP != nil {
		*pP = accsum(p.P, s12)
	}
	acccopy(p.A, &t)
	accadd(&t, S12)
	if pA != nil {
		*pA = areareduceA(&t, 4*pi*g.c2,
			p.crossings+transit(p.lon, p.lon0),
			reverse, sign)
	}
	return p.num
}

func transitdirect(lon1, lon2 float64) int {
	/* Compute exactly the parity of
	   int(ceil(lon2 / 360)) - int(ceil(lon1 / 360)) */
	lon1 = remainder(lon1, (real)(720))
	lon2 = remainder(lon2, (real)(720))
	var x, y int
	if lon2 <= 0 && lon2 > -360 {
		x = 1
	}
	if lon1 <= 0 && lon1 > -360 {
		y = 1
	}
	return x - y
}

func geodPolygonAddEdge(g *geodGeodesic, p *geodPolygon, azi, s float64) {
	if p.num != 0 {
		var lat, lon, S12 float64
		if p.polyline {
			geodGenDirect(g, p.lat, p.lon, azi, geodLongUnroll, s,
				&lat, &lon, nil, nil, nil, nil, nil, &S12)
		} else {
			geodGenDirect(g, p.lat, p.lon, azi, geodLongUnroll, s,
				&lat, &lon, nil, nil, nil, nil, nil, nil)
		}

		accadd(&p.P, s)
		if !p.polyline {
			accadd(&p.A, S12)
			p.crossings += transitdirect(p.lon, lon)
		}
		p.lat = lat
		p.lon = lon
		p.num++
	}
}

const radians = math.Pi / 180
const degrees = 180 / math.Pi

func sphericalDestination(radius float64, lat1, lon1, meters, bearingDegrees float64) (lat2, lon2 float64) {
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

func sphericalDistance(radius float64, lat1, lon1, lat2, lon2 float64) float64 {
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

func sphericalBearing(lat1, lon1, lat2, lon2 float64) float64 {
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

func sphericalInverse(
	radius float64,
	lat1 float64, lon1 float64,
	lat2 float64, lon2 float64,
	s12 *float64, azi1 *float64, azi2 *float64,
) {
	if s12 != nil {
		*s12 = sphericalDistance(radius, lat1, lon1, lat2, lon2)
	}
	if azi1 != nil {
		*azi1 = sphericalBearing(lat1, lon1, lat2, lon2)
	}
	if azi2 != nil {
		*azi2 = wrap180(sphericalBearing(lat2, lon2, lat1, lon1) + 180)
	}
}

func sphericalDirect(
	radius float64,
	lat1 float64, lon1 float64, azi1 float64,
	s12 float64,
	lat2 *float64, lon2 *float64, azi2 *float64,
) {
	la2, lo2 := sphericalDestination(radius, lat1, lon1, s12, azi1)
	if lat2 != nil {
		*lat2 = la2
	}
	if lon2 != nil {
		*lon2 = lo2
	}
	if azi2 != nil {
		*azi2 = wrap180(sphericalBearing(la2, lo2, lat1, lon1) + 180)
	}
}

// WGS84 conforming ellispoid
// https://en.wikipedia.org/wiki/World_Geodetic_System
var WGS84 = NewEllipsoid(6378137, float64(1.)/298.257223563)

// Globe is a pre-initialized spherical representing Earth as a
// terrestrial globe.
var Globe = NewSpherical(6378137)

// Ellipsoid is an object for performing geodesic operations.
type Ellipsoid struct {
	g          geodGeodesic
	radius     float64
	flattening float64
	spherical  bool
}

// NewEllipsoid initializes a new geodesic ellipsoid object.
//
// Param radius is the equatorial radius (meters).
// Param flattening is the flattening factor of the ellipsoid.
//
// The WGS84 pacakge-level variable is a pre-initialized ellipsoid
// representing Earth.
func NewEllipsoid(radius, flattening float64) *Ellipsoid {
	e := &Ellipsoid{radius: radius, flattening: flattening}
	geodInit(&e.g, radius, flattening)
	return e
}

// NewSpherical initializes a new geodesic ellipsoid object that uses
// simplified operations on a sphere.
//
// The Inverse and Direct operations will often be more computationally
// efficient than NewEllipsoid because it uses simplier great-circle
// calculations such as the Haversine formula.
//
// Param radius is the equatorial radius (meters).
//
// The Globe package-level variable is a pre-initialized spherical
// representing Earth as a terrestrial globe.
func NewSpherical(radius float64) *Ellipsoid {
	e := NewEllipsoid(radius, 0)
	e.spherical = true
	return e
}

// Radius of the Ellipsoid
func (e *Ellipsoid) Radius() float64 {
	return e.radius
}

// Flattening of the Ellipsoid
func (e *Ellipsoid) Flattening() float64 {
	return e.flattening
}

// Spherical returns true if the ellipsoid was initialized using NewSpherical.
func (e *Ellipsoid) Spherical() bool {
	return e.spherical
}

// Inverse solve the inverse geodesic problem.
//
// Param lat1 is latitude of point 1 (degrees).
// Param lon1 is longitude of point 1 (degrees).
// Param lat2 is latitude of point 2 (degrees).
// Param lon2 is longitude of point 2 (degrees).
// Out param ps12 is a pointer to the distance from point 1 to point 2 (meters).
// Out param pazi1 is a pointer to the azimuth at point 1 (degrees).
// Out param pazi2 is a pointer to the (forward) azimuth at point 2 (degrees).
//
// lat1 and lat2 should be in the range [-90,+90].
// The values of azi1 and azi2 returned are in the range [-180,+180].
// Any of the "return" arguments, ps12, etc., may be replaced with nil, if you
// do not need some quantities computed.
//
// The solution to the inverse problem is found using Newton's method.  If
// this fails to converge (this is very unlikely in geodetic applications
// but does occur for very eccentric ellipsoids), then the bisection method
// is used to refine the solution.
func (e *Ellipsoid) Inverse(
	lat1, lon1, lat2, lon2 float64,
	s12, azi1, azi2 *float64,
) {
	if e.spherical {
		sphericalInverse(e.radius, lat1, lon1, lat2, lon2, s12, azi1, azi2)
	} else {
		geodInverse(&e.g, lat1, lon1, lat2, lon2, s12, azi1, azi2)
	}
}

// Direct solves the direct geodesic problem.
//
// Param g is a pointer to the geod_geodesic object specifying the ellipsoid.
// Param lat1 is the latitude of point 1 (degrees).
// Param lon1 is the longitude of point 1 (degrees).
// Param azi1 is the azimuth at point 1 (degrees).
// Param s12 is the distance from point 1 to point 2 (meters). negative is ok.
// Out param plat2 is a pointer to the latitude of point 2 (degrees).
// Out param plon2 is a pointer to the longitude of point 2 (degrees).
// Out param pazi2 is a pointer to the (forward) azimuth at point 2 (degrees).
//
// lat1 should be in the range [-90,+90].
// The values of lon2 and azi2 returned are in the range [-180,+180].
// Any of the "return" arguments, plat2, etc., may be replaced with nil, if you
// do not need some quantities computed.
func (e *Ellipsoid) Direct(
	lat1, lon1, azi1, s12 float64,
	lat2, lon2, azi2 *float64,
) {
	if e.spherical {
		sphericalDirect(e.radius, lat1, lon1, azi1, s12, lat2, lon2, azi2)
	} else {
		geodDirect(&e.g, lat1, lon1, azi1, s12, lat2, lon2, azi2)
	}
}

// Polygon struct for accumulating information about a geodesic polygon.
// Used for computing the perimeter and area of a polygon.
// This must be initialized from Ellipsoid.PolygonInit before use.
type Polygon struct {
	e *Ellipsoid
	p geodPolygon
}

// PolygonInit initializes a polygon.
// Param polyline for polyline instead of a polygon.
//
// If polyline is not set, then the sequence of vertices and edges added by
// Polygon.AddPoint() and Polygon.AddEdge() define a polygon and
// the perimeter and area are returned by Polygon.Compute().
// If polyline is set, then the vertices and edges define a polyline and
// only the perimeter is returned by Polygon.Compute().
//
// The area and perimeter are accumulated at two times the standard floating
// point precision to guard against the loss of accuracy with many-sided
// polygons.  At any point you can ask for the perimeter and area so far.
func (e *Ellipsoid) PolygonInit(polyline bool) Polygon {
	var p Polygon
	p.e = e
	geodPolygonInit(&p.p, polyline)
	return p
}

// AddPoint adds a point to the polygon or polyline.
//
// Param lat is the latitude of the point (degrees).
// Param lon is the longitude of the point (degrees).
func (p *Polygon) AddPoint(lat, lon float64) {
	geodPolygonAddPoint(&p.e.g, &p.p, lat, lon)
}

// Compute the results for a polygon
//
// Param reverse, if set then clockwise (instead of
//   counter-clockwise) traversal counts as a positive area.
// Param sign, if set then return a signed result for the area if
//   the polygon is traversed in the "wrong" direction instead of returning
//   the area for the rest of the earth.
// Out param pA is a pointer to the area of the polygon (meters-squared);
// Out param pP is a pointer to the perimeter of the polygon or length of the
//   polyline (meters).
// Returns the number of points.
//
// The area and perimeter are accumulated at two times the standard floating
// point precision to guard against the loss of accuracy with many-sided
// polygons.  Arbitrarily complex polygons are allowed.  In the case of
// self-intersecting polygons the area is accumulated "algebraically", e.g.,
// the areas of the 2 loops in a figure-8 polygon will partially cancel.
// There's no need to "close" the polygon by repeating the first vertex.  Set
// pA or pP to nil, if you do not want the corresponding quantity returned.
//
// More points can be added to the polygon after this call.
func (p *Polygon) Compute(clockwise, sign bool, area, perimeter *float64) int {
	return int(geodPolygonCompute(&p.e.g, &p.p, clockwise, sign, area, perimeter))
}

// AddEdge adds an edge to the polygon or polyline.
//
// Param azi is the azimuth at current point (degrees).
// Param s is the distance from current point to next point (meters).
func (p *Polygon) AddEdge(azi, s float64) {
	geodPolygonAddEdge(&p.e.g, &p.p, azi, s)
}

// Clear the polygon, allowing a new polygon to be started.
func (p *Polygon) Clear() {
	geodPolygonClear(&p.p)
}

func sign(a float64) float64 {
	if a < 0 {
		return -1
	}
	if a > 0 {
		return 1
	}
	return 0
}

// crossAlongSuffix performs both cross-track and/or along-track operations
func crossAlongSuffix(
	e *Ellipsoid,
	// position of point
	lat, lon float64,
	// start and azimuth of segment
	startLat, startLon float64, azi1 float64,
	// return values
	cross, along *float64,
) {
	r := e.Radius()
	θ12 := azi1 * radians
	var a13 float64
	var d13 float64
	e.Inverse(startLat, startLon, lat, lon, &d13, &a13, nil)
	δ13 := d13 / r
	θ13 := a13 * radians
	sδ13, cδ13 := math.Sincos(δ13)
	δxt := math.Asin(sδ13 * math.Sin(θ13-θ12))
	if cross != nil {
		dXt := δxt * r
		*cross = dXt
	}
	if along != nil {
		δat := math.Acos(cδ13 / math.Abs(math.Cos(δxt)))
		*along = δat * sign(math.Cos(θ12-θ13)) * r
	}
}

// CrossAlong provides the cross-track and along-track distances from
// point (lat/lon) to segment (start/end).
//
// Out param cross is the cross-track distance from point to segment
// Out param along is the along-track distance from point to segment
func (e *Ellipsoid) CrossAlong(
	lat, lon float64,
	startLat, startLon, endLat, endLon float64,
	cross, along *float64,
) {
	var azi1 float64
	e.Inverse(startLat, startLon, endLat, endLon, nil, &azi1, nil)
	if cross != nil || along != nil {
		crossAlongSuffix(e, lat, lon, startLat, startLon, azi1, cross, along)
	}
}
