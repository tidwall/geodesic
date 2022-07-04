package geodesic

import (
	_ "embed"
	"encoding/binary"
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"

	"github.com/stretchr/testify/require"
)

// The "test.data" file has been generated from the
// https://github.com/tidwall/geodesic_cgo project.

//go:embed test.data
var testData []byte

func eqish(x, y float64, prec int) bool {
	return math.Abs(x-y) < float64(1.0)/math.Pow10(prec)
}

func readFloats(src []byte, count int) []float64 {
	vals := make([]float64, count)
	for i := 0; i < count; i++ {
		vals[i] = math.Float64frombits(binary.LittleEndian.Uint64(src[i*8:]))
	}
	return vals
}

func TestInput(t *testing.T) {
	data := testData
	for i := 0; i < len(data); {
		if data[i] == 'I' {
			v := readFloats(data[i+1:], 7)
			i += 1 + 7*8
			testInverse(t, v[0], v[1], v[2], v[3], v[4], v[5], v[6])
			testDirect(t, v[0], v[1], v[2], v[3], v[4], v[5], v[6])
		} else if data[i] == 'P' {
			i++
			npoints := int(data[i])
			i++
			points := readFloats(data[i:], npoints*2)
			i += npoints * 2 * 8
			vals := readFloats(data[i:], 8)
			i += 8 * 8
			testPolygon(t, points, vals)
		} else {
			t.Fatalf("invalid test data")
		}
	}
}

func testPolygon(t *testing.T, points []float64, vals []float64) {
	t.Helper()

	p := WGS84.PolygonInit(false)
	for i := 0; i < len(points); i += 2 {
		p.AddPoint(points[i+0], points[i+1])
	}
	retvals := make([]float64, 8)
	p.Compute(false, false, &retvals[0], &retvals[1])
	p.Compute(true, false, &retvals[2], &retvals[3])
	p.Compute(true, true, &retvals[4], &retvals[5])
	p.Compute(false, true, &retvals[6], &retvals[7])
	for i := 0; i < len(vals); i++ {
		if !eqish(vals[i], retvals[i], 3) {
			t.Fatalf("expected %f, got %f", vals, retvals)
		}
	}
}

func testInverse(t *testing.T, lat1, lon1, lat2, lon2, s12, azi1, azi2 float64) {
	t.Helper()
	var s12ret, azi1ret, azi2ret float64
	WGS84.Inverse(lat1, lon1, lat2, lon2, &s12ret, &azi1ret, &azi2ret)
	if !eqish(s12ret, s12, 7) || !eqish(azi1ret, azi1, 7) || !eqish(azi2ret, azi2, 7) {
		t.Fatalf("expected '%f, %f, %f', got '%f, %f, %f'",
			s12, azi1, azi2, s12ret, azi1ret, azi2ret)
	}
}

func testDirect(t *testing.T, lat1, lon1, lat2, lon2, s12, azi1, azi2 float64) {
	t.Helper()
	var lat2ret, lon2ret, azi2ret float64
	WGS84.Direct(lat1, lon1, azi1, s12, &lat2ret, &lon2ret, &azi2ret)
	if !eqish(lat2ret, lat2, 7) || !eqish(lon2ret, lon2, 7) || !eqish(azi2ret, azi2, 7) {
		t.Logf("direct   'lat1: %f, lon1: %f, azi1: %f, s12: %f'\n", lat1, lon1, azi1, s12)
		t.Logf("expected 'lat2: %f, lon2: %f, azi2: %f'\n", lat2, lon2, azi2)
		t.Logf("got      'lat2: %f, lon2: %f, azi2: %f'", lat2ret, lon2ret, azi2ret)
		t.FailNow()
	}
}

func TestSpherical(t *testing.T) {
	if !Globe.Spherical() {
		t.Fatal()
	}
	if Globe.Flattening() != 0 {
		t.Fatal()
	}
	if wrap180(-181) != 179 {
		t.Fatal()
	}
	if wrap180(+181) != -179 {
		t.Fatal()
	}

	rng := rand.New(rand.NewSource(time.Now().UnixNano()))

	e := NewEllipsoid(Globe.Radius(), 0)
	for i := 0; i < 1_000_000; i++ {
		lat1 := rng.Float64()*180 - 90
		lon1 := rng.Float64()*360 - 180
		lat2 := rng.Float64()*180 - 90
		lon2 := rng.Float64()*360 - 180

		var s12, azi1, azi2 float64
		e.Inverse(lat1, lon1, lat2, lon2, &s12, &azi1, &azi2)

		var ret [3]float64
		Globe.Inverse(lat1, lon1, lat2, lon2, &ret[0], &ret[1], &ret[2])
		if !eqish(ret[0], s12, 4) ||
			!eqish(ret[1], azi1, 4) ||
			!eqish(ret[2], azi2, 4) {
			t.Fatalf("inverse failure (%f %f %f %f %f %f %f)",
				lat1, lon1, lat2, lon2, s12, azi1, azi2)
		}
		Globe.Direct(lat1, lon1, azi1, s12, &ret[0], &ret[1], &ret[2])
		if !eqish(ret[0], lat2, 4) ||
			!eqish(ret[1], lon2, 4) ||
			!eqish(ret[2], azi2, 4) {
			t.Fatalf("direct failure (%f %f %f %f %f %f %f)",
				lat1, lon1, lat2, lon2, s12, azi1, azi2)
		}
	}
}

func TestLineInitFlags(t *testing.T) {
	lat0, lon0, azi0 := 37.750000, -17.500000, 5.793761558385
	line := WGS84.LineInit(lat0, lon0, azi0,
		Distance|DistanceIn|ReducedLength|GeodesicScale|Area,
	)

	s := 1899641.643510
	var lat2, lon2, azi2, s12, m12, M12, M21, S12 float64
	line.GenPosition(NoFlags,
		s, &lat2, &lon2, &azi2, &s12, &m12, &M12, &M21, &S12,
	)

	if lat2 == 0 {
		t.Error("gen pos: missing lat2")
	}
	if lon0 == 0 {
		t.Error("gen pos: missing lon0")
	}
	if azi2 == 0 {
		t.Error("gen pos: missing azi2")
	}
	if s12 == 0 {
		t.Error("gen pos: missing s12")
	}
	if m12 == 0 {
		t.Error("gen pos: missing m12")
	}
	if M12 == 0 {
		t.Error("gen pos: missing M12")
	}
	if M21 == 0 {
		t.Error("gen pos: missing M21")
	}
	if S12 == 0 {
		t.Error("gen pos: missing S12")
	}
}

func TestRemquo(t *testing.T) {
	tests := []struct {
		numer, rem float64
		quo        int
	}{
		{numer: 10, rem: 10, quo: 0},
		{numer: 50, rem: 5, quo: 1},
		{numer: 100, rem: 10, quo: 2},
		{numer: 150, rem: 15, quo: 3},
	}

	for _, tc := range tests {
		t.Run(fmt.Sprintf("%g", tc.numer), func(t *testing.T) {
			r, q := remquo(tc.numer, 45)
			require.Equal(t, tc.rem, r)
			require.Equal(t, tc.quo, q)
		})
	}
}
