package goitm

import (
	"math"
)

const pi = 3.141592653589793

func sin2(x float64) float64 { return math.Sin(x) * math.Sin(x) }
func cos2(x float64) float64 { return math.Cos(x) * math.Cos(x) }
func tan2(x float64) float64 { return math.Tan(x) * math.Tan(x) }
func tan4(x float64) float64 { return tan2(x) * tan2(x) }

type gGrid int

const (
	gICS gGrid = iota
	gITM
)

// GRID ok
type GRID struct {
	lon0   float64
	lat0   float64
	k0     float64
	falseE float64
	falseN float64
}

type eDatum int

const (
	eWGS84 eDatum = iota
	eGRS80
	eCLARK80M
)

var grid = []GRID{

	// ICS data
	GRID{
		0.6145667421719,     // lon0 = central meridian in radians of 35.12'43.490"
		0.55386447682762762, // lat0 = central latitude in radians of 31.44'02.749"
		1.00000,             // k0 = scale factor
		170251.555,          // false_easting
		2385259.0,           // false_northing
	},

	// ITM data
	GRID{
		0.61443473225468920, // lon0 = central meridian in radians 35.12'16.261"
		0.55386965463774187, // lat0 = central latitude in radians 31.44'03.817"
		1.0000067,           // k0 = scale factor
		219529.584,          // false_easting
		2885516.9488,        // false_northing = 3512424.3388-626907.390
		// MAPI says the false northing is 626907.390, and in another place
		// that the meridional arc at the central latitude is 3512424.3388
	},
}

type datum struct {
	a   float64 // a  Equatorial earth radius
	b   float64 // b  Polar earth radius
	f   float64 // f= (a-b)/a  Flatenning
	esq float64 // esq = 1-(b*b)/(a*a)  Eccentricity Squared
	e   float64 // sqrt(esq)  Eccentricity
	// deltas to WGS84
	dX float64
	dY float64
	dZ float64
}

// Datum ok
var Datum = [3]datum{

	// WGS84 data
	datum{
		6378137.0,            // a
		6356752.3142,         // b
		0.00335281066474748,  // f = 1/298.257223563
		0.006694380004260807, // esq
		0.0818191909289062,   // e
		// deltas to WGS84
		0,
		0,
		0},

	// GRS80 data
	datum{
		6378137.0,           // a
		6356752.3141,        // b
		0.0033528106811823,  // f = 1/298.257222101
		0.00669438002290272, // esq
		0.0818191910428276,  // e
		// deltas to WGS84
		-48,
		55,
		52},

	// Clark 1880 Modified data
	datum{
		6378300.789,          // a
		6356566.4116309,      // b
		0.003407549767264,    // f = 1/293.466
		0.006803488139112318, // esq
		0.08248325975076590,  // e
		// deltas to WGS84
		-235,
		-85,
		264}}

// Wgs842itm does WGS84 to Israel New Grid (ITM) conversion
func Wgs842itm(lat, lon float64) (int, int) {
	latr := lat * pi / 180
	lonr := lon * pi / 180

	// 1. Molodensky WGS84 -> GRS80

	lat80, lon80 := molodensky(latr, lonr, eWGS84, eGRS80)

	// 2. Lat/Lon (GRS80) -> Local Grid (ITM)
	return latLon2Grid(lat80, lon80, eGRS80, gITM)
}

// Wgs842ics does WGS84 to Israel Old Grid (ICS) conversion
func Wgs842ics(lat, lon float64) (int, int) {
	latr := lat * pi / 180
	lonr := lon * pi / 180

	// 1. molodensky WGS84 -> Clark_1880_modified
	lat80, lon80 := molodensky(latr, lonr, eWGS84, eCLARK80M)

	// 2. Lat/Lon (Clark_1880_modified) -> Local Grid (ICS)
	N, E := latLon2Grid(lat80, lon80, eCLARK80M, gICS)
	return E, N
}

// Itm2wgs84 converts Isreal ITM to WGS*$
func Itm2wgs84(N, E float64) (float64, float64) {
	// 1. Local grid (ITM) -> GRS80
	// var lat80, lon80 float64
	lat80, lon80 := Grid2LatLon(N, E, gITM, eGRS80)

	// 2. molodensky GRS80->WGS84

	lat84, lon84 := molodensky(lat80, lon80, eGRS80, eWGS84)

	// final results
	return lat84 * 180.0 / pi, lon84 * 180.0 / pi

}

// Ics2wgs84 converts Israel Old Grid (ICS) to WGS84 conversion
func Ics2wgs84(N, E float64) (float64, float64) {
	// 1. Local Grid (ICS) -> Clark_1880_modified
	lat80, lon80 := Grid2LatLon(N, E, gICS, eCLARK80M)

	// 2. molodensky Clark_1880_modified -> WGS84
	lat84, lon84 := molodensky(lat80, lon80, eCLARK80M, eWGS84)

	// final results
	lat := lat84 * 180 / pi
	lon := lon84 * 180 / pi
	return lat, lon
}

// Grid2LatLon Local grid to Lat/Lon conversion
//====================================
func Grid2LatLon(N, E float64, from gGrid, to eDatum) (float64, float64) {
	//================
	// GRID -> Lat/Lon
	//================
	y := N + grid[from].falseN
	x := E - grid[from].falseE
	M := y / grid[from].k0
	a := Datum[to].a
	b := Datum[to].b
	e := Datum[to].e
	esq := Datum[to].esq
	mu := M / (a * (1 - e*e/4 - 3*math.Pow(e, 4)/64 - 5*math.Pow(e, 6)/256))

	ee := math.Sqrt(1 - esq)
	e1 := (1 - ee) / (1 + ee)
	j1 := 3*e1/2 - 27*e1*e1*e1/32
	j2 := 21*e1*e1/16 - 55*e1*e1*e1*e1/32
	j3 := 151 * e1 * e1 * e1 / 96
	j4 := 1097 * e1 * e1 * e1 * e1 / 512
	// Footprint Latitude
	fp := mu + j1*math.Sin(2*mu) + j2*math.Sin(4*mu) + j3*math.Sin(6*mu) + j4*math.Sin(8*mu)

	sinfp := math.Sin(fp)
	cosfp := math.Cos(fp)
	tanfp := sinfp / cosfp
	eg := (e * a / b)
	eg2 := eg * eg
	C1 := eg2 * cosfp * cosfp
	T1 := tanfp * tanfp
	R1 := a * (1 - e*e) / math.Pow(1-(e*sinfp)*(e*sinfp), 1.5)
	N1 := a / math.Sqrt(1-(e*sinfp)*(e*sinfp))
	D := x / (N1 * grid[from].k0)

	Q1 := N1 * tanfp / R1
	Q2 := D * D / 2
	Q3 := (5 + 3*T1 + 10*C1 - 4*C1*C1 - 9*eg2*eg2) * (D * D * D * D) / 24
	Q4 := (61 + 90*T1 + 298*C1 + 45*T1*T1 - 3*C1*C1 - 252*eg2*eg2) * (D * D * D * D * D * D) / 720
	// result lat
	lat := fp - Q1*(Q2-Q3+Q4)

	Q5 := D
	Q6 := (1 + 2*T1 + C1) * (D * D * D) / 6
	Q7 := (5 - 2*C1 + 28*T1 - 3*C1*C1 + 8*eg2*eg2 + 24*T1*T1) * (D * D * D * D * D) / 120
	// result lon
	lon := grid[from].lon0 + (Q5-Q6+Q7)/cosfp
	return lat, lon
}

// molodensky  Abridged transformation between 2 datums
//======================================================
func molodensky(ilat, ilon float64, from, to eDatum) (float64, float64) {
	// from->WGS84 - to->WGS84 = from->WGS84 + WGS84->to = from->to
	dX := Datum[from].dX - Datum[to].dX
	dY := Datum[from].dY - Datum[to].dY
	dZ := Datum[from].dZ - Datum[to].dZ
	slat := math.Sin(ilat)
	clat := math.Cos(ilat)
	slon := math.Sin(ilon)
	clon := math.Cos(ilon)
	ssqlat := slat * slat

	//dlat = ((-dx * slat * clon - dy * slat * slon + dz * clat)
	//        + (da * rn * fromEsq * slat * clat / fromA)
	//        + (df * (rm * adb + rn / adb )* slat * clat))
	//       / (rm + from.h);

	fromF := Datum[from].f
	df := Datum[to].f - fromF
	fromA := Datum[from].a
	da := Datum[to].a - fromA
	fromEsq := Datum[from].esq
	adb := 1.0 / (1.0 - fromF)
	rn := fromA / math.Sqrt(1-fromEsq*ssqlat)
	rm := fromA * (1 - fromEsq) / math.Pow((1-fromEsq*ssqlat), 1.5)
	fromH := 0.0 // we're flat!
	dlat := (-dX*slat*clon - dY*slat*slon + dZ*clat + da*rn*fromEsq*slat*clat/fromA +
		+df*(rm*adb+rn/adb)*slat*clat) /
		(rm + fromH)

	// result lat (radians)
	olat := ilat + dlat

	// dlon = (-dx * slon + dy * clon) / ((rn + from.h) * clat);
	dlon := (-dX*slon + dY*clon) / ((rn + fromH) * clat)
	// result lon (radians)
	olon := ilon + dlon
	return olat, olon
}

//====================================
// Lat/Lon to Local Grid conversion
//====================================
func latLon2Grid(lat, lon float64, from eDatum, to gGrid) (int, int) {
	// Datum data for Lat/Lon to TM conversion
	a := Datum[from].a
	e := Datum[from].e // sqrt(esq);
	b := Datum[from].b

	//===============
	// Lat/Lon -> TM
	//===============
	slat1 := math.Sin(lat)
	clat1 := math.Cos(lat)
	clat1sq := clat1 * clat1
	tanlat1sq := slat1 * slat1 / clat1sq
	e2 := e * e
	e4 := e2 * e2
	e6 := e4 * e2
	eg := (e * a / b)
	eg2 := eg
	l1 := 1 - e2/4 - 3*e4/64 - 5*e6/256
	l2 := 3*e2/8 + 3*e4/32 + 45*e6/1024
	l3 := 15*e4/256 + 45*e6/1024
	l4 := 35 * e6 / 3072
	M := a * (l1*lat - l2*math.Sin(2*lat) + l3*math.Sin(4*lat) - l4*math.Sin(6*lat))
	//double rho = a*(1-e2) / pow((1-(e*slat1)*(e*slat1)),1.5);
	nu := a / math.Sqrt(1-(e*slat1)*(e*slat1))
	p := lon - grid[to].lon0
	k0 := grid[to].k0
	// y = northing = K1 + K2p2 + K3p4, where
	K1 := M * k0
	K2 := k0 * nu * slat1 * clat1 / 2
	K3 := (k0 * nu * slat1 * clat1 * clat1sq / 24) * (5 - tanlat1sq + 9*eg2*clat1sq + 4*eg2*eg2*clat1sq*clat1sq)
	// ING north
	Y := K1 + K2*p*p + K3*p*p*p*p - grid[to].falseN

	// x = easting = K4p + K5p3, where
	K4 := k0 * nu * clat1
	K5 := (k0 * nu * clat1 * clat1sq / 6) * (1 - tanlat1sq + eg2*clat1*clat1)
	// ING east
	X := K4*p + K5*p*p*p + grid[to].falseE

	// final rounded results
	E := int(X + 0.5)
	N := int(Y + 0.5)
	return E, N
}
