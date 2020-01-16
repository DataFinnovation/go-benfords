// Package benfords generates Benford's distributions.  It is done
// roughly in the style of a Gonum single-variable distribution.
// The key difference is that the domain of this version is
// a set of ints.
package benfords

import (
	"log"
	"math"
	"strconv"

	"golang.org/x/exp/rand"
	"gonum.org/v1/gonum/stat"
)

// Benfords represents a Benford's distribution https://en.wikipedia.org/wiki/Benford%27s_law
type Benfords struct {
	// Base is the base in which numbers are considered
	// Standard decimal has Base = 10
	Base int

	Src rand.Source
}

// Prob returns the probability distribution function for a given digit
func (b Benfords) Prob(x int) float64 {
	if x < 1 {
		return 0
	}
	if x > (b.Base - 1.0) {
		return 0
	}
	return math.Log(1.0+1.0/float64(x)) / math.Log(float64(b.Base))
}

// LogProb returns the natural log of the probability at x
func (b Benfords) LogProb(x int) float64 {
	if x < 1 || x > (b.Base-1) {
		return math.Inf(-1)
	}
	p := b.Prob(x)
	if p == 0.0 {
		return math.Inf(-1)
	}
	return p
}

// CDF returns the cumulative distribution up to x
func (b Benfords) CDF(x int) float64 {
	if x < 1 {
		return 0
	}
	if x > (b.Base - 1) {
		return 1.0
	}
	tot := 0.0
	for i := 1; i < b.Base; i++ {
		tot += b.Prob(i)
	}
	return tot
}

// NumParameters returns the number of parameters in the distribution.
func (Benfords) NumParameters() int {
	return 1
}

// FullPDF returns the full discrete PDF of the distribution
func (b Benfords) FullPDF() []float64 {
	res := make([]float64, b.Base-1)
	for i := range res {
		res[i] = b.Prob(i + 1)
	}
	return res
}

// FullCDF returns the full discrete CDF of the distribution
func (b Benfords) FullCDF() []float64 {
	res := make([]float64, b.Base-1)
	runningTotal := 0.0
	for i := range res {
		runningTotal += b.Prob(i + 1)
		res[i] = runningTotal
	}
	return res
}

// Domain returns all possible non-zero-probability values
func (b Benfords) Domain() []int {
	res := make([]int, b.Base-1)
	for i := range res {
		res[i] = i + 1
	}
	return res
}

// Rand returns a random sample from the distribution.
func (b Benfords) Rand() int {
	runif := rand.Float64
	if b.Src != nil {
		rnd := rand.New(b.Src)
		runif = rnd.Float64
	}

	p := runif()
	cdf := b.FullCDF()
	domain := b.Domain()
	for i, v := range cdf {
		if p < v {
			return domain[i]
		}
	}
	return (b.Base - 1)
}

// ChoGainesStat returns the Cho-Gaines statistic
func (b Benfords) ChoGainesStat(nSamples int, realisedDist []float64) float64 {
	idealPDF := b.FullPDF()
	if len(idealPDF) != len(realisedDist) {
		log.Panic("length mismatch")
	}
	totDiff := 0.0
	for i, r := range realisedDist {
		totDiff += math.Pow(r-idealPDF[i], 2.0)
	}
	return math.Sqrt(float64(nSamples) * totDiff)
}

// LeemisStat returns the Leemis statistic
func (b Benfords) LeemisStat(nSamples int, realisedDist []float64) float64 {
	idealPDF := b.FullPDF()
	if len(idealPDF) != len(realisedDist) {
		log.Panic("length mismatch")
	}
	maxDiff := 0.0
	for i, r := range realisedDist {
		maxDiff = math.Max(maxDiff, math.Abs(r-idealPDF[i]))
	}
	return math.Sqrt(float64(nSamples)) * maxDiff
}

// LeadDigit returns the leading digit of n in base base
func LeadDigit(n float64, base int) int {
	if n < 0 {
		n = -n
	}
	if n < 1 {
		for n < 1 {
			n = float64(base) * n
		}
	}
	resid := int(n)
	for resid >= base {
		resid = resid / base
	}
	return resid
}

// ComputeLeadDigitDistribution takes a vector of values and computes
// the first-digit distribution in the given base
// this function pulls out 0s and NaNs
// it returns both the distribution and the number of useful data points
// that were found
func ComputeLeadDigitDistribution(values []float64, base int) ([]float64, int) {
	validValues := 0.0
	dist := make([]float64, base-1)
	for _, v := range values {
		if v != 0.0 && !math.IsNaN(v) {
			validValues++
			leadDigit := LeadDigit(v, base)
			dist[leadDigit-1]++
		}
	}
	for i, v := range dist {
		dist[i] = v / validValues
	}
	return dist, int(validValues)
}

// ComputeLeadDigitDistributionFromStrings does what it says
// relying on the strconv.ParseFloat function
func ComputeLeadDigitDistributionFromStrings(strings []string, base int) ([]float64, int) {
	numericValues := make([]float64, len(strings))
	count := 0
	for _, v := range strings {
		asFloat, err := strconv.ParseFloat(v, 64)
		if err == nil && asFloat != 0.0 && !math.IsNaN(asFloat) {
			numericValues[count] = asFloat
			count++
		}
	}
	return ComputeLeadDigitDistribution(numericValues[0:count], base)
}

// ChiSquarePValue for the given distribution
func (b Benfords) ChiSquarePValue(dist []float64) float64 {
	if len(dist) != (b.Base - 1) {
		log.Panic("distribution has unexpected length")
	}
	idealPDF := b.FullPDF()
	return stat.ChiSquare(idealPDF, dist)
}

// eof
