package benfords_test

import (
	"math"
	"math/rand"
	"testing"

	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

func TestBenfords(t *testing.T) {
	const samples = 1000000
	bases := []int{4, 8, 10, 16, 32}

	for _, base := range bases {
		thisDist := distuv.Benfords{Base: base}
		pdf := thisDist.FullPDF()

		realised := make([]float64, base-1)
		for i := 0; i < samples; i++ {
			thisSample := thisDist.Rand()
			realised[thisSample-1] += 1.0 / float64(samples)
		}
		//		fmt.Println("ideal:", pdf)
		//		fmt.Println("realised:", realised)

		pstat := stat.ChiSquare(realised, pdf)
		//		fmt.Println("pstat:", pstat)
		if pstat > 0.01 {
			t.Error("Chi Squared stat looks invalid base ", base)
		}
		if base == 10 {
			// critical values only know for base 10
			cg := thisDist.ChoGainesStat(samples, realised)
			//		fmt.Println("cg:", cg)
			if cg > 2 {
				t.Error("Cho-Gaines stat looks invalid: ", cg)

			}
			leemis := thisDist.LeemisStat(samples, realised)
			//		fmt.Println("leemis:", leemis)
			if leemis > 1.5 {
				t.Error("Leemis stat looks invalid: ", leemis)

			}
		}
	}
}
func TestLeadDigit(t *testing.T) {
	const samples = 1000
	const maxDigits = 12
	for _, test := range []struct {
		base        int
		numbers     []float64
		firstDigits []int
	}{
		{
			base:        10,
			numbers:     []float64{999, 888, 777, 666, 987, 876, 765, 0.002, -321},
			firstDigits: []int{9, 8, 7, 6, 9, 8, 7, 2, 3},
		},
		{
			base:        9,
			numbers:     []float64{8, 80},
			firstDigits: []int{8, 8},
		},
		{
			base:        16,
			numbers:     []float64{1 * 16, 3 * 16, 5*16 + 4*16*16},
			firstDigits: []int{1, 3, 4},
		},
		{base: 5, numbers: nil, firstDigits: nil},
		{base: 7, numbers: nil, firstDigits: nil},
		{base: 9, numbers: nil, firstDigits: nil},
		{base: 13, numbers: nil, firstDigits: nil},
		{base: 19, numbers: nil, firstDigits: nil},
	} {
		base := test.base
		numbers := test.numbers
		firstDigits := test.firstDigits

		for i, n := range numbers {
			expectedFirst := firstDigits[i]
			foundFirst := distuv.LeadDigit(n, base)
			if expectedFirst != foundFirst {
				t.Error("lead digit issue:", expectedFirst, " , ", foundFirst)
			}
		}

		for i := 0; i < samples; i++ {
			digits := rand.Intn(maxDigits) + 1
			number := 0.0
			lastDigit := 0
			for d := 0; d < digits; d++ {
				lastDigit = 1 + rand.Intn(test.base-1)
				number += math.Pow(float64(test.base), float64(d)) * float64(lastDigit)
			}
			foundFirst := distuv.LeadDigit(number, test.base)
			if lastDigit != foundFirst {
				t.Error("lead digit issue:", lastDigit, " , ", foundFirst)
			}
		}

	}
}

// eof
