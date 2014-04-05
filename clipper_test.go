package clipper

import (
	"math"
	"math/rand"
	"testing"
	"time"
)

func init() {
	rand.Seed(time.Now().UTC().UnixNano())
}

func RandomPoly(maxWidth, maxHeight, vertCnt int) []*Point {
	result := make([]*Point, vertCnt)
	for i := 0; i < vertCnt; i++ {
		result[i] = &Point{rand.Intn(maxWidth), rand.Intn(maxHeight)}
	}
	return result
}

func TestRandom(t *testing.T) {
	scale := int(1e0)

	for i := 0; i < 1000; i++ {
		subj, clip := make([][]*Point, 0), make([][]*Point, 0)

		// Generate random subject and clip polygons ...
		subj = append(subj, RandomPoly(640*scale, 480*scale, 100))
		clip = append(clip, RandomPoly(640*scale, 480*scale, 100))

		//subjArea := AreaCombined(subj)
		//clipArea := AreaCombined(clip)

		c := NewClipper()
		pft := EvenOdd

		clipTypes := map[string]ClipType{"intersection": Intersection, "union": Union, "xor": Xor}
		areas := make(map[string]float64)
		// Load the polygons into Clipper and execute the boolean clip op ...
		c.AddPolygons(subj, Subject)
		c.AddPolygons(clip, Clip)

		for clipType, ct := range clipTypes {
			solution, result := c.Execute(ct, pft, pft)
			if !result {
				panic("Problem!")
			}
			areas[clipType] = AreaCombined(solution)
		}

		if different(areas["union"], areas["intersection"]+areas["xor"]) {
			t.Logf("%v\t%10.1f%10.1f\tFail",i,areas["union"], 
				areas["intersection"]+areas["xor"])
			t.FailNow()
		} else {
			t.Logf("%v\t%10.1f%10.1f\tPass",i,areas["union"], 
				areas["intersection"]+areas["xor"])
		}
	}
}

func different(a, b float64) bool {
	if math.Abs(a-b)/b > 0.01 {
		return true
	} else {
		return false
	}
}
