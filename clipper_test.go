package clipper

import (
	"fmt"
	"github.com/ctessum/carto"
	"github.com/twpayne/gogeom/geom"
	"image/color"
	"math"
	"math/rand"
	"os"
	"testing"
	"time"
)

func init() {
	rand.Seed(time.Now().UTC().UnixNano())
}

func RandomPoly(maxWidth, maxHeight, vertCnt int) Path {
	result := Path(make([]*IntPoint, vertCnt))
	for i := 0; i < vertCnt-1; i++ {
		result[i] = &IntPoint{cInt(rand.Intn(maxWidth)), cInt(rand.Intn(maxHeight))}
	}
	result[vertCnt-1] = result[0]
	return result
}

func TestRandom(t *testing.T) {
	scale := int(1e0)

	for i := 0; i < 1000; i++ {
		fmt.Println("--------------------",i,"-----------------")
		subj, clip := Paths(make([]Path, 0)), Paths(make([]Path, 0))

		// Generate random subject and clip polygons ...
		subj = append(subj, RandomPoly(640*scale, 480*scale, 4))
		clip = append(clip, RandomPoly(640*scale, 480*scale, 4))

		//subjArea := AreaCombined(subj)
		//clipArea := AreaCombined(clip)

		c := NewClipper(ioNone)
		pft := pftEvenOdd

		clipTypes := map[string]ClipType{"intersection": ctIntersection, "union": ctUnion, "xor": ctXor}
		areas := make(map[string]float64)
		// Load the polygons into Clipper and execute the boolean clip op ...
		c.AddPaths(subj, ptSubject, true)
		c.AddPaths(clip, ptClip, true)
		fmt.Println("Subj",subj)
		fmt.Println("clip",clip)

		for clipType, ct := range clipTypes {
			fmt.Println("Running" + clipType)
			solution, result := c.Execute1(ct, pft, pft)
			if !result {
				panic("Problem!")
			}
			fmt.Println("Finished Running" + clipType)
			areas[clipType] = AreaCombined(solution)
			fmt.Println("solution",solution)

			if clipType == "intersection" {
				f, err := os.Create(fmt.Sprintf("%v.png", i))
				if err != nil {
					panic(err)
				}
				carto.DrawShapes(f,
					[]color.NRGBA{{0, 0, 0, 255}, {0, 0, 0, 255}, {0, 0, 0, 255}},
					[]color.NRGBA{{255, 0, 0, 127}, {0, 255, 0, 127},
						{0, 0, 0, 200}},
					5, 0, clipper2geom(subj), clipper2geom(clip),
					clipper2geom(solution))
				f.Close()
			}
		}

		if different(areas["union"], areas["intersection"]+areas["xor"]) {
			t.Logf("%v\t%10.1f%10.1f\tFail", i, areas["union"],
				areas["intersection"]+areas["xor"])
			t.FailNow()
		} else {
			t.Logf("%v\t%10.1f%10.1f\tPass", i, areas["union"],
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
func clipper2geom(data []Path) geom.T {
	var out geom.T
	var temp geom.Polygon
	temp.Rings = make([][]geom.Point, len(data))
	for i, r := range data {
		temp.Rings[i] = make([]geom.Point, len(r))
		for j, p := range r {
			temp.Rings[i][j] = *geom.NewPoint(float64(p.X),
				float64(p.Y))
		}
	}
	if len(temp.Rings) != 0 {
		out = temp
	}
	return out
}
