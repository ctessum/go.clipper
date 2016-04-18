package clipper

import (
	"fmt"
	"image/color"
	"math"
	"math/rand"
	"os"
	"testing"
	"time"

	"github.com/ctessum/geom"
	"github.com/ctessum/geom/carto"
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

func runOps(t *testing.T, testName string, subj, clip Paths) {
	// fmt.Printf("--------------------%v-----------------\n", testName)
	//subjArea := AreaCombined(subj)
	//clipArea := AreaCombined(clip)

	c := NewClipper(ioNone)
	pft := pftEvenOdd

	clipTypes := map[string]ClipType{"intersection": ctIntersection, "union": ctUnion, "xor": ctXor}
	areas := make(map[string]float64)
	// Load the polygons into Clipper and execute the boolean clip op ...
	c.AddPaths(subj, ptSubject, true)
	c.AddPaths(clip, ptClip, true)
	// fmt.Println("Subj", subj)
	// fmt.Println("clip", clip)

	var subjGeom, clipGeom, solutionGeom geom.T
	for clipType, ct := range clipTypes {
		// fmt.Println("Running " + clipType)
		solution, ok := c.Execute1(ct, pft, pft)
		if !ok {
			t.Fatal("Problem!")
		}
		// fmt.Println("Finished Running" + clipType)
		areas[clipType] = AreaCombined(solution)
		// fmt.Println("solution", solution)

		if clipType == "intersection" {
			subjGeom, clipGeom, solutionGeom = clipper2geom(subj), clipper2geom(clip), clipper2geom(solution)
		}
	}

	if different(areas["union"], areas["intersection"]+areas["xor"]) {
		t.Logf("writing file %v.png", testName)
		f, err := os.Create(fmt.Sprintf("%v.png", testName))
		if err != nil {
			t.Fatal(err)
		}
		carto.DrawShapes(f,
			[]color.NRGBA{{0, 0, 0, 255}, {0, 0, 0, 255}, {0, 0, 0, 255}},
			[]color.NRGBA{{255, 0, 0, 127}, {0, 255, 0, 127},
				{0, 0, 0, 200}},
			5, 0, subjGeom, clipGeom, solutionGeom)
		f.Close()
		t.Logf("%v\t%10.1f%10.1f\tFail", testName, areas["union"], areas["intersection"]+areas["xor"])
		t.FailNow()
	} else {
		t.Logf("%v\t%10.1f%10.1f\tPass", testName, areas["union"], areas["intersection"]+areas["xor"])
	}
}

func TestRandom(t *testing.T) {
	scale := int(1e0)

	for i := 0; i < 1000; i++ {
		// Generate random subject and clip polygons ...
		subj, clip := Paths(make([]Path, 0)), Paths(make([]Path, 0))
		subj = append(subj, RandomPoly(640*scale, 480*scale, 4))
		clip = append(clip, RandomPoly(640*scale, 480*scale, 4))
		runOps(t, fmt.Sprintf("rnd%03d", i), subj, clip)
	}
}

func different(a, b float64) bool {
	if math.Abs(a-b)/b > 0.01 {
		return true
	}
	return false
}

func clipper2geom(data []Path) geom.T {
	var out geom.T
	var temp geom.Polygon = make([][]geom.Point, len(data))
	for i, r := range data {
		temp[i] = make([]geom.Point, len(r))
		for j, p := range r {
			temp[i][j] = *geom.NewPoint(float64(p.X), float64(p.Y))
		}
	}
	if len(temp) != 0 {
		out = temp
	}
	return out
}

type testCase struct {
	subj, clip Paths
	solutions  []Paths
}

var testCases = []testCase{
	testCase{
		subj: Paths{{{343, 243}, {97, 337}, {413, 41}, {343, 243}}},
		clip: Paths{{{319, 178}, {68, 48}, {183, 197}, {319, 178}}},
		solutions: []Paths{
			Paths{{{319, 178}, {258, 187}, {285, 161}}},                                             // intersection
			Paths{{{343, 243}, {97, 337}, {258, 187}, {183, 197}, {68, 48}, {285, 161}, {413, 41}}}, // union
			Paths{{{343, 243}, {97, 337}, {258, 187}, {319, 178}, {285, 161}, {413, 41}}, // xor
				{{285, 161}, {258, 187}, {183, 197}, {68, 48}}},
		},
	},
	testCase{
		subj: Paths{{{460, 233}, {98, 400}, {147, 32}, {460, 233}}},
		clip: Paths{{{104, 33}, {494, 436}, {613, 347}, {104, 33}}},
		solutions: []Paths{
			Paths{{{442, 241}, {348, 285}, {142, 72}, {144, 57}}},                                                                      // intersection
			Paths{{{460, 233}, {442, 241}, {613, 347}, {494, 436}, {348, 285}, {98, 400}, {142, 72}, {104, 33}, {144, 57}, {147, 32}}}, // union
			Paths{{{460, 233}, {442, 241}, {613, 347}, {494, 436}, {348, 285}, {442, 241}, {144, 57}, {147, 32}}, // xor
				{{144, 57}, {142, 72}, {348, 285}, {98, 400}, {142, 72}, {104, 33}}},
		},
	},
	testCase{
		subj: Paths{{{615, 17}, {597, 282}, {232, 151}, {615, 17}}},
		clip: Paths{{{392, 414}, {177, 230}, {67, 230}, {392, 414}}},
		solutions: []Paths{
			Paths{}, // intersection
			Paths{{{392, 414}, {67, 230}, {177, 230}}, // union
				{{597, 282}, {232, 151}, {615, 17}}},
			Paths{{{392, 414}, {67, 230}, {177, 230}}, // xor
				{{597, 282}, {232, 151}, {615, 17}}},
		},
	},
}

func TestExecute1(t *testing.T) {
	for i, v := range testCases {
		testName := fmt.Sprintf("fix%03d", i)
		// fmt.Printf("--------------------%v-----------------\n", testName)
		c := NewClipper(ioNone)
		pft := pftEvenOdd

		clipName := []string{"intersection", "union", "xor"}
		clipTypes := []ClipType{ctIntersection, ctUnion, ctXor}
		// Load the polygons into Clipper and execute the boolean clip op ...
		c.AddPaths(v.subj, ptSubject, true)
		c.AddPaths(v.clip, ptClip, true)
		// fmt.Println("Subj", v.subj)
		// fmt.Println("clip", v.clip)

		for j, ct := range clipTypes {
			// fmt.Println("Running " + clipName[j])
			solution, ok := c.Execute1(ct, pft, pft)
			if !ok {
				t.Errorf("error running Execute1 %v: ok=false", clipName[j])
			}
			if got, want := solution.String(), v.solutions[j].String(); got != want {
				t.Errorf("bad %v %v solution, got=%v, want=%v", testName, clipName[j], got, want)
			}
		}
	}
}
