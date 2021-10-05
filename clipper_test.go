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
		result[i] = &IntPoint{CInt(rand.Intn(maxWidth)), CInt(rand.Intn(maxHeight))}
	}
	result[vertCnt-1] = result[0]
	return result
}

func runOps(t *testing.T, testName string, subj, clip Paths) {
	// fmt.Printf("--------------------%v-----------------\n", testName)
	//subjArea := AreaCombined(subj)
	//clipArea := AreaCombined(clip)

	c := NewClipper(IoNone)
	pft := PftEvenOdd

	clipTypes := map[string]ClipType{"intersection": CtIntersection, "union": CtUnion, "xor": CtXor}
	areas := make(map[string]float64)
	// Load the polygons into Clipper and execute the boolean clip op ...
	c.AddPaths(subj, PtSubject, true)
	c.AddPaths(clip, PtClip, true)
	// fmt.Println("Subj", subj)
	// fmt.Println("clip", clip)

	var subjGeom, clipGeom, solutionGeom geom.Polygon
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

func clipper2geom(data []Path) geom.Polygon {
	var out geom.Polygon
	var temp geom.Polygon = make([]geom.Path, len(data))
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
		c := NewClipper(IoNone)
		pft := PftEvenOdd

		clipName := []string{"intersection", "union", "xor"}
		clipTypes := []ClipType{CtIntersection, CtUnion, CtXor}
		// Load the polygons into Clipper and execute the boolean clip op ...
		c.AddPaths(v.subj, PtSubject, true)
		c.AddPaths(v.clip, PtClip, true)
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

var testCasesPftNonZero = []testCase{
	{
		subj: Paths{
			{
				{X: 53000000000, Y: 180000000000},
				{X: 68000000000, Y: 200000000000},
				{X: 44000000000, Y: 199000000000},
			},
			{
				{X: 65000000000, Y: 160000000000},
				{X: 58000000000, Y: 189000000000},
				{X: 30000000000, Y: 190000000000},
			},
			{
				{X: 61000000000, Y: 189000000000},
				{X: 52000000000, Y: 195000000000},
				{X: 48000000000, Y: 187000000000},
			},
		},
		clip: Paths{},
		solutions: []Paths{
			// union
			{{{58426086957, 187234782609}, {59586956522, 188782608696}, {X: 61000000000, Y: 189000000000}, {X: 60166666667, Y: 189555555556}, {X: 68000000000, Y: 200000000000}, {X: 44000000000, Y: 199000000000}, {X: 48577437859, Y: 189336520076}, {X: 30000000000, Y: 190000000000}, {X: 65000000000, Y: 160000000000}}},
		},
	},
	{
		subj: Paths{
			{
				{X: 10606601714, Y: 115606601718},
				{X: 0, Y: 112071067812},
				{X: 17071067809, Y: 95000000001},
			},
			{
				{X: 42022005724, Y: 111666666668},
				{X: 21415404007, Y: 98131132763},
				{X: 28486471819, Y: 91060064951},
			},
			{
				{X: 9999999996, Y: 102071067814},
				{X: 17071067808, Y: 95000000002},
				{X: 27071067807, Y: 95000000002},
			},

			{
				{X: 20606601713, Y: 105606601717},
				{X: 7071067809, Y: 95000000000},
				{X: 17071067808, Y: 95000000000},
			},
			{
				{X: 58083558727, Y: 111818317094},
				{X: 36870355293, Y: 104747249284},
				{X: 51012490915, Y: 90605113660},
			},
			{
				{X: 5606601714, Y: 110606601717},
				{X: 2071067809, Y: 90000000000},
				{X: 12071067808, Y: 90000000000},
			},
		},
		clip: Paths{},
		solutions: []Paths{
			// union
			{{{28486471819, 91060064951}, {37236471817, 104381132760}, {51012490915, 90605113660}, {58083558727, 111818317094}, {37647005723, 105006132761}, {42022005724, 111666666668}, {21415404007, 98131132763}, {22761423751, 96785113019}, {18284271244, 98639610308}, {20606601713, 105606601717}, {15097873315, 101289898106}, {10606601714, 115606601718}, {0, 112071067812}, {4999999997, 107071067814}, {2071067809, 90000000000}, {12071067808, 90000000000}, {10502525313, 95000000000}, {17071067808, 95000000000}, {17071067809, 95000000002}, {24546536768, 95000000002}}, {{9825432008, 97158351806}, {7499999995, 104571067817}, {9999999997, 102071067813}, {12677669527, 99393398283}}},
		},
	},
}

func TestExecute1_Union_PftNonZero(t *testing.T) {
	for i, v := range testCasesPftNonZero {
		testName := fmt.Sprintf("fix%03d", i)
		c := NewClipper(IoNone)
		pft := PftNonZero

		clipName := []string{"union"}
		clipTypes := []ClipType{CtUnion}
		c.AddPaths(v.subj, PtSubject, true)
		c.AddPaths(v.clip, PtClip, true)

		for j, ct := range clipTypes {
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
