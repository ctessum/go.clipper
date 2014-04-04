//===============================================================================
//                                                                              //
// Author    :  Angus Johnson                                                   //
// Version   :  5.1.6(b)                                                        //
// Date      :  1 June 2013                                                     //
// Website   :  http://www.angusj.com                                           //
// Copyright :  Angus Johnson 2010-2013                                         //
//                                                                              //
// License:                                                                     //
// Use, modification & distribution is subject to Boost Software License Ver 1. //
// http://www.boost.org/LICENSE_1_0.txt                                         //
//                                                                              //
// Attributions:                                                                //
// The code in this library is an extension of Bala Vatti's clipping algorithm: //
// "A generic solution to polygon clipping"                                     //
// Communications of the ACM, Vol 35, Issue 7 (July 1992) PP 56-63.             //
// http://portal.acm.org/citation.cfm?id=129906                                 //
//                                                                              //
// Computer graphics && geometric modeling: implementation && algorithms      //
// By Max K. Agoston                                                            //
// Springer; 1 edition (January 4, 2005)                                        //
// http://books.google.com/books?q=vatti+clipping+agoston                       //
//                                                                              //
// See also:                                                                    //
// "Polygon Offsetting by Computing Winding Numbers"                            //
// Paper no. DETC2005-85513 PP. 565-575                                         //
// ASME 2005 International Design Engineering Technical Conferences             //
// && Computers && Information in Engineering Conference (IDETC/CIE2005)      //
// September 24-28, 2005 , Long Beach, California, USA                          //
// http://www.me.berkeley.edu/~mcmains/pubs/DAC05OffsetPolygon.pdf              //
//                                                                              //
//===============================================================================

package clipper

import (
	"fmt"
	"math"
	"sort"
)

var (
	horizontal = math.Inf(-1)
)

type ClipType int

const (
	Intersection ClipType = iota
	Union
	Difference
	Xor
)

type PolyType int

const (
	Subject PolyType = iota
	Clip
)

type PolyFillType int

const (
	EvenOdd PolyFillType = iota
	NonZero
	Positive
	Negative
)

type JoinType int

const (
	SquareJoin JoinType = iota
	RoundJoin
	MiterJoin
)

type EndType int

const (
	ClosedEnd EndType = iota
	ButtEnd
	SquareEnd
	RoundEnd
)

type edgeSide int

const (
	leftEdge edgeSide = iota
	rightEdge
)

type protects int

const (
	protectsNeither protects = iota
	protectsLeft
	protectsRight
	protectsBoth
)

type direction int

const (
	leftToRight direction = iota
	rightToLeft
)

type Point struct{ X, Y int }
type FloatPoint struct{ X, Y float64 }
type rect struct{ left, top, right, bottom int }

type localMinima struct {
	Y                     int
	leftBound, rightBound *edge
	nextLm                *localMinima
}

func newLocalMinima(y int, leftBound, rightBound *edge) *localMinima {
	out := new(localMinima)
	out.Y = y
	out.leftBound, out.rightBound = leftBound, rightBound
	return out
}

type scanbeam struct {
	Y      int
	nextSb *scanbeam
}

func (self *scanbeam) String() string {
	s := "nil"
	if self.nextSb != nil {
		s = "<obj>"
	}
	return fmt.Sprintf("(y:%d, nextSb:%s)", self.Y, s)
}

type intersectNode struct {
	e1, e2 *edge
	pt     *Point
	nextIn *intersectNode
}

type outPt struct {
	idx            int
	pt             *Point
	prevOp, nextOp *outPt
}

func newOutPt(idx int, pt *Point) *outPt {
	out := new(outPt)
	out.idx = idx
	out.pt = pt
	return out
}

type outRec struct {
	idx       int
	bottomPt  *outPt
	isHole    bool
	FirstLeft *outRec
	pts       *outPt
	polyNode  *PolyNode
}

func newOutRec(idx int) *outRec {
	out := new(outRec)
	out.idx = idx
	return out
}

type joinRec struct {
	pt1a, pt1b *Point
	poly1Idx   int
	pt2a, pt2b *Point
	poly2Idx   int
}

type horzJoin struct {
	edge           *edge
	savedIdx       int
	prevHj, nextHj *horzJoin
}

func newHorzJoin(e *edge, idx int) *horzJoin {
	hj := new(horzJoin)
	hj.edge = e
	hj.savedIdx = idx
	return hj
}

//===============================================================================
// Unit global functions ...
//===============================================================================

func IntsToPoints(ints []int) []*Point {
	result := make([]*Point, len(ints)/2)
	for i := 0; i < len(ints); i += 2 {
		result[i] = &Point{ints[i], ints[i+1]}
	}
	return result
}

// see http://www.mathopenref.com/coordpolygonarea2.html
func Area(polygon []*Point) int {
	highI := len(polygon) - 1
	A := (polygon[highI].X + polygon[0].X) * (polygon[0].Y - polygon[highI].Y)
	for i := 0; i < highI; i++ {
		A += (polygon[i].X + polygon[i+1].X) * (polygon[i+1].Y - polygon[i].Y)
	}
	return A / 2
}

func orientation(polygon []*Point) bool {
	return Area(polygon) > 0.0
}

//===============================================================================
// PolyNode & PolyTree classes (+ ancilliary functions)
//===============================================================================

// Node of PolyTree
type PolyNode struct {
	Contour           []*Point
	Childs            []*PolyNode
	Parent            *PolyNode
	Index, ChildCount int
}

func (self *PolyNode) IsHole() bool {
	result := true
	p := self.Parent
	for {
		if p == nil {
			break
		}
		result = !result
		p = p.Parent
	}
	return result
}

func (self *PolyNode) GetNext() *PolyNode {
	if self.ChildCount > 0 {
		return self.Childs[0]
	} else {
		return self.getNextSiblingUp()
	}
}

func (self *PolyNode) addChild(node *PolyNode) {
	self.Childs = append(self.Childs, node)
	node.Index = self.ChildCount
	node.Parent = self
	self.ChildCount += 1
}

func (self *PolyNode) getNextSiblingUp() *PolyNode {
	if self.Parent == nil {
		return nil
	} else if self.Index == self.Parent.ChildCount-1 {
		return self.Parent.getNextSiblingUp()
	} else {
		return self.Parent.Childs[self.Index+1]
	}
}

// Container for PolyNodes
type PolyTree struct {
	PolyNode
	allNodes []*PolyNode
}

func (self *PolyTree) toPolyNode() *PolyNode {
	node := new(PolyNode)
	node.Contour = self.Contour
	node.Childs = self.Childs
	node.Parent = self.Parent
	node.Index, node.ChildCount = self.Index, self.ChildCount
	return node
}

func (self *PolyTree) Clear() {
	self.allNodes = nil
	self.Childs = nil
	self.ChildCount = 0
}

func (self *PolyTree) GetFirst() *PolyNode {
	if self.ChildCount > 0 {
		return self.Childs[0]
	} else {
		return nil
	}
}

func (self *PolyTree) Total() int {
	return len(self.allNodes)
}

// Internal function for PolyTreeToPolygons()
func addPolyNodeToPolygons(polynode *PolyNode, polygons [][]*Point) {
	if len(polynode.Contour) > 0 {
		polygons = append(polygons, polynode.Contour)
	}
	for i := 0; i < polynode.ChildCount; i++ {
		addPolyNodeToPolygons(polynode.Childs[i], polygons)
	}
}

func PolyTreeToPolygons(polyTree *PolyTree) [][]*Point {
	result := make([][]*Point, 0)
	addPolyNodeToPolygons(polyTree.toPolyNode(), result)
	return result
}

//===============================================================================
// Edge class
//===============================================================================

type edge struct {
	Bot, Curr, Top, Delta                      *Point
	dx                                         float64
	polyType                                   PolyType
	side                                       edgeSide
	windDelta, windCnt, windCnt2, outIdx       int
	nextE, prevE, nextInLML                    *edge
	prevInAEL, nextInAEL, prevInSEL, nextInSEL *edge
}

func NewEdge() *edge {
	self := new(edge)
	self.Bot = new(Point)
	self.Curr = new(Point)
	self.Top = new(Point)
	self.Delta = new(Point)
	self.outIdx = -1
	self.polyType = Subject
	self.side = leftEdge
	return self
}

func (self *edge) String() string {
	return fmt.Sprintf("(%d,%d . %d,%d {dx:%0.2f} %i)",
		self.Bot.X, self.Bot.Y, self.Top.X, self.Top.Y, self.dx, self.outIdx)
}

//===============================================================================
// ClipperBase class (+ data structs & ancilliary functions)
//===============================================================================

func pointsEqual(pt1, pt2 *Point) bool {
	return (pt1.X == pt2.X) && (pt1.Y == pt2.Y)
}

func slopesEqual(pt1, pt2, pt3, pt4 *Point) bool {
	if pt4 == nil {
		return (pt1.Y-pt2.Y)*(pt2.X-pt3.X) == (pt1.X-pt2.X)*(pt2.Y-pt3.Y)
	} else {
		return (pt1.Y-pt2.Y)*(pt3.X-pt4.X) == (pt1.X-pt2.X)*(pt3.Y-pt4.Y)
	}
}

func slopesEqual2(e1, e2 *edge) bool {
	return e1.Delta.Y*e2.Delta.X == e1.Delta.X*e2.Delta.Y
}

func setDx(e *edge) {
	e.Delta = &Point{e.Top.X - e.Bot.X, e.Top.Y - e.Bot.Y}
	if e.Delta.Y == 0 {
		e.dx = horizontal
	} else {
		e.dx = float64(e.Delta.X) / float64(e.Delta.Y)
	}
}

func swapSides(e1, e2 *edge) {
	side := e1.side
	e1.side = e2.side
	e2.side = side
}

func swapPolyIndexes(e1, e2 *edge) {
	idx := e1.outIdx
	e1.outIdx = e2.outIdx
	e2.outIdx = idx
}

func initEdge(e, eNext, ePrev *edge, pt *Point, polyType PolyType) {
	e.nextE = eNext
	e.prevE = ePrev
	e.Curr = pt
	if e.Curr.Y >= e.nextE.Curr.Y {
		e.Bot = e.Curr
		e.Top = e.nextE.Curr
		e.windDelta = 1
	} else {
		e.Top = e.Curr
		e.Bot = e.nextE.Curr
		e.windDelta = -1
	}
	setDx(e)
	e.outIdx = -1
	e.polyType = polyType
}

func swapX(e *edge) {
	e.Curr = &Point{e.Top.X, e.Curr.Y}
	e.Top = &Point{e.Bot.X, e.Top.Y}
	e.Bot = &Point{e.Curr.X, e.Bot.Y}
}

type ClipperBase struct {
	edgeList      [][]*edge    // 2D array
	localMinList  *localMinima // single-linked list of LocalMinima
	currentLocMin *localMinima
}

func (self *ClipperBase) insertLocalMinima(lm *localMinima) {
	if self.localMinList == nil {
		self.localMinList = lm
	} else if lm.Y >= self.localMinList.Y {
		lm.nextLm = self.localMinList
		self.localMinList = lm
	} else {
		tmp := self.localMinList
		for tmp.nextLm != nil && lm.Y < tmp.nextLm.Y {
			tmp = tmp.nextLm
		}
		lm.nextLm = tmp.nextLm
		tmp.nextLm = lm
	}
}

func (self *ClipperBase) addBoundsToLML(e *edge) *edge {
	e.nextInLML = nil
	e = e.nextE
	for {
		if e.dx == horizontal {
			if (e.nextE.Top.Y < e.Top.Y) && (e.nextE.Bot.X > e.prevE.Bot.X) {
				break
			}
			if e.Top.X != e.prevE.Bot.X {
				swapX(e)
			}
			e.nextInLML = e.prevE
		} else if e.Bot.Y == e.prevE.Bot.Y {
			break
		} else {
			e.nextInLML = e.prevE
		}
		e = e.nextE
	}

	var lm *localMinima
	if e.dx == horizontal {
		if e.Bot.X != e.prevE.Bot.X {
			swapX(e)
		}
		lm = newLocalMinima(e.prevE.Bot.Y, e.prevE, e)
	} else if e.dx < e.prevE.dx {
		lm = newLocalMinima(e.prevE.Bot.Y, e.prevE, e)
	} else {
		lm = newLocalMinima(e.prevE.Bot.Y, e, e.prevE)
	}
	lm.leftBound.side = leftEdge
	lm.rightBound.side = rightEdge
	self.insertLocalMinima(lm)
	for {
		if e.nextE.Top.Y == e.Top.Y && e.nextE.dx != horizontal {
			break
		}
		e.nextInLML = e.nextE
		e = e.nextE
		if e.dx == horizontal && e.Bot.X != e.prevE.Top.X {
			swapX(e)
		}
	}
	return e.nextE
}

func (self *ClipperBase) resetBase() {
	lm := self.localMinList
	if lm != nil {
		self.currentLocMin = lm
	}
	for lm != nil {
		e := lm.leftBound
		for e != nil {
			e.Curr = e.Bot
			e.side = leftEdge
			e.outIdx = -1
			e = e.nextInLML
		}
		e = lm.rightBound
		for e != nil {
			e.Curr = e.Bot
			e.side = rightEdge
			e.outIdx = -1
			e = e.nextInLML
		}
		lm = lm.nextLm
	}
}

func (self *ClipperBase) AddPolygon(polygon []*Point, polyType PolyType) bool {
	ln := len(polygon)
	if ln < 3 {
		return false
	}
	pg := make([]*Point, len(polygon))
	copy(pg, polygon)
	j := 0
	// remove duplicate points && co-linear points
	for i := 1; i < len(polygon); i++ {
		if pointsEqual(pg[j], polygon[i]) {
			continue
		} else if (j > 0) && slopesEqual(pg[j-1], pg[j], polygon[i], nil) {
			if pointsEqual(pg[j-1], polygon[i]) {
				j -= 1
			}
		} else {
			j += 1
		}
		pg[j] = polygon[i]
	}
	if j < 2 {
		return false
	}
	// remove duplicate points && co-linear edges at the loop around
	// of the start && end coordinates ...
	ln = j + 1
	for ln > 2 {
		if pointsEqual(pg[j], pg[0]) {
			j -= 1
		} else if pointsEqual(pg[0], pg[1]) ||
			slopesEqual(pg[j], pg[0], pg[1], nil) {
			pg[0] = pg[j]
			j -= 1
		} else if slopesEqual(pg[j-1], pg[j], pg[0], nil) {
			j -= 1
		} else if slopesEqual(pg[0], pg[1], pg[2], nil) {
			for i := 2; i < j+1; i++ {
				pg[i-1] = pg[i]
			}
			j -= 1
		} else {
			break
		}
		ln -= 1
	}
	if ln < 3 {
		return false
	}
	edges := make([]*edge, 0)
	for i := 0; i < ln; i++ {
		edges = append(edges, NewEdge())
	}
	edges[0].Curr = pg[0]
	initEdge(edges[ln-1], edges[0], edges[ln-2], pg[ln-1], polyType)
	for i := ln - 2; i > 0; i-- {
		initEdge(edges[i], edges[i+1], edges[i-1], pg[i], polyType)
	}
	initEdge(edges[0], edges[1], edges[ln-1], pg[0], polyType)
	e := edges[0]
	eHighest := e
	for {
		e.Curr = e.Bot
		if e.Top.Y < eHighest.Top.Y {
			eHighest = e
		}
		e = e.nextE
		if e == edges[0] {
			break
		}
	}
	// make sure eHighest is positioned so the following loop works safely ...
	if eHighest.windDelta > 0 {
		eHighest = eHighest.nextE
	}
	if eHighest.dx == horizontal {
		eHighest = eHighest.nextE
	}
	// finally insert each local minima ...
	e = eHighest
	for {
		e = self.addBoundsToLML(e)
		if e == eHighest {
			break
		}
	}
	self.edgeList = append(self.edgeList, edges)
	return true
}

func (self *ClipperBase) AddPolygons(polygons [][]*Point, polyType PolyType) bool {
	result := false
	for _, p := range polygons {
		if self.AddPolygon(p, polyType) {
			result = true
		} else {
			break
		}
	}
	return result
}

func (self *ClipperBase) clearBase() {
	self.edgeList = make([][]*edge, 0)
	self.localMinList = nil
	self.currentLocMin = nil
}

func (self *ClipperBase) popLocalMinima() {
	if self.currentLocMin != nil {
		self.currentLocMin = self.currentLocMin.nextLm
	}
}

//===============================================================================
// Clipper class (+ data structs & ancilliary functions)
//===============================================================================
func intersectPoint(edge1, edge2 *edge) (*Point, bool) {
	var x, y int
	if slopesEqual2(edge1, edge2) {
		if edge2.Bot.Y > edge1.Bot.Y {
			y = edge2.Bot.Y
		} else {
			y = edge1.Bot.Y
		}
		return &Point{0, y}, false
	}
	if edge1.dx == 0 {
		x = edge1.Bot.X
		if edge2.dx == horizontal {
			y = edge2.Bot.Y
		} else {
			b2 := float64(edge2.Bot.Y) - float64(edge2.Bot.X)/edge2.dx
			y = round(float64(x)/edge2.dx + b2)
		}
	} else if edge2.dx == 0 {
		x = edge2.Bot.X
		if edge1.dx == horizontal {
			y = edge1.Bot.Y
		} else {
			b1 := float64(edge1.Bot.Y) - float64(edge1.Bot.X)/edge1.dx
			y = round(float64(x)/edge1.dx + b1)
		}
	} else {
		b1 := float64(edge1.Bot.X) - float64(edge1.Bot.Y)*edge1.dx
		b2 := float64(edge2.Bot.X) - float64(edge2.Bot.Y)*edge2.dx
		m := (b2 - b1) / (edge1.dx - edge2.dx)
		y = round(m)
		if math.Abs(edge1.dx) < math.Abs(edge2.dx) {
			x = round(edge1.dx*m + b1)
		} else {
			x = round(edge2.dx*m + b2)
		}
	}
	if (y < edge1.Top.Y) || (y < edge2.Top.Y) {
		if edge1.Top.Y > edge2.Top.Y {
			return edge1.Top, topX(edge2, edge1.Top.Y) < edge1.Top.X
		} else {
			return edge2.Top, topX(edge1, edge2.Top.Y) > edge2.Top.X
		}
	} else {
		return &Point{x, y}, true
	}
}

func topX(e *edge, currentY int) int {
	if currentY == e.Top.Y {
		return e.Top.X
	} else if e.Top.X == e.Bot.X {
		return e.Bot.X
	} else {
		return e.Bot.X + round(e.dx*float64(currentY-e.Bot.Y))
	}
}

func e2InsertsBeforeE1(e1, e2 *edge) bool {
	if e2.Curr.X == e1.Curr.X {
		if e2.Top.Y > e1.Top.Y {
			return e2.Top.X < topX(e1, e2.Top.Y)
		}
		return e1.Top.X > topX(e2, e1.Top.Y)
	} else {
		return e2.Curr.X < e1.Curr.X
	}
}

func isMinima(e *edge) bool {
	return (e != nil) && (e.prevE.nextInLML != e) && (e.nextE.nextInLML != e)
}

func isMaxima(e *edge, y int) bool {
	return (e != nil) && (e.Top.Y == y) && (e.nextInLML == nil)
}

func isIntermediate(e *edge, y int) bool {
	return e.Top.Y == y && e.nextInLML != nil
}

func getMaximaPair(e *edge) *edge {
	if !isMaxima(e.nextE, e.Top.Y) || e.nextE.Top.X != e.Top.X {
		return e.prevE
	} else {
		return e.nextE
	}
}

func getnextInAEL(e *edge, dir direction) *edge {
	if dir == leftToRight {
		return e.nextInAEL
	} else {
		return e.prevInAEL
	}
}

func protectLeft(val bool) protects {
	if val {
		return protectsBoth
	} else {
		return protectsRight
	}
}

func protectRight(val bool) protects {
	if val {
		return protectsBoth
	} else {
		return protectsLeft
	}
}

func getDx(pt1, pt2 *Point) float64 {
	if pt1.Y == pt2.Y {
		return horizontal
	} else {
		return float64(pt2.X-pt1.X) / float64(pt2.Y-pt1.Y)
	}
}

func param1RightOfParam2(outRec1, outRec2 *outRec) bool {
	for outRec1 != nil {
		outRec1 = outRec1.FirstLeft
		if outRec1 == outRec2 {
			return true
		}
	}
	return false
}

func firstParamIsbottomPt(btmPt1, btmPt2 *outPt) bool {
	p := btmPt1.prevOp
	for pointsEqual(p.pt, btmPt1.pt) && (p != btmPt1) {
		p = p.prevOp
	}
	dx1p := math.Abs(getDx(btmPt1.pt, p.pt))
	p = btmPt1.nextOp
	for pointsEqual(p.pt, btmPt1.pt) && (p != btmPt1) {
		p = p.nextOp
	}
	dx1n := math.Abs(getDx(btmPt1.pt, p.pt))

	p = btmPt2.prevOp
	for pointsEqual(p.pt, btmPt2.pt) && (p != btmPt2) {
		p = p.prevOp
	}
	dx2p := math.Abs(getDx(btmPt2.pt, p.pt))
	p = btmPt2.nextOp
	for pointsEqual(p.pt, btmPt2.pt) && (p != btmPt2) {
		p = p.nextOp
	}
	dx2n := math.Abs(getDx(btmPt2.pt, p.pt))
	return (dx1p >= dx2p && dx1p >= dx2n) || (dx1n >= dx2p && dx1n >= dx2n)
}

func getBottomPt(pp *outPt) *outPt {
	var dups *outPt
	p := pp.nextOp
	for p != pp {
		if p.pt.Y > pp.pt.Y {
			pp = p
			dups = nil
		} else if p.pt.Y == pp.pt.Y && p.pt.X <= pp.pt.X {
			if p.pt.X < pp.pt.X {
				dups = nil
				pp = p
			} else {
				if p.nextOp != pp && p.prevOp != pp {
					dups = p
				}
			}
		}
		p = p.nextOp
	}
	if dups != nil {
		for dups != p {
			if !firstParamIsbottomPt(p, dups) {
				pp = dups
			}
			dups = dups.nextOp
			for !pointsEqual(dups.pt, pp.pt) {
				dups = dups.nextOp
			}
		}
	}
	return pp
}

func getLowermostRec(outRec1, outRec2 *outRec) *outRec {
	var outPt1, outPt2 *outPt
	if outRec1.bottomPt == nil {
		outPt1 = getBottomPt(outRec1.pts)
	} else {
		outPt1 = outRec1.bottomPt
	}
	if outRec2.bottomPt == nil {
		outPt2 = getBottomPt(outRec2.pts)
	} else {
		outPt2 = outRec2.bottomPt
	}
	if outPt1.pt.Y > outPt2.pt.Y {
		return outRec1
	} else if outPt1.pt.Y < outPt2.pt.Y {
		return outRec2
	} else if outPt1.pt.X < outPt2.pt.X {
		return outRec1
	} else if outPt1.pt.X > outPt2.pt.X {
		return outRec2
	} else if outPt1.nextOp == outPt1 {
		return outRec2
	} else if outPt2.nextOp == outPt2 {
		return outRec1
	} else if firstParamIsbottomPt(outPt1, outPt2) {
		return outRec1
	} else {
		return outRec2
	}
}

func setHoleState(e *edge, outRec *outRec, polyOutList []*outRec) {
	isHole := false
	e2 := e.prevInAEL
	for e2 != nil {
		if e2.outIdx >= 0 {
			isHole = !isHole
			if outRec.FirstLeft == nil {
				outRec.FirstLeft = polyOutList[e2.outIdx]
			}
		}
		e2 = e2.prevInAEL
	}
	outRec.isHole = isHole
}

func pointCount(pts *outPt) int {
	if pts == nil {
		return 0
	}
	p := pts
	result := 0
	for {
		result++
		p = p.nextOp
		if p == pts {
			break
		}
	}
	return result
}

func pointIsVertex(pt *Point, outPts *outPt) bool {
	op := outPts
	for {
		if pointsEqual(op.pt, pt) {
			return true
		}
		op = op.nextOp
		if op == outPts {
			break
		}
	}
	return false
}

func reversePolyPtLinks(pp *outPt) {
	if pp == nil {
		return
	}
	pp1 := pp
	for {
		pp2 := pp1.nextOp
		pp1.nextOp = pp1.prevOp
		pp1.prevOp = pp2
		pp1 = pp2
		if pp1 == pp {
			break
		}
	}
}

func fixupOutPolygon(outRec *outRec) {
	var lastOK *outPt
	outRec.bottomPt = nil
	pp := outRec.pts
	for {
		if pp.prevOp == pp || pp.nextOp == pp.prevOp {
			outRec.pts = nil
			return
		}
		if pointsEqual(pp.pt, pp.nextOp.pt) ||
			slopesEqual(pp.prevOp.pt, pp.pt, pp.nextOp.pt, nil) {
			lastOK = nil
			pp.prevOp.nextOp = pp.nextOp
			pp.nextOp.prevOp = pp.prevOp
			pp = pp.prevOp
		} else if pp == lastOK {
			break
		} else {
			if lastOK == nil {
				lastOK = pp
			}
			pp = pp.nextOp
		}
	}
	outRec.pts = pp
}

func fixHoleLinkage(outrec *outRec) {
	if outrec.FirstLeft == nil ||
		(outrec.isHole != outrec.FirstLeft.isHole &&
			outrec.FirstLeft.pts != nil) {
		return
	}
	orfl := outrec.FirstLeft
	for orfl != nil &&
		(orfl.isHole == outrec.isHole || orfl.pts == nil) {
		orfl = orfl.FirstLeft
	}
	outrec.FirstLeft = orfl
}

func getOverlapSegment(pt1a, pt1b, pt2a, pt2b *Point) (*Point, *Point, bool) {
	// precondition: segments are co-linear
	var pt1, pt2 *Point
	if intAbs(pt1a.X-pt1b.X) > intAbs(pt1a.Y-pt1b.Y) {
		if pt1a.X > pt1b.X {
			tmp := pt1a
			pt1a = pt1b
			pt1b = tmp
		}
		if pt2a.X > pt2b.X {
			tmp := pt2a
			pt2a = pt2b
			pt2b = tmp
		}
		if pt1a.X > pt2a.X {
			pt1 = pt1a
		} else {
			pt1 = pt2a
		}
		if pt1b.X < pt2b.X {
			pt2 = pt1b
		} else {
			pt2 = pt2b
		}
		return pt1, pt2, pt1.X < pt2.X
	} else {
		if pt1a.Y < pt1b.Y {
			tmp := pt1a
			pt1a = pt1b
			pt1b = tmp
		}
		if pt2a.Y < pt2b.Y {
			tmp := pt2a
			pt2a = pt2b
			pt2b = tmp
		}
		if pt1a.Y < pt2a.Y {
			pt1 = pt1a
		} else {
			pt1 = pt2a
		}
		if pt1b.Y > pt2b.Y {
			pt2 = pt1b
		} else {
			pt2 = pt2b
		}
		return pt1, pt2, pt1.Y > pt2.Y
	}
}

func findSegment(outPt *outPt, pt1, pt2 *Point) (*outPt, *Point, *Point, bool) {
	if outPt == nil {
		return outPt, pt1, pt2, false
	}
	pt1a := pt1
	pt2a := pt2
	outPt2 := outPt
	for {
		if slopesEqual(pt1a, pt2a, outPt.pt, outPt.prevOp.pt) && slopesEqual(pt1a, pt2a, outPt.pt, nil) {
			pt1, pt2, overlap := getOverlapSegment(pt1a, pt2a, outPt.pt, outPt.prevOp.pt)
			if overlap {
				return outPt, pt1, pt2, true
			}
		}
		outPt = outPt.nextOp
		if outPt == outPt2 {
			return outPt, pt1, pt2, false
		}
	}
}

func pt3IsBetweenPt1AndPt2(pt1, pt2, pt3 *Point) bool {
	if pointsEqual(pt1, pt3) || pointsEqual(pt2, pt3) {
		return true
	} else if pt1.X != pt2.X {
		return (pt1.X < pt3.X) == (pt3.X < pt2.X)
	} else {
		return (pt1.Y < pt3.Y) == (pt3.Y < pt2.Y)
	}
}

func insertPolyPtBetween(outPt1, outPt2 *outPt, pt *Point) *outPt {
	if outPt1 == outPt2 {
		panic("JoinError")
	}
	result := newOutPt(outPt1.idx, pt)
	if outPt2 == outPt1.nextOp {
		outPt1.nextOp = result
		outPt2.prevOp = result
		result.nextOp = outPt2
		result.prevOp = outPt1
	} else {
		outPt2.nextOp = result
		outPt1.prevOp = result
		result.nextOp = outPt1
		result.prevOp = outPt2
	}
	return result
}

func PointOnLineSegment(pt, linePt1, linePt2 *Point) bool {
	return ((pt.X == linePt1.X) && (pt.Y == linePt1.Y)) ||
		((pt.X == linePt2.X) && (pt.Y == linePt2.Y)) ||
		(((pt.X > linePt1.X) == (pt.X < linePt2.X)) &&
			((pt.Y > linePt1.Y) == (pt.Y < linePt2.Y)) &&
			((pt.X-linePt1.X)*(linePt2.Y-linePt1.Y) ==
				(linePt2.X-linePt1.X)*(pt.Y-linePt1.Y)))
}

func PointOnPolygon(pt *Point, pp *outPt) bool {
	pp2 := pp
	for {
		if PointOnLineSegment(pt, pp2.pt, pp2.nextOp.pt) {
			return false
		}
		pp2 = pp2.nextOp
		if pp2 == pp {
			return false
		}
	}
}

func PointInPolygon(pt *Point, outPt *outPt) bool {
	result := false
	outPt2 := outPt
	for {
		if (((outPt2.pt.Y <= pt.Y) && (pt.Y < outPt2.prevOp.pt.Y)) ||
			((outPt2.prevOp.pt.Y <= pt.Y) && (pt.Y < outPt2.pt.Y))) &&
			(pt.X < (outPt2.prevOp.pt.X-outPt2.pt.X)*(pt.Y-outPt2.pt.Y)/
				(outPt2.prevOp.pt.Y-outPt2.pt.Y)+outPt2.pt.X) {
			result = !result
		}
		outPt2 = outPt2.nextOp
		if outPt2 == outPt {
			break
		}
	}
	return result
}

func poly2ContainsPoly1(outPt1, outPt2 *outPt) bool {
	pt := outPt1
	if PointOnPolygon(pt.pt, outPt2) {
		pt = pt.nextOp
		for pt != outPt1 && PointOnPolygon(pt.pt, outPt2) {
			pt = pt.nextOp
		}
		if pt == outPt1 {
			return true
		}
	}
	return PointInPolygon(pt.pt, outPt2)
}

func edgesAdjacent(inode *intersectNode) bool {
	return (inode.e1.nextInSEL == inode.e2) ||
		(inode.e1.prevInSEL == inode.e2)
}

func updateOutPtIdxs(outrec *outRec) {
	op := outrec.pts
	for {
		op.idx = outrec.idx
		op = op.prevOp
		if op == outrec.pts {
			break
		}
	}
}

type Clipper struct {
	ClipperBase
	ReverseSolution bool
	ForceSimple     bool

	polyOutList    []*outRec
	clipType       ClipType
	scanbeam       *scanbeam
	activeEdges    *edge
	sortedEdges    *edge
	intersectNodes *intersectNode
	clipFillType   PolyFillType
	subjFillType   PolyFillType
	executeLocked  bool
	usingPolyTree  bool
	joinList       []*joinRec
	horzJoins      *horzJoin
}

func NewClipper() *Clipper {
	self := new(Clipper)
	self.ReverseSolution = false
	self.ForceSimple = false
	self.polyOutList = make([]*outRec, 0)
	self.clipType = Intersection
	self.clipFillType = EvenOdd
	self.subjFillType = EvenOdd
	return self
}

func (self *Clipper) reset() {
	self.resetBase()
	self.scanbeam = nil
	self.polyOutList = make([]*outRec, 0)
	lm := self.localMinList
	for lm != nil {
		self.insertScanbeam(lm.Y)
		lm = lm.nextLm
	}
}

func (self *Clipper) Clear() {
	self.polyOutList = make([]*outRec, 0)
	self.clearBase()
}

func (self *Clipper) insertScanbeam(y int) {
	if self.scanbeam == nil {
		self.scanbeam = &scanbeam{Y: y}
	} else if y > self.scanbeam.Y {
		self.scanbeam = &scanbeam{y, self.scanbeam}
	} else {
		sb := self.scanbeam
		for sb.nextSb != nil && y <= sb.nextSb.Y {
			sb = sb.nextSb
		}
		if y == sb.Y {
			return
		}
		newSb := &scanbeam{y, sb.nextSb}
		sb.nextSb = newSb
	}
}

func (self *Clipper) popScanbeam() int {
	result := self.scanbeam.Y
	self.scanbeam = self.scanbeam.nextSb
	return result
}

func (self *Clipper) setWindingCount(edge *edge) {
	e := edge.prevInAEL
	for e != nil && e.polyType != edge.polyType {
		e = e.prevInAEL
	}
	if e == nil {
		edge.windCnt = edge.windDelta
		edge.windCnt2 = 0
		e = self.activeEdges
	} else if self.isEvenOddFillType(edge) {
		edge.windCnt = 1
		edge.windCnt2 = e.windCnt2
		e = e.nextInAEL
	} else {
		if e.windCnt*e.windDelta < 0 {
			if intAbs(e.windCnt) > 1 {
				if e.windDelta*edge.windDelta < 0 {
					edge.windCnt = e.windCnt
				} else {
					edge.windCnt = e.windCnt + edge.windDelta
				}
			} else {
				edge.windCnt = e.windCnt + e.windDelta + edge.windDelta
			}
		} else if (intAbs(e.windCnt) > 1) && (e.windDelta*edge.windDelta < 0) {
			edge.windCnt = e.windCnt
		} else if e.windCnt+edge.windDelta == 0 {
			edge.windCnt = e.windCnt
		} else {
			edge.windCnt = e.windCnt + edge.windDelta
		}
		edge.windCnt2 = e.windCnt2
		e = e.nextInAEL
	}
	// update windCnt2 ...
	if self.isEvenOddAltFillType(edge) {
		for e != edge {
			if edge.windCnt2 == 0 {
				edge.windCnt2 = 1
			} else {
				edge.windCnt2 = 0
			}
			e = e.nextInAEL
		}
	} else {
		for e != edge {
			edge.windCnt2 += e.windDelta
			e = e.nextInAEL
		}
	}
}

func (self *Clipper) isEvenOddFillType(edge *edge) bool {
	if edge.polyType == Subject {
		return self.subjFillType == EvenOdd
	} else {
		return self.clipFillType == EvenOdd
	}
}

func (self *Clipper) isEvenOddAltFillType(edge *edge) bool {
	if edge.polyType == Subject {
		return self.clipFillType == EvenOdd
	} else {
		return self.subjFillType == EvenOdd
	}
}

func (self *Clipper) isContributing(edge *edge) bool {
	var pft, pft2 PolyFillType
	if edge.polyType == Subject {
		pft = self.subjFillType
		pft2 = self.clipFillType
	} else {
		pft = self.clipFillType
		pft2 = self.subjFillType
	}
	if pft == EvenOdd || pft == NonZero {
		if intAbs(edge.windCnt) != 1 {
			return false
		}
	} else if pft == Positive {
		if edge.windCnt != 1 {
			return false
		}
	} else if pft == Negative {
		if edge.windCnt != -1 {
			return false
		}
	}

	if self.clipType == Intersection { //////////////////////
		if pft2 == EvenOdd || pft2 == NonZero {
			return edge.windCnt2 != 0
		} else if pft2 == Positive {
			return edge.windCnt2 > 0
		} else {
			return edge.windCnt2 < 0 // Negative
		}
	} else if self.clipType == Union { //////////////////////
		if pft2 == EvenOdd || pft2 == NonZero {
			return edge.windCnt2 == 0
		} else if pft2 == Positive {
			return edge.windCnt2 <= 0
		} else {
			return edge.windCnt2 >= 0 // Negative
		}
	} else if self.clipType == Difference { //////////////////////
		if edge.polyType == Subject {
			if pft2 == EvenOdd || pft2 == NonZero {
				return edge.windCnt2 == 0
			} else if pft2 == Positive {
				return edge.windCnt2 <= 0
			} else {
				return edge.windCnt2 >= 0
			}
		} else {
			if pft2 == EvenOdd || pft2 == NonZero {
				return edge.windCnt2 != 0
			} else if pft2 == Positive {
				return edge.windCnt2 > 0
			} else {
				return edge.windCnt2 < 0
			}
		}
	} else { // self._ClipType == ClipType.XOR:     //////////////////////
		return true
	}
}

func (self *Clipper) addEdgeToSEL(edge *edge) {
	if self.sortedEdges == nil {
		self.sortedEdges = edge
		edge.prevInSEL = nil
		edge.nextInSEL = nil
	} else {
		// add edge to front of list ...
		edge.nextInSEL = self.sortedEdges
		edge.prevInSEL = nil
		self.sortedEdges.prevInSEL = edge
		self.sortedEdges = nil
	}
}

func (self *Clipper) copyAELToSEL() {
	e := self.activeEdges
	self.sortedEdges = e
	for e != nil {
		e.prevInSEL = e.prevInAEL
		e.nextInSEL = e.nextInAEL
		e = e.nextInAEL
	}
}

func (self *Clipper) insertEdgeIntoAEL(edge *edge) {
	edge.prevInAEL = nil
	edge.nextInAEL = nil
	if self.activeEdges == nil {
		self.activeEdges = edge
	} else if e2InsertsBeforeE1(self.activeEdges, edge) {
		edge.nextInAEL = self.activeEdges
		self.activeEdges.prevInAEL = edge
		self.activeEdges = edge
	} else {
		e := self.activeEdges
		for e.nextInAEL != nil &&
			!e2InsertsBeforeE1(e.nextInAEL, edge) {
			e = e.nextInAEL
		}
		edge.nextInAEL = e.nextInAEL
		if e.nextInAEL != nil {
			e.nextInAEL.prevInAEL = edge
		}
		edge.prevInAEL = e
		e.nextInAEL = edge
	}
}

func (self *Clipper) insertLocalMinimaIntoAEL(botY int) {
	for self.currentLocMin != nil &&
		self.currentLocMin.Y == botY {
		lb := self.currentLocMin.leftBound
		rb := self.currentLocMin.rightBound
		self.insertEdgeIntoAEL(lb)
		self.insertScanbeam(lb.Top.Y)
		self.insertEdgeIntoAEL(rb)
		if self.isEvenOddFillType(lb) {
			lb.windDelta = 1
			rb.windDelta = 1
		} else {
			rb.windDelta = -lb.windDelta
		}
		self.setWindingCount(lb)
		rb.windCnt = lb.windCnt
		rb.windCnt2 = lb.windCnt2
		if rb.dx == horizontal {
			self.addEdgeToSEL(rb)
			self.insertScanbeam(rb.nextInLML.Top.Y)
		} else {
			self.insertScanbeam(rb.Top.Y)
		}
		if self.isContributing(lb) {
			self.addLocalMinPoly(lb, rb, &Point{lb.Curr.X, self.currentLocMin.Y})
		}

		if rb.outIdx >= 0 && rb.dx == horizontal && self.horzJoins != nil {
			hj := self.horzJoins
			for {
				_, _, overlap := getOverlapSegment(hj.edge.Bot, hj.edge.Top, rb.Bot, rb.Top)
				if overlap {
					self.addJoin(hj.edge, rb, hj.savedIdx, -1)
				}
				hj = hj.nextHj
				if hj == self.horzJoins {
					break
				}
			}
		}
		if lb.nextInAEL != rb {
			if rb.outIdx >= 0 && rb.prevInAEL.outIdx >= 0 && slopesEqual2(rb.prevInAEL, rb) {
				self.addJoin(rb, rb.prevInAEL, -1, -1)
			}
			e := lb.nextInAEL
			pt := lb.Curr
			for e != rb {
				self.intersectEdges(rb, e, pt, protectsNeither)
				e = e.nextInAEL
			}
		}
		self.popLocalMinima()
	}
}

func (self *Clipper) swapPositionsInAEL(e1, e2 *edge) {
	var nextE, prevE *edge
	if e1.nextInAEL == e2 {
		nextE = e2.nextInAEL
		if nextE != nil {
			nextE.prevInAEL = e1
		}
		prevE = e1.prevInAEL
		if prevE != nil {
			prevE.nextInAEL = e2
		}
		e2.prevInAEL = prevE
		e2.nextInAEL = e1
		e1.prevInAEL = e2
		e1.nextInAEL = nextE
	} else if e2.nextInAEL == e1 {
		nextE = e1.nextInAEL
		if nextE != nil {
			nextE.prevInAEL = e2
		}
		prevE = e2.prevInAEL
		if prevE != nil {
			prevE.nextInAEL = e1
		}
		e1.prevInAEL = prevE
		e1.nextInAEL = e2
		e2.prevInAEL = e1
		e2.nextInAEL = nextE
	} else {
		nextE = e1.nextInAEL
		prevE = e1.prevInAEL
		e1.nextInAEL = e2.nextInAEL
		if e1.nextInAEL != nil {
			e1.nextInAEL.prevInAEL = e1
		}
		e1.prevInAEL = e2.prevInAEL
		if e1.prevInAEL != nil {
			e1.prevInAEL.nextInAEL = e1
		}
		e2.nextInAEL = nextE
		if e2.nextInAEL != nil {
			e2.nextInAEL.prevInAEL = e2
		}
		e2.prevInAEL = prevE
		if e2.prevInAEL != nil {
			e2.prevInAEL.nextInAEL = e2
		}
	}
	if e1.prevInAEL == nil {
		self.activeEdges = e1
	} else if e2.prevInAEL == nil {
		self.activeEdges = e2
	}
}

func (self *Clipper) swapPositionsInSEL(e1, e2 *edge) {
	var nextE, prevE *edge
	if e1.nextInSEL == e2 {
		nextE = e2.nextInSEL
		if nextE != nil {
			nextE.prevInSEL = e1
		}
		prevE = e1.prevInSEL
		if prevE != nil {
			prevE.nextInSEL = e2
		}
		e2.prevInSEL = prevE
		e2.nextInSEL = e1
		e1.prevInSEL = e2
		e1.nextInSEL = nextE
	} else if e2.nextInSEL == e1 {
		nextE = e1.nextInSEL
		if nextE != nil {
			nextE.prevInSEL = e2
		}
		prevE = e2.prevInSEL
		if prevE != nil {
			prevE.nextInSEL = e1
		}
		e1.prevInSEL = prevE
		e1.nextInSEL = e2
		e2.prevInSEL = e1
		e2.nextInSEL = nextE
	} else {
		nextE = e1.nextInSEL
		prevE = e1.prevInSEL
		e1.nextInSEL = e2.nextInSEL
		e1.nextInSEL = e2.nextInSEL
		if e1.nextInSEL != nil {
			e1.nextInSEL.prevInSEL = e1
		}
		e1.prevInSEL = e2.prevInSEL
		if e1.prevInSEL != nil {
			e1.prevInSEL.nextInSEL = e1
		}
		e2.nextInSEL = nextE
		if e2.nextInSEL != nil {
			e2.nextInSEL.prevInSEL = e2
		}
		e2.prevInSEL = prevE
		if e2.prevInSEL != nil {
			e2.prevInSEL.nextInSEL = e2
		}
		if e1.prevInSEL == nil {
			self.sortedEdges = e1
		} else if e2.prevInSEL == nil {
			self.sortedEdges = e2
		}
	}
}

func (self *Clipper) isTopHorz(xPos int) bool {
	e := self.sortedEdges
	for e != nil {
		if (xPos >= min(e.Curr.X, e.Top.X)) && (xPos <= max(e.Curr.X, e.Top.X)) {
			return false
		}
		e = e.nextInSEL
	}
	return true
}

func (self *Clipper) processHorizontal(horzEdge *edge) {
	var horzLeft, horzRight int
	var direction direction
	if horzEdge.Curr.X < horzEdge.Top.X {
		horzLeft = horzEdge.Curr.X
		horzRight = horzEdge.Top.X
		direction = leftToRight
	} else {
		horzLeft = horzEdge.Top.X
		horzRight = horzEdge.Curr.X
		direction = rightToLeft
	}
	var eMaxPair *edge
	if horzEdge.nextInLML == nil {
		eMaxPair = getMaximaPair(horzEdge)
	}
	e := getnextInAEL(horzEdge, direction)
	for e != nil {
		if (e.Curr.X == horzEdge.Top.X) && eMaxPair == nil {
			if slopesEqual2(e, horzEdge.nextInLML) {
				if horzEdge.outIdx >= 0 && e.outIdx >= 0 {
					self.addJoin(horzEdge.nextInLML, e, horzEdge.outIdx, -1)
				}
				break
			} else if e.dx < horzEdge.nextInLML.dx {
				break
			}
		}
		eNext := getnextInAEL(e, direction)
		if eMaxPair != nil ||
			((direction == leftToRight) && (e.Curr.X < horzRight)) ||
			((direction == rightToLeft) && (e.Curr.X > horzLeft)) {
			if e == eMaxPair {
				if direction == leftToRight {
					self.intersectEdges(horzEdge, e,
						&Point{e.Curr.X, horzEdge.Curr.Y}, protectsNeither)
				} else {
					self.intersectEdges(e, horzEdge,
						&Point{e.Curr.X, horzEdge.Curr.Y}, protectsNeither)
				}
				return
			} else if e.dx == horizontal && !isMinima(e) && e.Curr.X <= e.Top.X {
				if direction == leftToRight {
					self.intersectEdges(horzEdge, e, &Point{e.Curr.X, horzEdge.Curr.Y},
						protectRight(!self.isTopHorz(e.Curr.X)))
				} else {
					self.intersectEdges(e, horzEdge, &Point{e.Curr.X, horzEdge.Curr.Y},
						protectLeft(!self.isTopHorz(e.Curr.X)))
				}
			} else if direction == leftToRight {
				self.intersectEdges(horzEdge, e, &Point{e.Curr.X, horzEdge.Curr.Y},
					protectRight(!self.isTopHorz(e.Curr.X)))
			} else {
				self.intersectEdges(e, horzEdge, &Point{e.Curr.X, horzEdge.Curr.Y},
					protectLeft(!self.isTopHorz(e.Curr.X)))
			}
			self.swapPositionsInAEL(horzEdge, e)
		} else if (direction == leftToRight && e.Curr.X >= horzRight) ||
			(direction == rightToLeft && e.Curr.X <= horzLeft) {
			break
		}
		e = eNext
	}
	if horzEdge.nextInLML != nil {
		if horzEdge.outIdx >= 0 {
			self.addOutPt(horzEdge, horzEdge.Top)
		}
		self.updateEdgeIntoAEL(horzEdge)
	} else {
		if horzEdge.outIdx >= 0 {
			self.intersectEdges(horzEdge, eMaxPair,
				&Point{horzEdge.Top.X, horzEdge.Curr.Y}, protectsBoth)
		}
		if eMaxPair.outIdx >= 0 {
			panic("Clipper: Horizontal Error")
			self.deleteFromAEL(eMaxPair)
			self.deleteFromAEL(horzEdge)
		}
	}
}

func (self *Clipper) processHorizontals() {
	for self.sortedEdges != nil {
		e := self.sortedEdges
		self.deleteFromSEL(e)
		self.processHorizontal(e)
	}
}

func (self *Clipper) addJoin(e1, e2 *edge, e1OutIdx, e2OutIdx int) {
	jr := new(joinRec)
	if e1OutIdx >= 0 {
		jr.poly1Idx = e1OutIdx
	} else {
		jr.poly1Idx = e1.outIdx
	}
	jr.pt1a = e1.Curr
	jr.pt1b = e1.Top
	if e2OutIdx >= 0 {
		jr.poly2Idx = e2OutIdx
	} else {
		jr.poly2Idx = e2.outIdx
	}
	jr.pt2a = e2.Curr
	jr.pt2b = e2.Top
	if self.joinList == nil {
		self.joinList = make([]*joinRec, 0)
	}
	self.joinList = append(self.joinList, jr)
}

func (self *Clipper) fixupJoinRecs(jr *joinRec, outPt *outPt, startIdx int) {
	for i := startIdx; i < len(self.joinList); i++ {
		jr2 := self.joinList[i]
		if jr2.poly1Idx == jr.poly1Idx && pointIsVertex(jr2.pt1a, outPt) {
			jr2.poly1Idx = jr.poly2Idx
		}
		if jr2.poly2Idx == jr.poly1Idx && pointIsVertex(jr2.pt2a, outPt) {
			jr2.poly2Idx = jr.poly2Idx
		}
	}
}

func (self *Clipper) addHorzJoin(e *edge, idx int) {
	hj := newHorzJoin(e, idx)
	if self.horzJoins == nil {
		self.horzJoins = hj
		hj.nextHj = hj
		hj.prevHj = hj
	} else {
		hj.nextHj = self.horzJoins
		hj.prevHj = self.horzJoins.prevHj
		self.horzJoins.prevHj.nextHj = hj
		self.horzJoins.prevHj = hj
	}
}

func (self *Clipper) insertIntersectNode(e1, e2 *edge, pt *Point) {
	newNode := &intersectNode{e1: e1, e2: e2, pt: pt}
	if self.intersectNodes == nil {
		self.intersectNodes = newNode
	} else if newNode.pt.Y > self.intersectNodes.pt.Y {
		newNode.nextIn = self.intersectNodes
		self.intersectNodes = newNode
	} else {
		node := self.intersectNodes
		for node.nextIn != nil &&
			newNode.pt.Y < node.nextIn.pt.Y {
			node = node.nextIn
		}
		newNode.nextIn = node.nextIn
		node.nextIn = newNode
	}
}

func (self *Clipper) processIntersections(botY, topY int) bool {
	defer func() {
		self.intersectNodes = nil
		self.sortedEdges = nil
	}()
	self.buildIntersectList(botY, topY)
	if self.intersectNodes == nil {
		return true
	}
	if self.intersectNodes.nextIn != nil &&
		!self.fixupIntersectionOrder() {
		return false
	}
	self.processIntersectList()
	return true
}

func (self *Clipper) buildIntersectList(botY, topY int) {
	e := self.activeEdges
	if e == nil {
		return
	}
	self.sortedEdges = e
	for e != nil {
		e.prevInSEL = e.prevInAEL
		e.nextInSEL = e.nextInAEL
		e.Curr = &Point{topX(e, topY), e.Curr.Y}
		e = e.nextInAEL
	}
	for {
		isModified := false
		e = self.sortedEdges
		for e.nextInSEL != nil {
			eNext := e.nextInSEL
			if e.Curr.X <= eNext.Curr.X {
				e = eNext
				continue
			}
			pt, intersected := intersectPoint(e, eNext)
			if !intersected && e.Curr.X > eNext.Curr.X+1 {
				panic("Intersect Error")
			}
			if pt.Y > botY {
				pt = &Point{topX(e, botY), botY}
			}
			self.insertIntersectNode(e, eNext, pt)
			self.swapPositionsInSEL(e, eNext)
			isModified = true
		}
		if e.prevInSEL != nil {
			e.prevInSEL.nextInSEL = nil
		} else {
			break
		}
		if !isModified {
			break
		}
	}
	self.sortedEdges = nil
	return
}

func (self *Clipper) processIntersectList() {
	for self.intersectNodes != nil {
		node := self.intersectNodes
		self.intersectEdges(node.e1, node.e2, node.pt, protectsBoth)
		self.swapPositionsInAEL(node.e1, node.e2)
		self.intersectNodes = node.nextIn
	}
}

func (self *Clipper) deleteFromAEL(e *edge) {
	aelPrev := e.prevInAEL
	aelNext := e.nextInAEL
	if aelPrev == nil && aelNext == nil && e != self.activeEdges {
		return
	}
	if aelPrev != nil {
		aelPrev.nextInAEL = aelNext
	} else {
		self.activeEdges = aelNext
	}
	if aelNext != nil {
		aelNext.prevInAEL = aelPrev
	}
	e.nextInAEL = nil
	e.prevInAEL = nil
}

func (self *Clipper) deleteFromSEL(e *edge) {
	SELPrev := e.prevInSEL
	SELNext := e.nextInSEL
	if SELPrev == nil && SELNext == nil && e != self.sortedEdges {
		return
	}
	if SELPrev != nil {
		SELPrev.nextInSEL = SELNext
	} else {
		self.sortedEdges = SELNext
	}
	if SELNext != nil {
		SELNext.prevInSEL = SELPrev
	}
	e.nextInSEL = nil
	e.prevInSEL = nil
}

func (self *Clipper) intersectEdges(e1, e2 *edge, pt *Point, protect protects) {
	e1stops := protect&protectsLeft == 0 &&
		e1.nextInLML == nil &&
		e1.Top.X == pt.X && e1.Top.Y == pt.Y
	e2stops := protect&protectsRight == 0 &&
		e2.nextInLML == nil &&
		e2.Top.X == pt.X && e2.Top.Y == pt.Y
	e1Contributing := e1.outIdx >= 0
	e2contributing := e2.outIdx >= 0

	if e1.polyType == e2.polyType {
		if self.isEvenOddFillType(e1) {
			e1Wc := e1.windCnt
			e1.windCnt = e2.windCnt
			e2.windCnt = e1Wc
		} else {
			if e1.windCnt+e2.windDelta == 0 {
				e1.windCnt = -e1.windCnt
			} else {
				e1.windCnt += e2.windDelta
			}
			if e2.windCnt-e1.windDelta == 0 {
				e2.windCnt = -e2.windCnt
			} else {
				e2.windCnt -= e1.windDelta
			}
		}
	} else {
		if !self.isEvenOddFillType(e2) {
			e1.windCnt2 += e2.windDelta
		} else if e1.windCnt2 == 0 {
			e1.windCnt2 = 1
		} else {
			e1.windCnt2 = 0
		}
		if !self.isEvenOddFillType(e1) {
			e2.windCnt2 -= e1.windDelta
		} else if e2.windCnt2 == 0 {
			e2.windCnt2 = 1
		} else {
			e2.windCnt2 = 0
		}
	}
	var e1FillType, e1FillType2, e2FillType, e2FillType2 PolyFillType
	if e1.polyType == Subject {
		e1FillType = self.subjFillType
		e1FillType2 = self.clipFillType
	} else {
		e1FillType = self.clipFillType
		e1FillType2 = self.subjFillType
	}
	if e2.polyType == Subject {
		e2FillType = self.subjFillType
		e2FillType2 = self.clipFillType
	} else {
		e2FillType = self.clipFillType
		e2FillType2 = self.subjFillType
	}
	var e1Wc, e2Wc int
	if e1FillType == Positive {
		e1Wc = e1.windCnt
	} else if e1FillType == Negative {
		e1Wc = -e1.windCnt
	} else {
		e1Wc = intAbs(e1.windCnt)
	}

	if e2FillType == Positive {
		e2Wc = e2.windCnt
	} else if e2FillType == Negative {
		e2Wc = -e2.windCnt
	} else {
		e2Wc = intAbs(e2.windCnt)
	}

	if e1Contributing && e2contributing {
		if e1stops || e2stops ||
			(e1Wc != 0 && e1Wc != 1) || (e2Wc != 0 && e2Wc != 1) ||
			(e1.polyType != e2.polyType && self.clipType != Xor) {
			self.addLocalMaxPoly(e1, e2, pt)
		} else {
			self.addOutPt(e1, pt)
			self.addOutPt(e2, pt)
			swapSides(e1, e2)
			swapPolyIndexes(e1, e2)
		}
	} else if e1Contributing {
		if e2Wc == 0 || e2Wc == 1 {
			self.addOutPt(e1, pt)
			swapSides(e1, e2)
			swapPolyIndexes(e1, e2)
		}
	} else if e2contributing {
		if e1Wc == 0 || e1Wc == 1 {
			self.addOutPt(e2, pt)
			swapSides(e1, e2)
			swapPolyIndexes(e1, e2)
		}
	} else if (e1Wc == 0 || e1Wc == 1) && (e2Wc == 0 || e2Wc == 1) &&
		!e1stops && !e2stops {

		e1FillType2, e2FillType2 = EvenOdd, EvenOdd
		var e1Wc2, e2Wc2 int
		if e1FillType2 == Positive {
			e1Wc2 = e1.windCnt2
		} else if e1FillType2 == Negative {
			e1Wc2 = -e1.windCnt2
		} else {
			e1Wc2 = intAbs(e1.windCnt2)
		}
		if e2FillType2 == Positive {
			e2Wc2 = e2.windCnt2
		} else if e2FillType2 == Negative {
			e2Wc2 = -e2.windCnt2
		} else {
			e2Wc2 = intAbs(e2.windCnt2)
		}

		if e1.polyType != e2.polyType {
			self.addLocalMinPoly(e1, e2, pt)
		} else if e1Wc == 1 && e2Wc == 1 {
			if self.clipType == Intersection {
				if e1Wc2 > 0 && e2Wc2 > 0 {
					self.addLocalMinPoly(e1, e2, pt)
				}
			} else if self.clipType == Union {
				if e1Wc2 <= 0 && e2Wc2 <= 0 {
					self.addLocalMinPoly(e1, e2, pt)
				}
			} else if self.clipType == Difference {
				if (e1.polyType == Clip && e1Wc2 > 0 && e2Wc2 > 0) ||
					(e1.polyType == Subject && e1Wc2 <= 0 && e2Wc2 <= 0) {
					self.addLocalMinPoly(e1, e2, pt)
				}
			} else {
				self.addLocalMinPoly(e1, e2, pt)
			}
		} else {
			swapSides(e1, e2) //, self.polyOutList)
		}
	}
	if e1stops != e2stops &&
		((e1stops && e1.outIdx >= 0) || (e2stops && e2.outIdx >= 0)) {
		swapSides(e1, e2) //, self.polyOutList)
		swapPolyIndexes(e1, e2)
	}
	if e1stops {
		self.deleteFromAEL(e1)
	}
	if e2stops {
		self.deleteFromAEL(e2)
	}
}

func (self *Clipper) doMaxima(e *edge, topY int) {
	eMaxPair := getMaximaPair(e)
	x := e.Top.X
	eNext := e.nextInAEL
	for eNext != eMaxPair {
		if eNext == nil {
			panic("DoMaxima error")
		}
		self.intersectEdges(e, eNext, &Point{x, topY}, protectsBoth)
		self.swapPositionsInAEL(e, eNext)
		eNext = e.nextInAEL
	}
	if e.outIdx < 0 && eMaxPair.outIdx < 0 {
		self.deleteFromAEL(e)
		self.deleteFromAEL(eMaxPair)
	} else if e.outIdx >= 0 && eMaxPair.outIdx >= 0 {
		self.intersectEdges(e, eMaxPair, &Point{x, topY}, protectsNeither)
	} else {
		panic("DoMaxima error")
	}
}

func (self *Clipper) updateEdgeIntoAEL(e *edge) *edge {
	if e.nextInLML == nil {
		panic("UpdateEdgeIntoAEL error")
	}
	aelPrev := e.prevInAEL
	aelNext := e.nextInAEL
	e.nextInLML.outIdx = e.outIdx
	if aelPrev != nil {
		aelPrev.nextInAEL = e.nextInLML
	} else {
		self.activeEdges = e.nextInLML
	}
	if aelNext != nil {
		aelNext.prevInAEL = e.nextInLML
	}
	e.nextInLML.side = e.side
	e.nextInLML.windDelta = e.windDelta
	e.nextInLML.windCnt = e.windCnt
	e.nextInLML.windCnt2 = e.windCnt2
	e = e.nextInLML
	e.prevInAEL = aelPrev
	e.nextInAEL = aelNext
	if e.dx != horizontal {
		self.insertScanbeam(e.Top.Y)
	}
	return e
}

func (self *Clipper) addLocalMinPoly(e1, e2 *edge, pt *Point) {
	var e, prevE *edge
	if e2.dx == horizontal || e1.dx > e2.dx {
		self.addOutPt(e1, pt)
		e2.outIdx = e1.outIdx
		e1.side = leftEdge
		e2.side = rightEdge
		e = e1
		if e.prevInAEL == e2 {
			prevE = e2.prevInAEL
		} else {
			prevE = e1.prevInAEL
		}
	} else {
		self.addOutPt(e2, pt)
		e1.outIdx = e2.outIdx
		e1.side = rightEdge
		e2.side = leftEdge
		e = e2
		if e.prevInAEL == e1 {
			prevE = e1.prevInAEL
		} else {
			prevE = e.prevInAEL
		}
	}
	if prevE != nil && prevE.outIdx >= 0 &&
		topX(prevE, pt.Y) == topX(e, pt.Y) &&
		slopesEqual2(e, prevE) {
		self.addJoin(e, prevE, -1, -1)
	}
	return
}

func (self *Clipper) addLocalMaxPoly(e1, e2 *edge, pt *Point) {
	self.addOutPt(e1, pt)
	if e1.outIdx == e2.outIdx {
		e1.outIdx = -1
		e2.outIdx = -1
	} else if e1.outIdx < e2.outIdx {
		self.appendPolygon(e1, e2)
	} else {
		self.appendPolygon(e2, e1)
	}
}

func (self *Clipper) createOutRec() *outRec {
	outrec := newOutRec(len(self.polyOutList))
	self.polyOutList = append(self.polyOutList, outrec)
	return outrec
}

func (self *Clipper) addOutPt(e *edge, pt *Point) {
	toFront := e.side == leftEdge
	var outrec *outRec
	if e.outIdx < 0 {
		outrec = self.createOutRec()
		e.outIdx = outrec.idx
		op := newOutPt(outrec.idx, pt)
		op.nextOp = op
		op.prevOp = op
		outrec.pts = op
		setHoleState(e, outrec, self.polyOutList)
	} else {
		outrec = self.polyOutList[e.outIdx]
		op := outrec.pts
		if (toFront && pointsEqual(pt, op.pt)) ||
			(!toFront && pointsEqual(pt, op.prevOp.pt)) {
			return
		}
		op2 := newOutPt(outrec.idx, pt)
		op2.nextOp = op
		op2.prevOp = op.prevOp
		op.prevOp.nextOp = op2
		op.prevOp = op2
		if toFront {
			outrec.pts = op2
		}
	}
}

func (self *Clipper) appendPolygon(e1, e2 *edge) {
	outRec1 := self.polyOutList[e1.outIdx]
	outRec2 := self.polyOutList[e2.outIdx]
	var holeStateRec *outRec
	if param1RightOfParam2(outRec1, outRec2) {
		holeStateRec = outRec2
	} else if param1RightOfParam2(outRec2, outRec1) {
		holeStateRec = outRec1
	} else {
		holeStateRec = getLowermostRec(outRec1, outRec2)
	}

	p1_lft := outRec1.pts
	p2_lft := outRec2.pts
	p1_rt := p1_lft.prevOp
	p2_rt := p2_lft.prevOp
	newSide := leftEdge

	if e1.side == leftEdge {
		if e2.side == leftEdge {
			// z y x a b c
			reversePolyPtLinks(p2_lft)
			p2_lft.nextOp = p1_lft
			p1_lft.prevOp = p2_lft
			p1_rt.nextOp = p2_rt
			p2_rt.prevOp = p1_rt
			outRec1.pts = p2_rt
		} else {
			// x y z a b c
			p2_rt.nextOp = p1_lft
			p1_lft.prevOp = p2_rt
			p2_lft.prevOp = p1_rt
			p1_rt.nextOp = p2_lft
			outRec1.pts = p2_lft
		}
	} else {
		newSide = rightEdge
		if e2.side == rightEdge {
			// a b c z y x
			reversePolyPtLinks(p2_lft)
			p1_rt.nextOp = p2_rt
			p2_rt.prevOp = p1_rt
			p2_lft.nextOp = p1_lft
			p1_lft.prevOp = p2_lft
		} else {
			// a b c x y z
			p1_rt.nextOp = p2_lft
			p2_lft.prevOp = p1_rt
			p1_lft.prevOp = p2_rt
			p2_rt.nextOp = p1_lft
		}
	}
	outRec1.bottomPt = nil
	if holeStateRec == outRec2 {
		if outRec2.FirstLeft != outRec1 {
			outRec1.FirstLeft = outRec2.FirstLeft
		}
		outRec1.isHole = outRec2.isHole
	}
	outRec2.pts = nil
	outRec2.bottomPt = nil
	outRec2.FirstLeft = outRec1
	OKIdx := outRec1.idx
	ObsoleteIdx := outRec2.idx

	e1.outIdx = -1
	e2.outIdx = -1

	e := self.activeEdges
	for e != nil {
		if e.outIdx == ObsoleteIdx {
			e.outIdx = OKIdx
			e.side = newSide
			break
		}
		e = e.nextInAEL
	}
	outRec2.idx = outRec1.idx
}

func (self *Clipper) fixupIntersectionOrder() bool {
	self.copyAELToSEL()
	inode := self.intersectNodes
	for inode != nil {
		if !edgesAdjacent(inode) {
			nextNode := inode.nextIn
			for nextNode != nil && !edgesAdjacent(nextNode) {
				nextNode = nextNode.nextIn
			}
			if nextNode == nil {
				return false
			}
			e1 := inode.e1
			e2 := inode.e2
			p := inode.pt
			inode.e1 = nextNode.e1
			inode.e2 = nextNode.e2
			inode.pt = nextNode.pt
			nextNode.e1 = e1
			nextNode.e2 = e2
			nextNode.pt = p
		}
		self.swapPositionsInSEL(inode.e1, inode.e2)
		inode = inode.nextIn
	}
	return true
}

func (self *Clipper) processEdgesAtTopOfScanbeam(topY int) {
	e := self.activeEdges
	var ePrev, eNext *edge
	for e != nil {
		if isMaxima(e, topY) && getMaximaPair(e).dx != horizontal {
			ePrev = e.prevInAEL
			self.doMaxima(e, topY)
			if ePrev == nil {
				e = self.activeEdges
			} else {
				e = ePrev.nextInAEL
			}
		} else {
			intermediateVert := isIntermediate(e, topY)
			if intermediateVert && e.nextInLML.dx == horizontal {
				if e.outIdx >= 0 {
					self.addOutPt(e, e.Top)
					hj := self.horzJoins
					if hj != nil {
						for {
							_, _, overlap := getOverlapSegment(
								hj.edge.Bot, hj.edge.Top, e.nextInLML.Bot, e.nextInLML.Top)
							if overlap {
								self.addJoin(hj.edge, e.nextInLML, hj.savedIdx, e.outIdx)
							}
							hj = hj.nextHj
						}
						if hj == self.horzJoins {
							break
						}
						self.addHorzJoin(e.nextInLML, e.outIdx)
					}
				}
				e = self.updateEdgeIntoAEL(e)
				self.addEdgeToSEL(e)
			} else {
				e.Curr = &Point{topX(e, topY), topY}
				if self.ForceSimple && e.prevInAEL != nil &&
					e.prevInAEL.Curr.X == e.Curr.X &&
					e.outIdx >= 0 && e.prevInAEL.outIdx >= 0 {
					if intermediateVert {
						self.addOutPt(e.prevInAEL, &Point{e.Curr.X, topY})
					} else {
						self.addOutPt(e, &Point{e.Curr.X, topY})
					}
				}
			}
			e = e.nextInAEL
		}
	}
	self.processHorizontals()

	e = self.activeEdges
	for e != nil {
		if isIntermediate(e, topY) {
			if e.outIdx >= 0 {
				self.addOutPt(e, e.Top)
			}
			e = self.updateEdgeIntoAEL(e)

			ePrev = e.prevInAEL
			eNext = e.nextInAEL
			if ePrev != nil && ePrev.Curr.X == e.Bot.X &&
				(ePrev.Curr.Y == e.Bot.Y) && (e.outIdx >= 0) &&
				(ePrev.outIdx >= 0) && (ePrev.Curr.Y > ePrev.Top.Y) &&
				slopesEqual2(e, ePrev) {
				self.addOutPt(ePrev, e.Bot)
				self.addJoin(e, ePrev, -1, -1)
			} else if eNext != nil && (eNext.Curr.X == e.Bot.X) &&
				(eNext.Curr.Y == e.Bot.Y) && (e.outIdx >= 0) &&
				(eNext.outIdx >= 0) && (eNext.Curr.Y > eNext.Top.Y) &&
				slopesEqual2(e, eNext) {
				self.addOutPt(eNext, e.Bot)
				self.addJoin(e, eNext, -1, -1)
			}
		}
		e = e.nextInAEL
	}
}

// see http://www.mathopenref.com/coordpolygonarea2.html
func (self *Clipper) area(pts *outPt) float64 {
	result := 0.0
	p := pts
	for {
		result += float64((p.pt.X + p.prevOp.pt.X) * (p.prevOp.pt.Y - p.pt.Y))
		p = p.nextOp
		if p == pts {
			break
		}
	}
	return result / 2.
}

func (self *Clipper) joinPoints(jr *joinRec) (*outPt, *outPt, bool) {
	var p1, p2, p3, p4 *outPt
	var result bool
	outRec1 := self.polyOutList[jr.poly1Idx]
	outRec2 := self.polyOutList[jr.poly2Idx]
	if outRec1 == nil || outRec2 == nil {
		return p1, p2, false
	}
	pp1a := outRec1.pts
	pp2a := outRec2.pts
	pt1 := jr.pt2a
	pt2 := jr.pt2b
	pt3 := jr.pt1a
	pt4 := jr.pt1b
	pp1a, pt1, pt2, result = findSegment(pp1a, pt1, pt2)
	if !result {
		return p1, p2, false
	}
	if outRec1 == outRec2 {
		pp2a = pp1a.nextOp
		pp2a, pt3, pt4, result = findSegment(pp2a, pt3, pt4)
		if !result || pp2a == pp1a {
			return p1, p2, false
		}
	} else {
		pp2a, pt3, pt4, result = findSegment(pp2a, pt3, pt4)
		if !result {
			return p1, p2, false
		}
	}
	pt1, pt2, result = getOverlapSegment(pt1, pt2, pt3, pt4)
	if !result {
		return p1, p2, false
	}

	prevOp := pp1a.prevOp
	if pointsEqual(pp1a.pt, pt1) {
		p1 = pp1a
	} else if pointsEqual(prevOp.pt, pt1) {
		p1 = prevOp
	} else {
		p1 = insertPolyPtBetween(pp1a, prevOp, pt1)
	}

	if pointsEqual(pp1a.pt, pt2) {
		p2 = pp1a
	} else if pointsEqual(prevOp.pt, pt2) {
		p2 = prevOp
	} else if (p1 == pp1a) || (p1 == prevOp) {
		p2 = insertPolyPtBetween(pp1a, prevOp, pt2)
	} else if pt3IsBetweenPt1AndPt2(pp1a.pt, p1.pt, pt2) {
		p2 = insertPolyPtBetween(pp1a, p1, pt2)
	} else {
		p2 = insertPolyPtBetween(p1, prevOp, pt2)
	}

	prevOp = pp2a.prevOp
	if pointsEqual(pp2a.pt, pt1) {
		p3 = pp2a
	} else if pointsEqual(prevOp.pt, pt1) {
		p3 = prevOp
	} else {
		p3 = insertPolyPtBetween(pp2a, prevOp, pt1)
	}
	if pointsEqual(pp2a.pt, pt2) {
		p4 = pp2a
	} else if pointsEqual(prevOp.pt, pt2) {
		p4 = prevOp
	} else if (p3 == pp2a) || (p3 == prevOp) {
		p4 = insertPolyPtBetween(pp2a, prevOp, pt2)
	} else if pt3IsBetweenPt1AndPt2(pp2a.pt, p3.pt, pt2) {
		p4 = insertPolyPtBetween(pp2a, p3, pt2)
	} else {
		p4 = insertPolyPtBetween(p3, prevOp, pt2)
	}

	if p1.nextOp == p2 && p3.prevOp == p4 {
		p1.nextOp = p3
		p3.prevOp = p1
		p2.prevOp = p4
		p4.nextOp = p2
		return p1, p2, true
	} else if p1.prevOp == p2 && p3.nextOp == p4 {
		p1.prevOp = p3
		p3.nextOp = p1
		p2.nextOp = p4
		p4.prevOp = p2
		return p1, p2, true
	}
	return p1, p2, false
}

func (self *Clipper) fixupFirstLefts1(oldOutRec, newOutRec *outRec) {
	for _, outRec := range self.polyOutList {
		if outRec.pts != nil && outRec.FirstLeft == oldOutRec {
			if poly2ContainsPoly1(outRec.pts, newOutRec.pts) {
				outRec.FirstLeft = newOutRec
			}
		}
	}
}

func (self *Clipper) fixupFirstLefts2(oldOutRec, newOutRec *outRec) {
	for _, outRec := range self.polyOutList {
		if outRec.FirstLeft == oldOutRec {
			outRec.FirstLeft = newOutRec
		}
	}
}

func (self *Clipper) getOutRec(idx int) *outRec {
	outrec := self.polyOutList[idx]
	for outrec != self.polyOutList[outrec.idx] {
		outrec = self.polyOutList[outrec.idx]
	}
	return outrec
}

func (self *Clipper) joinCommonEdges() {
	for i := 0; i < len(self.joinList); i++ {
		jr := self.joinList[i]
		outRec1 := self.getOutRec(jr.poly1Idx)
		outRec2 := self.getOutRec(jr.poly2Idx)
		if outRec1.pts == nil || outRec2.pts == nil {
			continue
		}

		var holeStateRec *outRec
		if outRec1 == outRec2 {
			holeStateRec = outRec1
		} else if param1RightOfParam2(outRec1, outRec2) {
			holeStateRec = outRec2
		} else if param1RightOfParam2(outRec2, outRec1) {
			holeStateRec = outRec1
		} else {
			holeStateRec = getLowermostRec(outRec1, outRec2)
		}

		p1, p2, result := self.joinPoints(jr)
		if !result {
			continue
		}

		if outRec1 == outRec2 {
			outRec1.pts = p1
			outRec1.bottomPt = nil
			outRec2 = self.createOutRec()
			outRec2.pts = p2
			jr.poly2Idx = outRec2.idx

			if poly2ContainsPoly1(outRec2.pts, outRec1.pts) {
				outRec2.isHole = !outRec1.isHole
				outRec2.FirstLeft = outRec1

				self.fixupJoinRecs(jr, p2, i+1)

				if self.usingPolyTree {
					self.fixupFirstLefts2(outRec2, outRec1)
				}

				fixupOutPolygon(outRec1)
				fixupOutPolygon(outRec2)

				if (outRec2.isHole != self.ReverseSolution) ==
					(self.area(outRec2.pts) > 0.0) {
					reversePolyPtLinks(outRec2.pts)
				}

			} else if poly2ContainsPoly1(outRec1.pts, outRec2.pts) {
				outRec2.isHole = outRec1.isHole
				outRec1.isHole = !outRec2.isHole
				outRec2.FirstLeft = outRec1.FirstLeft
				outRec1.FirstLeft = outRec2

				self.fixupJoinRecs(jr, p2, i+1)

				if self.usingPolyTree {
					self.fixupFirstLefts2(outRec1, outRec2)
				}

				fixupOutPolygon(outRec1)
				fixupOutPolygon(outRec2)

				if (outRec1.isHole != self.ReverseSolution) ==
					(self.area(outRec1.pts) > 0.0) {
					reversePolyPtLinks(outRec1.pts)
				}
			} else {
				outRec2.isHole = outRec1.isHole
				outRec2.FirstLeft = outRec1.FirstLeft

				self.fixupJoinRecs(jr, p2, i+1)
				if self.usingPolyTree {
					self.fixupFirstLefts1(outRec1, outRec2)
				}

				fixupOutPolygon(outRec1)
				fixupOutPolygon(outRec2)
			}
		} else {
			fixupOutPolygon(outRec1)
			outRec2.pts = nil
			outRec2.bottomPt = nil
			outRec2.idx = outRec1.idx

			outRec1.isHole = holeStateRec.isHole
			if holeStateRec == outRec2 {
				outRec1.FirstLeft = outRec2.FirstLeft
			}
			outRec2.FirstLeft = outRec1

			if self.usingPolyTree {
				self.fixupFirstLefts2(outRec2, outRec1)
			}
		}
	}
	return
}

func (self *Clipper) doSimplePolygons() {
	i := 0
	for i < len(self.polyOutList) {
		outrec := self.polyOutList[i]
		i += 1
		op := outrec.pts
		if op == nil {
			continue
		}
		for {
			op2 := op.nextOp
			for op2 != outrec.pts {
				if pointsEqual(op.pt, op2.pt) && op2.nextOp != op && op2.prevOp != op {
					//split the polygon into two ...
					op3 := op.prevOp
					op4 := op2.prevOp
					op.prevOp = op4
					op4.nextOp = op
					op2.prevOp = op3
					op3.nextOp = op2

					outrec.pts = op
					outrec2 := self.createOutRec()
					outrec2.pts = op2
					updateOutPtIdxs(outrec2)
					if poly2ContainsPoly1(outrec2.pts, outrec.pts) {
						//OutRec2 is contained by OutRec1 ...
						outrec2.isHole = !outrec.isHole
						outrec2.FirstLeft = outrec
					} else if poly2ContainsPoly1(outrec.pts, outrec2.pts) {
						//OutRec1 is contained by OutRec2 ...
						outrec2.isHole = outrec.isHole
						outrec.isHole = !outrec2.isHole
						outrec2.FirstLeft = outrec.FirstLeft
						outrec.FirstLeft = outrec2
					} else {
						//the 2 polygons are separate ...
						outrec2.isHole = outrec.isHole
						outrec2.FirstLeft = outrec.FirstLeft
					}
					op2 = op // ie get ready for the next iteration
				}
				op2 = op2.nextOp
			}
			op = op.nextOp
			if op == outrec.pts {
				break
			}
		}
	}
	return
}

func (self *Clipper) executeInternal() bool {
	defer func() {
		self.joinList = nil
		self.horzJoins = nil
	}()
	self.reset()
	if self.scanbeam == nil {
		return true
	}
	botY := self.popScanbeam()
	for {
		self.insertLocalMinimaIntoAEL(botY)
		self.horzJoins = nil
		self.processHorizontals()
		topY := self.popScanbeam()
		if !self.processIntersections(botY, topY) {
			return false
		}
		self.processEdgesAtTopOfScanbeam(topY)
		botY = topY
		if self.scanbeam == nil && self.currentLocMin == nil {
			break
		}
	}

	for _, outRec := range self.polyOutList {
		if outRec.pts == nil {
			continue
		}
		fixupOutPolygon(outRec)
		if outRec.pts == nil {
			continue
		}
		if (outRec.isHole != self.ReverseSolution) == (self.area(outRec.pts) > 0.0) {
			reversePolyPtLinks(outRec.pts)
		}
	}

	if self.joinList != nil {
		self.joinCommonEdges()
	}
	if self.ForceSimple {
		self.doSimplePolygons()
	}

	return true
}

func (self *Clipper) Execute(
	clipType ClipType,
	solution [][]*Point,
	subjFillType,
	clipFillType PolyFillType) bool {
	if self.executeLocked {
		return false
	}
	defer func() {
		self.executeLocked = false
		self.usingPolyTree = false
	}()

	self.executeLocked = true
	self.usingPolyTree = true
	solution = nil
	self.subjFillType = subjFillType
	self.clipFillType = clipFillType
	self.clipType = clipType
	result := self.executeInternal()
	if result {
		self.buildResult(solution)
	}

	return result
}

func (self *Clipper) Execute2(
	clipType ClipType,
	solutionTree *PolyTree,
	subjFillType,
	clipFillType PolyFillType) bool {
	if self.executeLocked {
		return false
	}
	defer func() {
		self.executeLocked = false
		self.usingPolyTree = false
	}()
	self.executeLocked = true
	self.usingPolyTree = true
	solutionTree.Clear()
	self.subjFillType = subjFillType
	self.clipFillType = clipFillType
	self.clipType = clipType
	result := self.executeInternal()
	if result {
		self.buildResult2(solutionTree)
	}

	return result
}

func (self *Clipper) buildResult(polygons [][]*Point) {
	for _, outRec := range self.polyOutList {
		if outRec == nil {
			continue
		}
		cnt := pointCount(outRec.pts)
		if cnt < 3 {
			continue
		}
		poly := make([]*Point, 0)
		op := outRec.pts
		for i := 0; i < cnt; i++ {
			poly = append(poly, op.pt)
			op = op.prevOp
		}
		polygons = append(polygons, poly)
	}
	return
}

func (self *Clipper) buildResult2(polyTree *PolyTree) {
	for _, outRec := range self.polyOutList {
		if outRec == nil {
			continue
		}
		cnt := pointCount(outRec.pts)
		if cnt < 3 {
			continue
		}
		fixHoleLinkage(outRec)

		// add nodes to allNodes list ...
		polyNode := new(PolyNode)
		polyTree.allNodes = append(polyTree.allNodes, polyNode)
		outRec.polyNode = polyNode
		op := outRec.pts
		for {
			polyNode.Contour = append(polyNode.Contour, op.pt)
			op = op.prevOp
			if op == outRec.pts {
				break
			}
		}
	}
	// build the tree ...
	for _, outRec := range self.polyOutList {
		if outRec.polyNode == nil {
			continue
		}
		if outRec.FirstLeft == nil {
			polyTree.addChild(outRec.polyNode)
		} else {
			outRec.FirstLeft.polyNode.addChild(outRec.polyNode)
		}
	}
	return
}

//===============================================================================
// OffsetPolygons (+ ancilliary functions)
//===============================================================================

func getUnitNormal(pt1, pt2 *Point) *FloatPoint {
	if pt2.X == pt1.X && pt2.Y == pt1.Y {
		return &FloatPoint{0.0, 0.0}
	}
	dx := float64(pt2.X - pt1.X)
	dy := float64(pt2.Y - pt1.Y)
	f := 1.0 / math.Hypot(dx, dy)
	dx = float64(dx) * f
	dy = float64(dy) * f
	return &FloatPoint{dy, -dx}
}

func getBounds(pts [][]*Point) *rect {
	var left, top, right, bottom int
	for _, poly := range pts {
		for _, pt := range poly {
			left = pt.X
			top = pt.Y
			right = pt.X
			bottom = pt.Y
			break
		}
		break
	}

	for _, poly := range pts {
		for _, pt := range poly {
			if pt.X < left {
				left = pt.X
			}
			if pt.X > right {
				right = pt.X
			}
			if pt.Y < top {
				top = pt.Y
			}
			if pt.Y > bottom {
				bottom = pt.Y
			}
		}
	}
	return &rect{left, top, right, bottom}
}

func getLowestPt(poly []*Point) *Point {
	// precondition: poly must not be empty
	result := poly[0]
	for _, pt := range poly {
		if pt.Y > result.Y || (pt.Y == result.Y && pt.X < result.X) {
			result = pt
		}
	}
	return result
}

func stripDupPts(poly []*Point) []*Point {
	if len(poly) == 0 {
		return poly
	}
	for i := 1; i < len(poly); i++ {
		if pointsEqual(poly[i-1], poly[i]) {
			poly = append(poly[:i], poly[i+1:]...) // remove item i
		}
	}
	i := len(poly) - 1
	for i > 0 && pointsEqual(poly[i], poly[0]) {
		poly = append(poly[:i], poly[i+1:]...) // remove item i
		i -= 1
	}
	return poly
}

func offsetInternal(polys [][]*Point, isPolygon bool, delta float64,
	jointype JoinType, endtype EndType, limit float64) [][]*Point {

	var sinA, mcos, msin, step360, miterLim float64
	var pts []*Point
	var Normals []*FloatPoint
	var k, j int
	var result polySorter
	doSquare := func(pt *Point) {
		// see offset_triginometry.svg in the documentation folder ...
		dx := math.Tan(math.Atan2(sinA,
			(Normals[k].X*Normals[j].X+Normals[k].Y*Normals[j].Y)) / 4)
		result = append(result, &Point{
			round(float64(pt.X) + delta*(Normals[k].X-
				Normals[k].Y*dx)),
			round(float64(pt.Y) + delta*(Normals[k].Y+Normals[k].X*dx))})
		result = append(result, &Point{
			round(float64(pt.X) + delta*(Normals[j].X+Normals[j].Y*dx)),
			round(float64(pt.Y) + delta*(Normals[j].Y-Normals[j].X*dx))})
		return
	}

	doMiter := func(pt *Point, r float64) {
		q := delta / r
		result = append(result, &Point{
			round(float64(pt.X) + (Normals[k].X+Normals[j].X)*q),
			round(float64(pt.Y) + (Normals[k].Y+Normals[j].Y)*q)})
		return
	}

	doRound := func(pt *Point) {
		a := math.Atan2(sinA,
			(Normals[k].X*Normals[j].X + Normals[k].Y*Normals[j].Y))
		steps := round(step360 * math.Abs(a))
		X, Y := Normals[k].X, Normals[k].Y
		for i := 0; i < steps; i++ {
			result = append(result, &Point{
				round(float64(pt.X) + X*delta),
				round(float64(pt.Y) + Y*delta)})
			X2 := X
			X = X*mcos - msin*Y
			Y = X2*msin + Y*mcos
		}
		result = append(result, &Point{round(float64(pt.X) +
			Normals[j].X*delta), round(float64(pt.Y) + Normals[j].Y*delta)})
		return
	}

	GetSin := func() float64 {
		result := (Normals[k].X*Normals[j].Y - Normals[j].X*Normals[k].Y)
		if result > 1.0 {
			result = 1.0
		} else if result < -1.0 {
			result = -1.0
		}
		return result
	}

	offsetPoint := func(jointype JoinType) (j int) {
		if sinA*delta < 0 {
			result = append(result, &Point{round(float64(pts[j].X) +
				Normals[k].X*delta),
				round(float64(pts[j].Y) + Normals[k].Y*delta)})
			result = append(result, pts[j])
			result = append(result, &Point{round(float64(pts[j].X) +
				Normals[j].X*delta),
				round(float64(pts[j].Y) + Normals[j].Y*delta)})
		} else if jointype == MiterJoin {
			r := 1.0 + (Normals[j].X*Normals[k].X + Normals[j].Y*Normals[k].Y)
			if r >= miterLim {
				doMiter(pts[j], r)
			} else {
				doSquare(pts[j])
			}
		} else if jointype == SquareJoin {
			doSquare(pts[j])
		} else {
			doRound(pts[j])
		}
		return j
	}

	if delta == 0 {
		return polys
	}
	if !isPolygon && delta < 0 {
		delta = -delta
	}

	if jointype == MiterJoin {
		// miterLim: see offset_triginometry3.svg in the documentation folder ...
		if limit > 2 {
			miterLim = 2 / (limit * limit)
		} else {
			miterLim = 0.5
		}
		if endtype == RoundEnd {
			limit = 0.25
		}
	}

	if jointype == RoundJoin || endtype == RoundEnd {
		if limit <= 0 {
			limit = 0.25
		} else if limit > math.Abs(delta)*0.25 {
			limit = math.Abs(delta) * 0.25
		}
		// step360: see offset_triginometry2.svg in the documentation folder ...
		step360 := math.Pi / math.Acos(1-limit/math.Abs(delta))
		msin := math.Sin(2 * math.Pi / step360)
		//mcos := math.Cos(2 * math.Pi / step360)
		step360 /= math.Pi * 2
		if delta < 0 {
			msin = -msin
		}
	}

	res := make([][]*Point, 0)
	ppts := make([][]*Point, len(polys))
	copy(ppts, polys)
	for _, pts := range ppts {
		Normals = make([]*FloatPoint, 0)
		result = make([]*Point, 0)
		cnt := len(pts)

		if cnt == 0 || cnt < 3 && delta <= 0 {
			continue
		}

		if cnt == 1 {
			if jointype == RoundJoin {
				X, Y := 1.0, 0.0
				for i := 0; i < round(step360*2*math.Pi); i++ {
					result = append(result, &Point{round(float64(pts[0].X) +
						X*delta), round(float64(pts[0].Y) + Y*delta)})
					X2 := X
					X = X*mcos - msin*Y
					Y = X2*msin + Y*mcos
				}
			} else {
				X, Y := -1.0, -1.0
				for i := 0; i < 4; i++ {
					result = append(result, &Point{round(
						float64(pts[0].X) + X*delta),
						round(float64(pts[0].Y) + Y*delta)})
					if X < 0 {
						X = 1
					} else if Y < 0 {
						Y = 1
					} else {
						X = -1
					}
				}
			}
			continue
		}

		forceClose := pointsEqual(pts[0], pts[cnt-1])
		if forceClose {
			cnt -= 1
		}

		for j := 0; j < cnt-1; j++ {
			Normals = append(Normals, getUnitNormal(pts[j], pts[j+1]))
		}
		if isPolygon || forceClose {
			Normals = append(Normals, getUnitNormal(pts[cnt-1], pts[0]))
		} else {
			Normals = append(Normals, Normals[cnt-2])
		}

		if isPolygon || forceClose {
			k = cnt - 1
			for j := 0; j < cnt; j++ {
				sinA = GetSin()
				k = offsetPoint(jointype)
			}
			res = append(res, result)

			if !isPolygon {
				result = make([]*Point, 0)
				delta = -delta
				k = cnt - 1
				for j := 0; j < cnt; j++ {
					sinA = GetSin()
					k = offsetPoint(jointype)
				}
				delta = -delta
				sort.Reverse(result)
				res = append(res, []*Point(result))
			}

		} else {
			// offset the polyline going forward ...
			var pt1 *Point
			k = 0
			for j := 1; j < cnt-1; j++ {
				sinA = GetSin()
				k = offsetPoint(jointype)
			}

			// handle the end (butt, round || square) ...
			if endtype == ButtEnd {
				j = cnt - 1
				pt1 = &Point{round(float64(pts[j].X) + Normals[j].X*delta),
					round(float64(pts[j].Y) + Normals[j].Y*delta)}
				result = append(result, pt1)
				pt1 = &Point{round(float64(pts[j].X) - Normals[j].X*delta),
					round(float64(pts[j].Y) - Normals[j].Y*delta)}
				result = append(result, pt1)
			} else {
				j = cnt - 1
				k = cnt - 2
				Normals[j] = &FloatPoint{-Normals[j].X, -Normals[j].Y}
				if endtype == SquareEnd {
					doSquare(pts[j])
				} else {
					doRound(pts[j])
				}
			}

			// re-build Normals ...
			for j := cnt - 1; j > 0; j-- {
				Normals[j] = &FloatPoint{-Normals[j-1].X, -Normals[j-1].Y}
			}
			Normals[0] = &FloatPoint{-Normals[1].X, -Normals[1].Y}

			// offset the polyline going backward ...
			k = cnt - 1
			for j := cnt - 2; j > 0; j-- {
				sinA = GetSin()
				k = offsetPoint(jointype)
			}

			// finally handle the start (butt, round || square) ...
			if endtype == ButtEnd {
				pt1 = &Point{round(float64(pts[0].X) - Normals[0].X*delta),
					round(float64(pts[0].Y) - Normals[0].Y*delta)}
				result = append(result, pt1)
				pt1 = &Point{round(float64(pts[0].X) + Normals[0].X*delta),
					round(float64(pts[0].Y) + Normals[0].Y*delta)}
				result = append(result, pt1)
			} else {
				j = 0
				k = 1
				if endtype == SquareEnd {
					doSquare(pts[0])
				} else {
					doRound(pts[0])
				}
			}
			res = append(res, result)
		}
	}
	c := NewClipper()
	c.AddPolygons(res, Subject)
	if delta > 0 {
		c.Execute(Union, res, Positive, Positive)
	} else {
		bounds := getBounds(res)
		outer := make([]*Point, 0)
		outer = append(outer, &Point{bounds.left - 10, bounds.bottom + 10})
		outer = append(outer, &Point{bounds.right + 10, bounds.bottom + 10})
		outer = append(outer, &Point{bounds.right + 10, bounds.top - 10})
		outer = append(outer, &Point{bounds.left - 10, bounds.top - 10})
		c.AddPolygon(outer, Subject)
		c.ReverseSolution = true
		c.Execute(Union, res, Negative, Negative)
		if len(res) > 0 {
			res = res[1:len(res)]
		}
	}
	return res
}

// defaults: jointype = SquareJoin, limit = 0.0, autoFix = true
func OffsetPolygons(polys [][]*Point, delta float64, jointype JoinType,
	limit float64, autoFix bool) [][]*Point {
	if !autoFix {
		return offsetInternal(polys, true, delta, jointype, ButtEnd, limit)
	}
	pts := make([][]*Point, len(polys))
	copy(pts, polys)
	botPoly := make([]*Point, 0)
	var botPt *Point
	for _, poly := range pts {
		poly = stripDupPts(poly)
		if len(poly) < 3 {
			continue
		}
		bot := getLowestPt(poly)
		if botPt == nil || (bot.Y > botPt.Y) ||
			(bot.Y == botPt.Y && bot.X < botPt.X) {
			botPt = bot
			botPoly = poly
		}
	}
	if botPt == nil {
		return nil
	}
	// if the outermost polygon has the wrong orientation,
	// reverse the orientation of all the polygons ...
	if Area(botPoly) < 0.0 {
		for i := 0; i < len(pts); i++ {
			ps := polySorter(pts[i])
			sort.Reverse(ps)
			pts[i] = []*Point(ps)
		}
	}
	return offsetInternal(pts, true, delta, jointype, ButtEnd, limit)
}

//defaults: jointype = SquareJoin, endtype = SquareEnd, limit = 0.0
func OffsetPolyLines(polys [][]*Point, delta float64, jointype JoinType,
	endtype EndType, limit float64) [][]*Point {
	polys2 := make([][]*Point, len(polys))
	copy(polys2, polys)
	for _, p := range polys2 {
		if len(p) == 0 {
			continue
		}
		for i := 1; i < len(p); i++ {
			if pointsEqual(p[i-1], p[i]) {
				p = append(p[:i], p[i+1:]...) // remove item i
			}
		}
	}

	if endtype == ClosedEnd {
		for i := 0; i < len(polys2); i++ {
			ps := polySorter(polys2[i])
			sort.Reverse(ps)
			polys2 = append(polys2, []*Point(ps))
		}
		return offsetInternal(polys2, true, delta, jointype, ButtEnd, limit)
	} else {
		return offsetInternal(polys2, false, delta, jointype, endtype, limit)
	}
}

func distanceSqrd(pt1, pt2 *Point) int {
	dx := (pt1.X - pt2.X)
	dy := (pt1.Y - pt2.Y)
	return (dx*dx + dy*dy)
}

func closestPointOnLine(pt, linePt1, linePt2 *Point) *Point {
	dx := linePt2.X - linePt1.X
	dy := linePt2.Y - linePt1.Y
	if dx == 0 && dy == 0 {
		return &Point{linePt1.X, linePt1.Y}
	}
	q := ((pt.X-linePt1.X)*dx + (pt.Y-linePt1.Y)*dy) / (dx*dx + dy*dy)
	return &Point{
		(1-q)*linePt1.X + q*linePt2.X,
		(1-q)*linePt1.Y + q*linePt2.Y}
}

func slopesNearColinear(pt1, pt2, pt3 *Point, distSqrd int) bool {
	if distanceSqrd(pt1, pt2) > distanceSqrd(pt1, pt3) {
		return false
	}
	cpol := closestPointOnLine(pt2, pt1, pt3)
	dx := pt2.X - cpol.X
	dy := pt2.Y - cpol.Y
	return (dx*dx + dy*dy) < distSqrd
}

func pointsAreClose(pt1, pt2 *Point, distSqrd int) bool {
	dx := pt1.X - pt2.X
	dy := pt1.Y - pt2.Y
	return (dx*dx)+(dy*dy) <= distSqrd
}

//defaults: distance float64 = 1.415
func CleanPolygon(poly []*Point, distance int) []*Point {
	distSqrd := distance * distance
	highI := len(poly) - 1
	for highI > 0 && pointsEqual(poly[highI], poly[0]) {
		highI -= 1
	}
	if highI < 2 {
		return nil
	}
	pt := poly[highI]
	result := make([]*Point, 0)
	i := 0
	for {
		for i < highI && pointsAreClose(pt, poly[i+1], distSqrd) {
			i += 2
		}
		i2 := i
		for i < highI && (pointsAreClose(poly[i], poly[i+1], distSqrd) ||
			slopesNearColinear(pt, poly[i], poly[i+1], distSqrd)) {
			i += 1
		}
		if i >= highI {
			break
		} else if i != i2 {
			continue
		}
		pt = poly[i]
		i += 1
		result = append(result, pt)
	}

	if i <= highI {
		result = append(result, poly[i])
	}
	j := len(result)
	if j > 2 && slopesNearColinear(result[j-2], result[j-1], result[0], distSqrd) {
		result = result[0 : j-1]
	}
	if len(result) < 3 {
		return nil
	} else {
		return result
	}
}

// defaults: distance float64 = 1.415
func CleanPolygons(polys [][]*Point, distance int) [][]*Point {
	result := make([][]*Point, 0)
	for _, poly := range polys {
		result = append(result, CleanPolygon(poly, distance))
	}
	return result
}

func SimplifyPolygon(poly []*Point, fillType PolyFillType) [][]*Point {
	var result [][]*Point
	c := NewClipper()
	c.ForceSimple = true
	c.AddPolygon(poly, Subject)
	c.Execute(Union, result, fillType, fillType)
	return result
}

func SimplifyPolygons(polys [][]*Point, fillType PolyFillType) [][]*Point {
	var result [][]*Point
	c := NewClipper()
	c.ForceSimple = true
	c.AddPolygons(polys, Subject)
	c.Execute(Union, result, fillType, fillType)
	return result
}

type polySorter []*Point

func (p polySorter) Len() int           { return len(p) }
func (p polySorter) Less(i, j int) bool { return false }
func (p polySorter) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

// convert float to int (rounding)
func round(f float64) int {
	return int(f + 0.5)
}
func intAbs(i int) int {
	if i > 0 {
		return i
	} else {
		return i * -1
	}
}
func min(a, b int) int {
	if a < b {
		return a
	} else {
		return b
	}
}
func max(a, b int) int {
	if a > b {
		return a
	} else {
		return b
	}
}
