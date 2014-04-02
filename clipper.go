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
	Miter
)

type EndType int

const (
	Closed EndType = iota
	ButtEnd
	SquareEnd
	Round
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
type rect struct{ left, top, right, bottom float64 }

type localMinima struct {
	y, leftBound, rightBound, nextLm float64
}

type scanbeam struct {
	y      float64
	nextSb *scanbeam
}

func (sb *scanbeam) String() string {
	s := "nil"
	if self.nextSb != nil {
		s = "<obj>"
	}
	return fmt.Sprintf("(y:%d, nextSb:%s)", sb.y, s)
}

type intersectNode struct {
	e1, e2, pt float64
	nextIn     *intersectNode
}

type outPt struct {
	idx, pt        int
	prevOp, nextOp *outPt
}

type outRec struct {
	idx       int
	bottomPt  *Point
	isHole    bool
	FirstLeft int
	pts       int
	PolyNode  int
}

type joinRec struct {
	pt1a, pt1b, poly1Idx, pt2a, pt2b, poly2Idx int
}

type horzJoin struct {
	edge           int
	savedIdx       int
	prevHj, nextHj *horzJoin
}

//===============================================================================
// Unit global functions ...
//===============================================================================

func IntsToPoints(ints []int) []*Point {
	result := make([]*Point, len(ints)/2)
	for i := 0; i < len(ints); i += 2 {
		result[i] = Point{ints[i], ints[i+1]}
	}
	return result
}

// see http://www.mathopenref.com/coordpolygonarea2.html
func Area(polygon []*Point) int {
	highI := len(polygon) - 1
	A = (polygon[highI].x + polygon[0].x) * (polygon[0].y - polygon[highI].y)
	for i := 0; i < highI; i++ {
		A += (polygon[i].x + polygon[i+1].x) * (polygon[i+1].y - polygon[i].y)
	}
	return A / 2
}

func orientation(polygon []*FloatPoint) bool {
	return Area(polygon) > 0.0
}

//===============================================================================
// PolyNode & PolyTree classes (+ ancilliary functions)
//===============================================================================

// Node of PolyTree
type PolyNode struct {
	Contour           []*FloatPoint
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

func (self *PolyNode) getNextSiblingUp() {
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

func (self *PolyTree) Clear() {
	self.allNodes = nil
	self.Childs = nil
	self.ChildCount = 0
}

func (self *PolyTree) GetFirst() {
	if self.ChildCount > 0 {
		return self.Childs[0]
	} else {
		return nil
	}
}

func (self *PolyTree) Total() {
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

func PolyTreeToPolygons(polyTree *PolyTree) {
	result = make([][]*Point, 0)
	addPolyNodeToPolygons(polyTree, result)
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
	self := new(Edge)
	self.Bot = Point(0, 0)
	self.Curr = Point(0, 0)
	self.Top = Point(0, 0)
	self.Delta = Point(0, 0)
	self.outIdx = -1
	self.polyType = Subject
	self.side = leftSide
	return self
}

func (self *edge) String() string {
	return fmt.Sprintf("(%d,%d . %d,%d {dx:%0.2f} %i)",
		self.Bot.x, self.Bot.y, self.Top.x, self.Top.y, self.dx, self.outIdx)
}

//===============================================================================
// ClipperBase class (+ data structs & ancilliary functions)
//===============================================================================

func pointsEqual(pt1, pt2 *Point) bool {
	return (pt1.x == pt2.x) && (pt1.y == pt2.y)
}

func slopesEqual(pt1, pt2, pt3, pt4 *Point) bool {
	if pt4 == nil {
		return (pt1.y-pt2.y)*(pt2.x-pt3.x) == (pt1.x-pt2.x)*(pt2.y-pt3.y)
	} else {
		return (pt1.y-pt2.y)*(pt3.x-pt4.x) == (pt1.x-pt2.x)*(pt3.y-pt4.y)
	}
}

func slopesEqual2(e1, e2 *edge) bool {
	return e1.Delta.y*e2.Delta.x == e1.Delta.x*e2.Delta.y
}

func setDx(e *edge) {
	e.Delta = Point{e.Top.x - e.Bot.x, e.Top.y - e.Bot.y}
	if e.Delta.y == 0 {
		e.dx = horizontal
	} else {
		e.dx = float(e.Delta.x) / float(e.Delta.y)
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

func initEdge(e, eNext, ePrev *edge, pt *Point, polyType *PolyType) {
	e.nextE = eNext
	e.prevE = ePrev
	e.Curr = pt
	if e.Curr.y >= e.nextE.Curr.y {
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
	e.PolyType = polyType
}

func swapX(e *edge) {
	e.Curr = Point{e.Top.x, e.Curr.y}
	e.Top = Point{e.Bot.x, e.Top.y}
	e.Bot = Point{e.Curr.x, e.Bot.y}
}

type ClipperBase struct {
	EdgeList      [][]*edge // 2D array
	LocalMinList  *edge     // single-linked list of LocalMinima
	CurrentLocMin *edge
}

func (self *ClipperBase) insertLocalMinima(lm *edge) {
	if self._LocalMinList == nil {
		self._LocalMinList = lm
	} else if lm.y >= self._LocalMinList.y {
		lm.nextLm = self._LocalMinList
		self._LocalMinList = lm
	} else {
		tmp := self._LocalMinList
		for tmp.nextLm != nil && lm.y < tmp.nextLm.y {
			tmp = tmp.nextLm
		}
		lm.nextLm = tmp.nextLm
		tmp.nextLm = lm
	}
}

func (self *ClipperBase) addBoundsToLML(e *edge) {
	e.nextInLML = nil
	e = e.nextE
	for {
		if e.dx == horizontal {
			if (e.nextE.Top.y < e.Top.y) && (e.nextE.Bot.x > e.prevE.Bot.x) {
				break
			}
			if e.Top.x != e.prevE.Bot.x {
				swapX(e)
			}
			e.nextInLML = e.prevE
		} else if e.Bot.y == e.prevE.Bot.y {
			break
		} else {
			e.nextInLML = e.prevE
		}
		e = e.nextE
	}

	if e.dx == horizontal {
		if e.Bot.x != e.prevE.Bot.x {
			swapX(e)
		}
		lm = LocalMinima(e.prevE.Bot.y, e.prevE, e)
	} else if e.dx < e.prevE.dx {
		lm = LocalMinima(e.prevE.Bot.y, e.prevE, e)
	} else {
		lm = LocalMinima(e.prevE.Bot.y, e, e.prevE)
	}
	lm.leftBound.side = EdgeSide.Left
	lm.rightBound.side = EdgeSide.Right
	self.insertLocalMinima(lm)
	for {
		if e.nextE.Top.y == e.Top.y && e.nextE.dx != horizontal {
			break
		}
		e.nextInLML = e.nextE
		e = e.nextE
		if e.dx == horizontal && e.Bot.x != e.prevE.Top.x {
			swapX(e)
		}
	}
	return e.nextE
}

func (self *ClipperBase) reset() {
	lm := self._LocalMinList
	if lm != nil {
		self._CurrentLocMin = lm
	}
	for lm != nil {
		e := lm.leftBound
		for e != nil {
			e.Curr = e.Bot
			e.side = EdgeSide.Left
			e.outIdx = -1
			e = e.nextInLML
		}
		e = lm.rightBound
		for e != nil {
			e.Curr = e.Bot
			e.side = EdgeSide.Right
			e.outIdx = -1
			e = e.nextInLML
		}
		lm = lm.nextLm
	}
}

func (self *ClipperBase) AddPolygon(polygon []*Point, polyType *PolyType) bool {
	ln := len(polygon)
	if ln < 3 {
		return false
	}
	pg := polygon[:]
	j := 0
	// remove duplicate points && co-linear points
	for i := 1; i < len(polygon); i++ {
		if pointsEqual(pg[j], polygon[i]) {
			continue
		} else if (j > 0) && slopesEqual(pg[j-1], pg[j], polygon[i]) {
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
		} else if pointsEqual(pg[0], pg[1]) || slopesEqual(pg[j], pg[0], pg[1]) {
			pg[0] = pg[j]
			j -= 1
		} else if slopesEqual(pg[j-1], pg[j], pg[0]) {
			j -= 1
		} else if slopesEqual(pg[0], pg[1], pg[2]) {
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
		if e.Top.y < eHighest.Top.y {
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

func (self *ClipperBase) Clear() {
	self.edgeList = make([]*edge, 0)
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
	var x, y float64
	if slopesEqual2(edge1, edge2) {
		if edge2.Bot.y > edge1.Bot.y {
			y = edge2.Bot.y
		} else {
			y = edge1.Bot.y
		}
		return Point(0, y), false
	}
	if edge1.dx == 0 {
		x = edge1.Bot.x
		if edge2.dx == horizontal {
			y = edge2.Bot.y
		} else {
			b2 := edge2.Bot.y - float(edge2.Bot.x)/edge2.dx
			y = round(float(x)/edge2.dx + b2)
		}
	} else if edge2.dx == 0 {
		x = edge2.Bot.x
		if edge1.dx == horizontal {
			y = edge1.Bot.y
		} else {
			b1 := edge1.Bot.y - float(edge1.Bot.x)/edge1.dx
			y = round(float(x)/edge1.dx + b1)
		}
	} else {
		b1 := float(edge1.Bot.x) - float(edge1.Bot.y)*edge1.dx
		b2 := float(edge2.Bot.x) - float(edge2.Bot.y)*edge2.dx
		m := (b2 - b1) / (edge1.dx - edge2.dx)
		y = round(m)
		if math.Abs(edge1.dx) < math.Abs(edge2.dx) {
			x = round(edge1.dx*m + b1)
		} else {
			x = round(edge2.dx*m + b2)
		}
	}
	if (y < edge1.Top.y) || (y < edge2.Top.y) {
		if edge1.Top.y > edge2.Top.y {
			return edge1.Top, topX(edge2, edge1.Top.y) < edge1.Top.x
		} else {
			return edge2.Top, topX(edge1, edge2.Top.y) > edge2.Top.x
		}
	} else {
		return Point(x, y), true
	}
}

func topX(e *edge, currentY float64) float64 {
	if currentY == e.Top.y {
		return e.Top.x
	} else if e.Top.x == e.Bot.x {
		return e.Bot.x
	} else {
		return e.Bot.x + round(e.dx*float(currentY-e.Bot.y))
	}
}

func e2InsertsBeforeE1(e1, e2 *edge) bool {
	if e2.Curr.x == e1.Curr.x {
		if e2.Top.y > e1.Top.y {
			return e2.Top.x < _TopX(e1, e2.Top.y)
		}
		return e1.Top.x > _TopX(e2, e1.Top.y)
	} else {
		return e2.Curr.x < e1.Curr.x
	}
}

func isMinima(e *edge) bool {
	return (e != nil) && (e.prevE.nextInLML != e) && (e.nextE.nextInLML != e)
}

func isMaxima(e *edge, y float64) bool {
	return (e != nil) && (e.Top.y == y) && (e.nextInLML == nil)
}

func isIntermediate(e *edge, y float64) bool {
	return e.Top.y == y && e.nextInLML != nil
}

func getMaximaPair(e *edge) {
	if !isMaxima(e.nextE, e.Top.y) || e.nextE.Top.x != e.Top.x {
		return e.prevE
	} else {
		return e.nextE
	}
}

func getnextInAEL(e *edge, dir *direction) {
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

func getDx(pt1, pt2 *Point) {
	if pt1.y == pt2.y {
		return horizontal
	} else {
		return float(pt2.x-pt1.x) / float(pt2.y-pt1.y)
	}
}

func param1RightOfParam2(outRec1, outRec2 *edge) bool {
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
	dx1p := abs(getDx(btmPt1.pt, p.pt))
	p = btmPt1.nextOp
	for pointsEqual(p.pt, btmPt1.pt) && (p != btmPt1) {
		p = p.nextOp
	}
	dx1n := abs(getDx(btmPt1.pt, p.pt))

	p = btmPt2.prevOp
	for pointsEqual(p.pt, btmPt2.pt) && (p != btmPt2) {
		p = p.prevOp
	}
	dx2p := abs(getDx(btmPt2.pt, p.pt))
	p = btmPt2.nextOp
	for pointsEqual(p.pt, btmPt2.pt) && (p != btmPt2) {
		p = p.nextOp
	}
	dx2n = abs(getDx(btmPt2.pt, p.pt))
	return (dx1p >= dx2p && dx1p >= dx2n) || (dx1n >= dx2p && dx1n >= dx2n)
}

func getBottomPt(pp *outPt) *Point {
	var dups *outPt
	p = pp.nextOp
	for p != pp {
		if p.pt.y > pp.pt.y {
			pp = p
			dups = nil
		} else if p.pt.y == pp.pt.y && p.pt.x <= pp.pt.x {
			if p.pt.x < pp.pt.x {
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

func getLowermostRec(outRec1, outRec2 *outPt) {
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
	if outPt1.pt.y > outPt2.pt.y {
		return outRec1
	} else if outPt1.pt.y < outPt2.pt.y {
		return outRec2
	} else if outPt1.pt.x < outPt2.pt.x {
		return outRec1
	} else if outPt1.pt.x > outPt2.pt.x {
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

func setHoleState(e *edge, outRec *outPt, polyOutList [][]*Point) {
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

func pointCount(pts []*Point) int {
	if pts == nil {
		return 0
	}
	p := pts
	result = 0
	for {
		result++
		p = p.nextOp
		if p == pts {
			break
		}
	}
	return result
}

func pointIsVertex(pt *Point, outPts []*outPt) bool {
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

func reversePolyPtLinks(pp *Point) {
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

func fixupOutPolygon(outRec []*Point) {
	lastOK = nil
	outRec.bottomPt = nil
	pp = outRec.pts
	for {
		if pp.prevOp == pp || pp.nextOp == pp.prevOp {
			outRec.pts = nil
			return
		}
		if pointsEqual(pp.pt, pp.nextOp.pt) ||
			slopesEqual(pp.prevOp.pt, pp.pt, pp.nextOp.pt) {
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

func fixHoleLinkage(outRec []*Point) {
	if outRec.FirstLeft == nil ||
		(outRec.isHole != outRec.FirstLeft.isHole &&
			outRec.FirstLeft.pts != nil) {
		return
	}
	orfl := outRec.FirstLeft
	for orfl != nil &&
		(orfl.isHole == outRec.isHole || orfl.pts == nil) {
		orfl = orfl.FirstLeft
	}
	outRec.FirstLeft = orfl
}

func getOverlapSegment(pt1a, pt1b, pt2a, pt2b *Point) (*Point, *Point, bool) {
	// precondition: segments are co-linear
	var pt1, pt2 *Point
	if abs(pt1a.x-pt1b.x) > abs(pt1a.y-pt1b.y) {
		if pt1a.x > pt1b.x {
			tmp := pt1a
			pt1a = pt1b
			pt1b = tmp
		}
		if pt2a.x > pt2b.x {
			tmp := pt2a
			pt2a = pt2b
			pt2b = tmp
		}
		if pt1a.x > pt2a.x {
			pt1 = pt1a
		} else {
			pt1 = pt2a
		}
		if pt1b.x < pt2b.x {
			pt2 = pt1b
		} else {
			pt2 = pt2b
		}
		return pt1, pt2, pt1.x < pt2.x
	} else {
		if pt1a.y < pt1b.y {
			tmp := pt1a
			pt1a = pt1b
			pt1b = tmp
		}
		if pt2a.y < pt2b.y {
			tmp := pt2a
			pt2a = pt2b
			pt2b = tmp
		}
		if pt1a.y < pt2a.y {
			pt1 = pt1a
		} else {
			pt1 = pt2a
		}
		if pt1b.y > pt2b.y {
			pt2 = pt1b
		} else {
			pt2 = pt2b
		}
		return pt1, pt2, pt1.y > pt2.y
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
		if slopesEqual(pt1a, pt2a, outPt.pt, outPt.prevOp.pt) && slopesEqual(pt1a, pt2a, outPt.pt) {
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
	} else if pt1.x != pt2.x {
		return (pt1.x < pt3.x) == (pt3.x < pt2.x)
	} else {
		return (pt1.y < pt3.y) == (pt3.y < pt2.y)
	}
}

func insertPolyPtBetween(outPt1, outPt2 *outPt, pt *Point) *outPt {
	if outPt1 == outPt2 {
		panic("JoinError")
	}
	result := &OutPt{outPt1.idx, pt}
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
	return ((pt.x == linePt1.x) && (pt.y == linePt1.y)) ||
		((pt.x == linePt2.x) && (pt.y == linePt2.y)) ||
		(((pt.x > linePt1.x) == (pt.x < linePt2.x)) &&
			((pt.y > linePt1.y) == (pt.y < linePt2.y)) &&
			((pt.x-linePt1.x)*(linePt2.y-linePt1.y) ==
				(linePt2.x-linePt1.x)*(pt.y-linePt1.y)))
}

func PointOnPolygon(pt, pp *Point) bool {
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

func PointInPolygon(pt []*Point, outPt *outPt) bool {
	result := false
	outPt2 := outPt
	for {
		if (((outPt2.pt.y <= pt.y) && (pt.y < outPt2.prevOp.pt.y)) ||
			((outPt2.prevOp.pt.y <= pt.y) && (pt.y < outPt2.pt.y))) &&
			(pt.x < (outPt2.prevOp.pt.x-outPt2.pt.x)*(pt.y-outPt2.pt.y)/
				(outPt2.prevOp.pt.y-outPt2.pt.y)+outPt2.pt.x) {
			result = !result
		}
		outPt2 = outPt2.nextOp
		if outPt2 == outPt {
			break
		}
	}
	return result
}

func Poly2ContainsPoly1(outPt1, outPt2 *outPt) bool {
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

func EdgesAdjacent(inode *intersectNode) bool {
	return (inode.e1.nextInSEL == inode.e2) ||
		(inode.e1.prevInSEL == inode.e2)
}

func UpdateOutPtIdxs(outrec *outRec) {
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

	polyOutList    [][]*Point
	clipType       ClipType
	scanbeam       *edge
	activeEdges    []*edge
	sortedEdges    []*edge
	intersectNodes []*intersectNode
	clipFillType   PolyFillType
	subjFillType   PolyFillType
	executeLocked  bool
	usingPolyTree  bool
	joinList       []*edge
	horzJoins      []*edge
}

func NewClipper() *Clipper {
	self := new(Clipper)
	self.ReverseSolution = false
	self.ForceSimple = false
	self.polyOutList = make([][]*Point, 0)
	self.clipType = ClipType.Intersection
	self.clipFillType = PolyFillType.EvenOdd
	self.subjFillType = PolyFillType.EvenOdd
	return self
}

func (self *Clipper) reset() {
	ClipperBase.reset(self)
	self.scanbeam = nil
	self.polyOutList = make([][]*Point, 0)
	lm = self._LocalMinList
	for lm != nil {
		self._InsertScanbeam(lm.y)
		lm = lm.nextLm
	}
}

func (self *Clipper) Clear() {
	self.polyOutList = make([][]*Point, 0)
	ClipperBase.Clear(self)
}

func (self *Clipper) insertScanbeam(y float64) {
	if self.scanbeam == nil {
		self.scanbeam = Scanbeam(y)
	} else if y > self.scanbeam.y {
		self.scanbeam = Scanbeam(y, self.scanbeam)
	} else {
		sb := self.scanbeam
		for sb.nextSb != nil && y <= sb.nextSb.y {
			sb = sb.nextSb
		}
		if y == sb.y {
			return
		}
		newSb := Scanbeam(y, sb.nextSb)
		sb.nextSb = newSb
	}
}

func (self *Clipper) popScanbeam() *Scanbeam {
	result := self.scanbeam.y
	self.scanbeam = self.scanbeam.nextSb
	return result
}

func (self *Clipper) setWindingCount(edge *edge) {
	e := edge.prevInAEL
	for e != nil && e.PolyType != edge.PolyType {
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
			if abs(e.windCnt) > 1 {
				if e.windDelta*edge.windDelta < 0 {
					edge.windCnt = e.windCnt
				} else {
					edge.windCnt = e.windCnt + edge.windDelta
				}
			} else {
				edge.windCnt = e.windCnt + e.windDelta + edge.windDelta
			}
		} else if (abs(e.windCnt) > 1) && (e.windDelta*edge.windDelta < 0) {
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

func (self *Clipper) isEvenOddFillType(edge *edge) {
	if edge.PolyType == PolyType.Subject {
		return self._SubjFillType == PolyFillType.EvenOdd
	} else {
		return self._ClipFillType == PolyFillType.EvenOdd
	}
}

func (self *Clipper) isEvenOddAltFillType(edge *edge) {
	if edge.PolyType == PolyType.Subject {
		return self._ClipFillType == PolyFillType.EvenOdd
	} else {
		return self._SubjFillType == PolyFillType.EvenOdd
	}
}

func (self *Clipper) isContributing(edge *edge) {
	var pft, pft2 PolyFillType
	if edge.PolyType == PolyType.Subject {
		pft = self._SubjFillType
		pft2 = self._ClipFillType
	} else {
		pft = self._ClipFillType
		pft2 = self._SubjFillType
	}
	if pft == PolyFillType.EvenOdd || pft == PolyFillType.NonZero {
		if abs(edge.windCnt) != 1 {
			return false
		}
	} else if pft == PolyFillType.Positive {
		if edge.windCnt != 1 {
			return false
		}
	} else if pft == PolyFillType.Negative {
		if edge.windCnt != -1 {
			return false
		}
	}

	if self._ClipType == ClipType.Intersection { //////////////////////
		if pft2 == PolyFillType.EvenOdd || pft2 == PolyFillType.NonZero {
			return edge.windCnt2 != 0
		} else if pft2 == PolyFillType.Positive {
			return edge.windCnt2 > 0
		} else {
			return edge.windCnt2 < 0 // Negative
		}
	} else if self._ClipType == ClipType.Union { //////////////////////
		if pft2 == PolyFillType.EvenOdd || pft2 == PolyFillType.NonZero {
			return edge.windCnt2 == 0
		} else if pft2 == PolyFillType.Positive {
			return edge.windCnt2 <= 0
		} else {
			return edge.windCnt2 >= 0 // Negative
		}
	} else if self._ClipType == ClipType.Difference { //////////////////////
		if edge.PolyType == PolyType.Subject {
			if pft2 == PolyFillType.EvenOdd || pft2 == PolyFillType.NonZero {
				return edge.windCnt2 == 0
			} else if edge.PolyType == PolyFillType.Positive {
				return edge.windCnt2 <= 0
			} else {
				return edge.windCnt2 >= 0
			}
		} else {
			if pft2 == PolyFillType.EvenOdd || pft2 == PolyFillType.NonZero {
				return edge.windCnt2 != 0
			} else if pft2 == PolyFillType.Positive {
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
		edge.nextInSEL = self._SortedEdges
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
		e := self._ActiveEdges
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

func (self *Clipper) insertLocalMinimaIntoAEL(botY *edge) {
	for self.currentLocMin != nil &&
		self.currentLocMin.y == botY {
		lb := self._CurrentLocMin.leftBound
		rb := self._CurrentLocMin.rightBound
		self.insertEdgeIntoAEL(lb)
		self.insertScanbeam(lb.Top.y)
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
			self.insertScanbeam(rb.nextInLML.Top.y)
		} else {
			self.insertScanbeam(rb.Top.y)
		}
		if self.isContributing(lb) {
			self.addLocalMinPoly(lb, rb, Point(lb.Curr.x, self.currentLocMin.y))
		}

		if rb.outIdx >= 0 && rb.dx == horizontal && self._HorzJoins != nil {
			hj = self._HorzJoins
			for {
				dummy1, dummy2, overlap := getOverlapSegment(hj.edge.Bot, hj.edge.Top, rb.Bot, rb.Top)
				if overlap {
					self._AddJoin(hj.edge, rb, hj.savedIdx, -1)
				}
				hj = hj.nextHj
				if hj == self.horzJoins {
					break
				}
			}
		}
		if lb.nextInAEL != rb {
			if rb.outIdx >= 0 && rb.prevInAEL.outIdx >= 0 && _SlopesEqual2(rb.prevInAEL, rb) {
				self.addJoin(rb, rb.prevInAEL, -1, -1)
			}
			e := lb.nextInAEL
			pt := lb.Curr
			for e != rb {
				self._IntersectEdges(rb, e, pt, Protects.Neither)
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

func (self *Clipper) isTopHorz(xPos float64) bool {
	e = self.sortedEdges
	for e != nil {
		if (xPos >= min(e.Curr.x, e.Top.x)) && (xPos <= max(e.Curr.x, e.Top.x)) {
			return false
		}
		e = e.nextInSEL
	}
	return true
}

func (self *Clipper) processHorizontal(horzEdge *edge) {
	var horzLeft, horzRight *edge
	var direction *Direction
	if horzEdge.Curr.x < horzEdge.Top.x {
		horzLeft = horzEdge.Curr.x
		horzRight = horzEdge.Top.x
		direction = Direction.LeftToRight
	} else {
		horzLeft = horzEdge.Top.x
		horzRight = horzEdge.Curr.x
		direction = Direction.RightToLeft
	}
	eMaxPair := nil
	if horzEdge.nextInLML == nil {
		eMaxPair = _GetMaximaPair(horzEdge)
	}
	e := getnextInAEL(horzEdge, direction)
	for e != nil {
		if (e.Curr.x == horzEdge.Top.x) && eMaxPair == nil {
			if slopesEqual2(e, horzEdge.nextInLML) {
				if horzEdge.outIdx >= 0 && e.outIdx >= 0 {
					self.addJoin(horzEdge.nextInLML, e, horzEdge.outIdx, -1)
				}
				break
			} else if e.dx < horzEdge.nextInLML.dx {
				break
			}
		}
		eNext = getnextInAEL(e, direction)
		if eMaxPair != nil ||
			((direction == Direction.LeftToRight) && (e.Curr.x < horzRight)) ||
			((direction == Direction.RightToLeft) && (e.Curr.x > horzLeft)) {
			if e == eMaxPair {
				if direction == Direction.LeftToRight {
					self.intersectEdges(horzEdge, e, Point(e.Curr.x, horzEdge.Curr.y), protectsNeither)
				} else {
					self.intersectEdges(e, horzEdge, Point(e.Curr.x, horzEdge.Curr.y), protectsNeither)
				}
				return
			} else if e.dx == horizontal && !isMinima(e) && e.Curr.x <= e.Top.x {
				if direction == Direction.LeftToRight {
					self.intersectEdges(horzEdge, e, Point(e.Curr.x, horzEdge.Curr.y),
						protectRight(!self.isTopHorz(e.Curr.x)))
				} else {
					self.intersectEdges(e, horzEdge, Point(e.Curr.x, horzEdge.Curr.y),
						protectLeft(!self.isTopHorz(e.Curr.x)))
				}
			} else if direction == Direction.LeftToRight {
				self.intersectEdges(horzEdge, e, Point(e.Curr.x, horzEdge.Curr.y),
					protectRight(!self.isTopHorz(e.Curr.x)))
			} else {
				self.intersectEdges(e, horzEdge, Point(e.Curr.x, horzEdge.Curr.y),
					protectLeft(!self.isTopHorz(e.Curr.x)))
			}
			self.swapPositionsInAEL(horzEdge, e)
		} else if (direction == Direction.LeftToRight && e.Curr.x >= horzRight) ||
			(direction == Direction.RightToLeft && e.Curr.x <= horzLeft) {
			break
		}
		e = eNext
	}
	if horzEdge.nextInLML != nil {
		if horzEdge.outIdx >= 0 {
			self._AddOutPt(horzEdge, horzEdge.Top)
		}
		self.updateEdgeIntoAEL(horzEdge)
	} else {
		if horzEdge.outIdx >= 0 {
			self._IntersectEdges(horzEdge, eMaxPair,
				Point(horzEdge.Top.x, horzEdge.Curr.y), protectsBoth)
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
		e = self.sortedEdges
		self.deleteFromSEL(e)
		self.processHorizontal(e)
	}
}

func (self *Clipper) addJoin(e1, e2 *edge, e1OutIdx, e2OutIdx int) {
	jr := JoinRec()
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
		self.joinList = make([]*JoinRec, 0)
	}
	self.joinList = append(self.joinList, jr)
}

func (self *Clipper) fixupJoinRecs(jr *JoinRec, outPt *outPt, startIdx int) {
	for i := startIdx; i < len(self._JoinList); i++ {
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
	hj := HorzJoin(e, idx)
	if self._HorzJoins == nil {
		self._HorzJoins = hj
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
	newNode := IntersectNode(e1, e2, pt)
	if self.intersectNodes == nil {
		self.intersectNodes = newNode
	} else if newNode.pt.y > self._IntersectNodes.pt.y {
		newNode.nextIn = self._IntersectNodes
		self.intersectNodes = newNode
	} else {
		node = self._IntersectNodes
		for node.nextIn != nil &&
			newNode.pt.y < node.nextIn.pt.y {
			node = node.nextIn
		}
		newNode.nextIn = node.nextIn
		node.nextIn = newNode
	}
}

func (self *Clipper) processIntersections(botY, topY float64) bool {
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

func (self *Clipper) buildIntersectList(botY, topY *edge) {
	e := self.activeEdges
	if e == nil {
		return
	}
	self.sortedEdges = e
	for e != nil {
		e.prevInSEL = e.prevInAEL
		e.nextInSEL = e.nextInAEL
		e.Curr = Point(_TopX(e, topY), e.Curr.y)
		e = e.nextInAEL
	}
	for {
		isModified := false
		e = self.sortedEdges
		for e.nextInSEL != nil {
			eNext := e.nextInSEL
			if e.Curr.x <= eNext.Curr.x {
				e = eNext
				continue
			}
			pt, intersected = intersectPoint(e, eNext)
			if !intersected && e.Curr.x > eNext.Curr.x+1 {
				panic("Intersect Error")
			}
			if pt.y > botY {
				pt = Point(_TopX(e, botY), botY)
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
		e1.Top.x == pt.x && e1.Top.y == pt.y
	e2stops := protect&protectsRight == 0 &&
		e2.nextInLML == nil &&
		e2.Top.x == pt.x && e2.Top.y == pt.y
	e1Contributing := e1.outIdx >= 0
	e2contributing := e2.outIdx >= 0

	if e1.PolyType == e2.PolyType {
		if self._IsEvenOddFillType(e1) {
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
	var e1FillType, e1FillType2 FillType
	if e1.PolyType == PolyType.Subject {
		e1FillType = self.subjFillType
		e1FillType2 = self.clipFillType
	} else {
		e1FillType = self.clipFillType
		e1FillType2 = self.subjFillType
	}
	if e2.PolyType == PolyType.Subject {
		e2FillType = self.subjFillType
		e2FillType2 = self.clipFillType
	} else {
		e2FillType = self.clipFillType
		e2FillType2 = self.subjFillType
	}
	if e1FillType == PolyFillType.Positive {
		e1Wc = e1.windCnt
	} else if e1FillType == PolyFillType.Negative {
		e1Wc = -e1.windCnt
	} else {
		e1Wc = abs(e1.windCnt)
	}

	if e2FillType == PolyFillType.Positive {
		e2Wc = e2.windCnt
	} else if e2FillType == PolyFillType.Negative {
		e2Wc = -e2.windCnt
	} else {
		e2Wc = abs(e2.windCnt)
	}

	if e1Contributing && e2contributing {
		if e1stops || e2stops ||
			(e1Wc != 0 && e1Wc != 1) || (e2Wc != 0 && e2Wc != 1) ||
			(e1.PolyType != e2.PolyType && self._ClipType != ClipType.Xor) {
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

		e1FillType2, e2FillType2 = PolyFillType.EvenOdd, PolyFillType.EvenOdd
		if e1FillType2 == PolyFillType.Positive {
			e1Wc2 = e1.windCnt2
		} else if e1FillType2 == PolyFillType.Negative {
			e1Wc2 = -e1.windCnt2
		} else {
			e1Wc2 = abs(e1.windCnt2)
		}
		if e2FillType2 == PolyFillType.Positive {
			e2Wc2 = e2.windCnt2
		} else if e2FillType2 == PolyFillType.Negative {
			e2Wc2 = -e2.windCnt2
		} else {
			e2Wc2 = abs(e2.windCnt2)
		}

		if e1.PolyType != e2.PolyType {
			self.addLocalMinPoly(e1, e2, pt)
		} else if e1Wc == 1 && e2Wc == 1 {
			if self._ClipType == ClipType.Intersection {
				if e1Wc2 > 0 && e2Wc2 > 0 {
					self.addLocalMinPoly(e1, e2, pt)
				}
			} else if self.clipType == ClipType.Union {
				if e1Wc2 <= 0 && e2Wc2 <= 0 {
					self._AddLocalMinPoly(e1, e2, pt)
				}
			} else if self._ClipType == ClipType.Difference {
				if (e1.PolyType == PolyType.Clip && e1Wc2 > 0 && e2Wc2 > 0) ||
					(e1.PolyType == PolyType.Subject && e1Wc2 <= 0 && e2Wc2 <= 0) {
					self.addLocalMinPoly(e1, e2, pt)
				}
			} else {
				self.addLocalMinPoly(e1, e2, pt)
			}
		} else {
			swapSides(e1, e2, self.polyOutList)
		}
	}
	if e1stops != e2stops &&
		((e1stops && e1.outIdx >= 0) || (e2stops && e2.outIdx >= 0)) {
		swapSides(e1, e2, self._PolyOutList)
		swapPolyIndexes(e1, e2)
	}
	if e1stops {
		self.deleteFromAEL(e1)
	}
	if e2stops {
		self.deleteFromAEL(e2)
	}
}

func (self *Clipper) doMaxima(e *edge, topY float64) {
	eMaxPair := getMaximaPair(e)
	x := e.Top.x
	eNext := e.nextInAEL
	for eNext != eMaxPair {
		if eNext == nil {
			panic("DoMaxima error")
		}
		self.intersectEdges(e, eNext, Point(x, topY), Protects.Both)
		self.swapPositionsInAEL(e, eNext)
		eNext = e.nextInAEL
	}
	if e.outIdx < 0 && eMaxPair.outIdx < 0 {
		self.deleteFromAEL(e)
		self.deleteFromAEL(eMaxPair)
	} else if e.outIdx >= 0 && eMaxPair.outIdx >= 0 {
		self.intersectEdges(e, eMaxPair, Point(x, topY), Protects.Neither)
	} else {
		panic("DoMaxima error")
	}
}

func (self *Clipper) updateEdgeIntoAEL(e *edge) {
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
		self.insertScanbeam(e.Top.y)
	}
	return e
}

func (self *Clipper) addLocalMinPoly(e1, e2 *edge, pt *Point) {
	var e *edge
	if e2.dx == horizontal || e1.dx > e2.dx {
		self.addOutPt(e1, pt)
		e2.outIdx = e1.outIdx
		e1.side = EdgeSide.Left
		e2.side = EdgeSide.Right
		e = e1
		if e.prevInAEL == e2 {
			prevE = e2.prevInAEL
		} else {
			prevE = e1.prevInAEL
		}
	} else {
		self.addOutPt(e2, pt)
		e1.outIdx = e2.outIdx
		e1.side = EdgeSide.Right
		e2.side = EdgeSide.Left
		e = e2
		if e.prevInAEL == e1 {
			prevE = e1.prevInAEL
		} else {
			prevE = e.prevInAEL
		}
	}
	if prevE != nil && prevE.outIdx >= 0 &&
		topX(prevE, pt.y) == topX(e, pt.y) &&
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

func (self *Clipper) createOutRec() *OutRec {
	outRec = OutRec(len(self.polyOutList))
	self.polyOutList = append(self.polyOutList, outRec)
	return outRec
}

func (self *Clipper) addOutPt(e *edge, pt *Point) {
	toFront := e.side == EdgeSide.Left
	if e.outIdx < 0 {
		outRec = self.createOutRec()
		e.outIdx = outRec.idx
		op := OutPt(outRec.idx, pt)
		op.nextOp = op
		op.prevOp = op
		outRec.pts = op
		setHoleState(e, outRec, self._PolyOutList)
	} else {
		outRec = self.polyOutList[e.outIdx]
		op := outRec.pts
		if (toFront && _PointsEqual(pt, op.pt)) ||
			(!toFront && _PointsEqual(pt, op.prevOp.pt)) {
			return
		}
		op2 := OutPt(outRec.idx, pt)
		op2.nextOp = op
		op2.prevOp = op.prevOp
		op.prevOp.nextOp = op2
		op.prevOp = op2
		if toFront {
			outRec.pts = op2
		}
	}
}

func (self *Clipper) appendPolygon(e1, e2 *edge) {
	outRec1 := self.polyOutList[e1.outIdx]
	outRec2 := self.polyOutList[e2.outIdx]
	var holeStateRec *OutRec
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
	newSide := EdgeSide.Left

	if e1.side == EdgeSide.Left {
		if e2.side == EdgeSide.Left {
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
		newSide = EdgeSide.Right
		if e2.side == EdgeSide.Right {
			// a b c z y x
			_ReversePolyPtLinks(p2_lft)
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
	OKIdx = outRec1.idx
	ObsoleteIdx = outRec2.idx

	e1.outIdx = -1
	e2.outIdx = -1

	e = self._ActiveEdges
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
			nextNode = inode.nextIn
			for nextNode && !edgesAdjacent(nextNode) {
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

func (self *Clipper) processEdgesAtTopOfScanbeam(topY float64) {
	e := self._ActiveEdges
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
				e.Curr = Point(topX(e, topY), topY)
				if self.ForceSimple && e.prevInAEL != nil &&
					e.prevInAEL.Curr.x == e.Curr.x &&
					e.outIdx >= 0 && e.prevInAEL.outIdx >= 0 {
					if intermediateVert {
						self.addOutPt(e.prevInAEL, Point(e.Curr.x, topY))
					} else {
						self.addOutPt(e, Point(e.Curr.x, topY))
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
			if ePrev != nil && ePrev.Curr.x == e.Bot.x &&
				(ePrev.Curr.y == e.Bot.y) && (e.outIdx >= 0) &&
				(ePrev.outIdx >= 0) && (ePrev.Curr.y > ePrev.Top.y) &&
				slopesEqual2(e, ePrev) {
				self.addOutPt(ePrev, e.Bot)
				self.addJoin(e, ePrev, -1, -1)
			} else if eNext != nil && (eNext.Curr.x == e.Bot.x) &&
				(eNext.Curr.y == e.Bot.y) && (e.outIdx >= 0) &&
				(eNext.outIdx >= 0) && (eNext.Curr.y > eNext.Top.y) &&
				slopesEqual2(e, eNext) {
				self.addOutPt(eNext, e.Bot)
				self.addJoin(e, eNext, -1, -1)
			}
		}
		e = e.nextInAEL
	}
}

// see http://www.mathopenref.com/coordpolygonarea2.html
func (self *Clipper) area(pts []*Point) {
	result := 0.0
	p := pts
	for {
		result += (p.pt.x + p.prevOp.pt.x) * (p.prevOp.pt.y - p.pt.y)
		p = p.nextOp
		if p == pts {
			break
		}
	}
	return result / 2.
}

func (self *Clipper) joinPoints(jr *JoinRecord) (*Point, *Point, bool) {
	var p1, p2 *Point
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

	prevOp = pp1a.prevOp
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

func (self *Clipper) fixupFirstLefts1(oldOutRec, newOutRec *OutRec) {
	for _, outRec := range self.polyOutList {
		if outRec.pts != nil && outRec.FirstLeft == oldOutRec {
			if poly2ContainsPoly1(outRec.pts, newOutRec.pts) {
				outRec.FirstLeft = newOutRec
			}
		}
	}
}

func (self *Clipper) fixupFirstLefts2(oldOutRec, newOutRec *OutRec) {
	for _, outRec := range self.polyOutList {
		if outRec.FirstLeft == oldOutRec {
			outRec.FirstLeft = newOutRec
		}
	}
}

func (self *Clipper) getOutRec(idx int) *OutRec {
	outrec := self.polyOutList[idx]
	for outrec != self._PolyOutList[outrec.idx] {
		outrec = self._PolyOutList[outrec.idx]
	}
	return outrec
}

func (self *Clipper) joinCommonEdges() {
	for i := 0; i < len(self.joinList); i++ {
		jr := self._JoinList[i]
		outRec1 := self._GetOutRec(jr.poly1Idx)
		outRec2 := self._GetOutRec(jr.poly2Idx)
		if outRec1.pts == nil || outRec2.pts == nil {
			continue
		}

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

				if (outRec2.isHole ^ self.ReverseSolution) == self.area(outRec2) > 0.0 {
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

				if (outRec1.isHole ^ self.ReverseSolution) == self.area(outRec1) > 0.0 {
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
	for i < len(self._PolyOutList) {
		outrec := self._PolyOutList[i]
		i += 1
		op = outrec.pts
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
					outrec2 = self.createOutRec()
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
		if (outRec.isHole ^ self.ReverseSolution) == (self._Area(outRec.pts) > 0.0) {
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
	solution []*Point,
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
	solutionTree []*Point,
	subjFillType,
	clipFillType PolyFillType) bool {
	if self.executeLocked {
		return false
	}
	defer func() {
		self._ExecuteLocked = false
		self._UsingPolyTree = false
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
		poly := make([][]*Point, 0)
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

		// add nodes to _AllNodes list ...
		polyNode := PolyNode()
		polyTree.allNodes = append(polyTree.allNodes, polyNode)
		outRec.PolyNode = polyNode
		op = outRec.pts
		for {
			polyNode.Contour = append(polyNode.Countour, op.pt)
			op = op.prevOp
			if op == outRec.pts {
				break
			}
		}
	}
	// build the tree ...
	for _, outRec := range self.polyOutList {
		if outRec.PolyNode == nil {
			continue
		}
		if outRec.FirstLeft == nil {
			polyTree.addChild(outRec.PolyNode)
		} else {
			outRec.FirstLeft.PolyNode.addChild(outRec.PolyNode)
		}
	}
	return
}

//===============================================================================
// OffsetPolygons (+ ancilliary functions)
//===============================================================================

func getUnitNormal(pt1, pt2 *Point) *FloatPoint {
	if pt2.x == pt1.x && pt2.y == pt1.y {
		return &FloatPoint{0.0, 0.0}
	}
	dx := float(pt2.x - pt1.x)
	dy := float(pt2.y - pt1.y)
	f := 1.0 / math.Hypot(dx, dy)
	dx = float(dx) * f
	dy = float(dy) * f
	return &FloatPoint{dy, -dx}
}

func getBounds(pts []*Point) *rect {
	var left, top, right, bottom *edge
	for _, poly := range pts {
		for _, pt := range poly {
			left = pt.x
			top = pt.y
			right = pt.x
			bottom = pt.y
			break
		}
		break
	}

	for _, poly := range pts {
		for _, pt := range poly {
			if pt.x < left {
				left = pt.x
			}
			if pt.x > right {
				right = pt.x
			}
			if pt.y < top {
				top = pt.y
			}
			if pt.y > bottom {
				bottom = pt.y
			}
		}
	}
	if left == nil {
		return Rect(0, 0, 0, 0)
	} else {
		return Rect(left, top, right, bottom)
	}
}

func getLowestPt(poly []*Point) *Point {
	// precondition: poly must not be empty
	result := poly[0]
	for _, pt := range poly {
		if pt.y > result.y || (pt.y == result.y && pt.x < result.x) {
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
	i = len(poly) - 1
	for i > 0 && pointsEqual(poly[i], poly[0]) {
		poly = append(poly[:i], poly[i+1:]...) // remove item i
		i -= 1
	}
	return poly
}

func offsetInternal(polys [][]*Point, isPolygon bool, delta int, jointype JoinType, endtype EndType, limit float64) []*Point {

	doSquare := func(pt *Point) {
		// see offset_triginometry.svg in the documentation folder ...
		dx := math.Tan(math.Atan2(sinA,
			Normals[k].x*Normals[j].x+Normals[k].y*Normals[j].y) / 4)
		result = append(result, Point(
			round(pt.x+delta*(Normals[k].x-Normals[k].y*dx)),
			round(pt.y+delta*(Normals[k].y+Normals[k].x*dx))))
		result = append(result, Point(
			round(pt.x+delta*(Normals[j].x+Normals[j].y*dx)),
			round(pt.y+delta*(Normals[j].y-Normals[j].x*dx))))
		return
	}

	doMiter := func(pt *Point, r float64) {
		q := delta / r
		result = append(result, Point(
			round(pt.x+(Normals[k].x+Normals[j].x)*q),
			round(pt.y+(Normals[k].y+Normals[j].y)*q)))
		return
	}

doRound:
	+func(pt) {
		a = math.Atan2(sinA,
			Normals[k].x*Normals[j].x+Normals[k].y*Normals[j].y)
		steps = round(step360 * abs(a))
		X, Y = Normals[k].x, Normals[k].y
		for i := 0; i < steps; i++ {
			result = append(result, Point(
				round(pt.x+X*delta), round(pt.y+Y*delta)))
			X2 = X
			X = X*mcos - msin*Y
			Y = X2*msin + Y*mcos
		}
		result = append(result, Point(round(pt.x+Normals[j].x*delta),
			round(pt.y+Normals[j].y*delta)))
		return
	}

	GetSin := func() float64 {
		result := (Normals[k].x*Normals[j].y - Normals[j].x*Normals[k].y)
		if result > 1.0 {
			result = 1.0
		} else if result < -1.0 {
			result = -1.0
		}
		return result
	}

	offsetPoint := func(jointype JoinType) (j int) {
		if sinA*delta < 0 {
			result = append(result, Point(round(pts[j].x+Normals[k].x*delta),
				round(pts[j].y+Normals[k].y*delta)))
			result = append(result, pts[j])
			result = append(result, Point(round(pts[j].x+Normals[j].x*delta),
				round(pts[j].y+Normals[j].y*delta)))
		} else if jointype == JoinType.Miter {
			r := 1.0 + (Normals[j].x*Normals[k].x + Normals[j].y*Normals[k].y)
			if r >= miterLim {
				doMiter(pts[j], r)
			} else {
				doSquare(pts[j])
			}
		} else if jointype == JoinType.Square {
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

	if jointype == JoinType.Miter {
		// miterLim: see offset_triginometry3.svg in the documentation folder ...
		if limit > 2 {
			miterLim = 2 / (limit * limit)
		} else {
			miterLim = 0.5
		}
		if endtype == EndType.Round {
			limit = 0.25
		}
	}

	if jointype == JoinType.Round || endtype == EndType.Round {
		if limit <= 0 {
			limit = 0.25
		} else if limit > abs(delta)*0.25 {
			limit = abs(delta) * 0.25
		}
		// step360: see offset_triginometry2.svg in the documentation folder ...
		step360 := math.Pi / math.Acos(1-limit/abs(delta))
		msin := math.Sin(2 * math.Pi / step360)
		mcos := math.Cos(2 * math.Pi / step360)
		step360 /= math.Pi * 2
		if delta < 0 {
			msin = -msin
		}
	}

	res := make([]*Point, 0)
	ppts := polys[:]
	for _, pts := range ppts {
		Normals = make([]*Point, 0)
		result = make([]*Point, 0)
		cnt := len(pts)

		if cnt == 0 || cnt < 3 && delta <= 0 {
			continue
		}

		if cnt == 1 {
			if jointype == JoinType.Round {
				X, Y := 1.0, 0.0
				for i := 0; i < round(step360*2*math.Pi); i++ {
					result = append(result, Point(round(pts[0].x+X*delta),
						round(pts[0].y+Y*delta)))
					X2 := X
					X = X*mcos - msin*Y
					Y = X2*msin + Y*mcos
				}
			} else {
				X, Y := -1.0, -1.0
				for i := 0; i < 4; i++ {
					result = append(result, Point(round(pts[0].x+X*delta),
						round(pts[0].y+Y*delta)))
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
			k := cnt - 1
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
				res = append(res, sort.Reverse(result).([]*Point))
			}

		} else {
			// offset the polyline going forward ...
			k = 0
			for j := 1; j < cnt-1; j++ {
				sinA = GetSin()
				k = offsetPoint(jointype)
			}

			// handle the end (butt, round || square) ...
			if endtype == EndType.Butt {
				j = cnt - 1
				pt1 = Point(round(float(pts[j].x)+Normals[j].x*delta),
					round(float(pts[j].y)+Normals[j].y*delta))
				result = append(result, pt1)
				pt1 = Point(round(float(pts[j].x)-Normals[j].x*delta),
					round(float(pts[j].y)-Normals[j].y*delta))
				result = append(result, pt1)
			} else {
				j = cnt - 1
				k = cnt - 2
				Normals[j] = FloatPoint(-Normals[j].x, -Normals[j].y)
				if endtype == EndType.Square {
					doSquare(pts[j])
				} else {
					doRound(pts[j])
				}
			}

			// re-build Normals ...
			for j := cnt - 1; j > 0; j-- {
				Normals[j] = FloatPoint(-Normals[j-1].x, -Normals[j-1].y)
			}
			Normals[0] = FloatPoint(-Normals[1].x, -Normals[1].y)

			// offset the polyline going backward ...
			k = cnt - 1
			for j := cnt - 2; j > 0; j-- {
				sinA = GetSin()
				k = offsetPoint(jointype)
			}

			// finally handle the start (butt, round || square) ...
			if endtype == EndType.Butt {
				pt1 = Point(round(float(pts[0].x)-Normals[0].x*delta),
					round(float(pts[0].y)-Normals[0].y*delta))
				result = append(result, pt1)
				pt1 = Point(round(float(pts[0].x)+Normals[0].x*delta),
					round(float(pts[0].y)+Normals[0].y*delta))
				result = append(result, pt1)
			} else {
				j = 0
				k = 1
				if endtype == EndType.Square {
					doSquare(pts[0])
				} else {
					doRound(pts[0])
				}
			}
			res.append(result)
		}
	}
	c := NewClipper()
	c.AddPolygons(res, PolyType.Subject)
	if delta > 0 {
		c.Execute(ClipType.Union, res, PolyFillType.Positive, PolyFillType.Positive)
	} else {
		bounds = getBounds(res)
		outer = make([]*Point, 0)
		outer = append(outer, Point(bounds.left-10, bounds.bottom+10))
		outer = append(outer, Point(bounds.right+10, bounds.bottom+10))
		outer = append(outer, Point(bounds.right+10, bounds.top-10))
		outer = append(outer, Point(bounds.left-10, bounds.top-10))
		c.AddPolygon(outer, PolyType.Subject)
		c.ReverseSolution = true
		c.Execute(ClipType.Union, res, PolyFillType.Negative, PolyFillType.Negative)
		if len(res) > 0 {
			res = res[1:len(res)]
		}
	}
	return res
}

// defaults: jointype = JoinType.Square, limit = 0.0, autoFix = true
func OffsetPolygons(polys [][]*Point, delta float64, jointype JoinType, limit float64, autoFix bool) [][]*Point {
	if !autoFix {
		return offsetInternal(polys, true, delta, jointype, EndType.Butt, limit)
	}
	pts := polys[:]
	botPoly := make([]*Point, 0)
	var botPt *Point
	for _, poly := range pts {
		poly = stripDupPts(poly)
		if len(poly) < 3 {
			continue
		}
		bot = getLowestPt(poly)
		if botPt == nil || (bot.y > botPt.y) ||
			(bot.y == botPt.y && bot.x < botPt.x) {
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
			pts[i] = sort.Reverse(pts[i]).([]*Point)
		}
	}
	return offsetInternal(pts, true, delta, jointype, EndType.Butt, limit)
}

//defaults: jointype = JoinType.Square, endtype = EndType.Square, limit = 0.0
func OffsetPolyLines(polys [][]*Point, delta float64, jointype JoinType, endtype EndType, limit float64) [][]*Point {
	polys2 := polys[:]
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

	if endtype == EndType.Closed {
		for i := 0; i < len(polys2); i++ {
			polys2 = append(polys2, sort.Reverse(polys2[i]).([]*Point))
		}
		return offsetInternal(polys2, true, delta, jointype, EndType.Butt, limit)
	} else {
		return offsetInternal(polys2, false, delta, jointype, endtype, limit)
	}
}

func distanceSqrd(pt1, pt2 *Point) int {
	dx := (pt1.x - pt2.x)
	dy := (pt1.y - pt2.y)
	return (dx*dx + dy*dy)
}

func closestPointOnLine(pt, linePt1, linePt2 *Point) *FloatPoint {
	dx := linePt2.x - linePt1.x
	dy := linePt2.y - linePt1.y
	if dx == 0 && dy == 0 {
		return FloatPoint(linePt1.x, linePt1.y)
	}
	q := ((pt.x-linePt1.x)*dx + (pt.Y-linePt1.Y)*dy) / (dx*dx + dy*dy)
	return FloatPoint(
		(1-q)*linePt1.X+q*linePt2.X,
		(1-q)*linePt1.Y+q*linePt2.Y)
}

func slopesNearColinear(pt1, pt2, pt3 *Point, distSqrd int) bool {
	if distanceSqrd(pt1, pt2) > distanceSqrd(pt1, pt3) {
		return false
	}
	cpol := closestPointOnLine(pt2, pt1, pt3)
	dx := pt2.x - cpol.x
	dy := pt2.y - cpol.y
	return (dx*dx + dy*dy) < distSqrd
}

func pointsAreClose(pt1, pt2 *Point, distSqrd int) bool {
	dx := pt1.x - pt2.x
	dy := pt1.y - pt2.y
	return (dx*dx)+(dy*dy) <= distSqrd
}

//defaults: distance float64 = 1.415
func CleanPolygon(poly []*Point, distance float64) []*Point {
	distSqrd := distance * distance
	highI := len(poly) - 1
	for highI > 0 && pointsEqual(poly[highI], poly[0]) {
		highI -= 1
	}
	if highI < 2 {
		return nil
	}
	pt := poly[highI]
	result = make([]*Point, 0)
	i = 0
	for {
		for i < highI && _PointsAreClose(pt, poly[i+1], distSqrd) {
			i += 2
		}
		i2 = i
		for i < highI && (_PointsAreClose(poly[i], poly[i+1], distSqrd) ||
			_SlopesNearColinear(pt, poly[i], poly[i+1], distSqrd)) {
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
	j = len(result)
	if j > 2 && _SlopesNearColinear(result[j-2], result[j-1], result[0], distSqrd) {
		result = result[0 : j-1]
	}
	if len(result) < 3 {
		return nil
	} else {
		return result
	}
}

// defaults: distance float64 = 1.415
func CleanPolygons(polys [][]*Point, distance float64) [][]*Point {
	result = make([][]*Point, 0)
	for _, poly := range polys {
		result = append(result, CleanPolygon(poly, distance))
	}
	return result
}

func SimplifyPolygon(poly []*Point, fillType FillType) []*Point {
	var result []*Point
	c := Clipper()
	c.ForceSimple = true
	c.AddPolygon(poly, PolyType.Subject)
	c.Execute(ClipType.Union, result, fillType, fillType)
	return result
}

func SimplifyPolygons(polys [][]*Point, fillType FillType) [][]*Point {
	var result [][]*Point
	c := Clipper()
	c.ForceSimple = true
	c.AddPolygons(polys, PolyType.Subject)
	c.Execute(ClipType.Union, result, fillType, fillType)
	return result
}
