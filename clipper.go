/*******************************************************************************
*                                                                              *
* Author    :  Angus Johnson                                                   *
* Version   :  6.1.5                                                           *
* Date      :  28 March 2014                                                   *
* Website   :  http://www.angusj.com                                           *
* Copyright :  Angus Johnson 2010-2014                                         *
*                                                                              *
* License:                                                                     *
* Use, modification & distribution is subject to Boost Software License Ver 1. *
* http://www.boost.org/LICENSE_1_0.txt                                         *
*                                                                              *
* Attributions:                                                                *
* The code in this library is an extension of Bala Vatti's clipping algorithm: *
* "A generic solution to polygon clipping"                                     *
* Communications of the ACM, Vol 35, Issue 7 (July 1992) pp 56-63.             *
* http://portal.acm.org/citation.cfm?id=129906                                 *
*                                                                              *
* Computer graphics and geometric modeling: implementation and algorithms      *
* By Max K. Agoston                                                            *
* Springer; 1 edition (January 4, 2005)                                        *
* http://books.google.com/books?q=vatti+clipping+agoston                       *
*                                                                              *
* See also:                                                                    *
* "Polygon Offsetting by Computing Winding Numbers"                            *
* Paper no. DETC2005-85513 pp. 565-575                                         *
* ASME 2005 International Design Engineering Technical Conferences             *
* and Computers and Information in Engineering Conference (IDETC/CIE2005)      *
* September 24-28, 2005 , Long Beach, California, USA                          *
* http://www.me.berkeley.edu/~mcmains/pubs/DAC05OffsetPolygon.pdf              *
*                                                                              *
*******************************************************************************/

/*******************************************************************************
*                                                                              *
* This is a translation of the Delphi Clipper library and the naming style     *
* used has retained a Delphi flavour.                                          *
*                                                                              *
*******************************************************************************/

//use_int32: When enabled 32bit ints are used instead of 64bit ints. This
//improve performance but coordinate values are limited to the range +/- 46340

package clipper

import (
	"fmt"
	"math"
	"math/big"
	"sort"
	"strings"
)

type Path []*IntPoint
type Paths []Path

func NewPath() Path {
	return Path(make([]*IntPoint, 0))
}
func NewPaths() Paths {
	return Paths(make([]Path, 0))
}

func (p Path) String() string {
	v := make([]string, len(p))
	for i, pp := range p {
		v[i] = pp.String()
	}
	return fmt.Sprintf("{%v}", strings.Join(v, ", "))
}

func (p Paths) String() string {
	s := "{"
	for i, pp := range p {
		s += pp.String()
		if i != len(p)-1 {
			s += "\n"
		}
	}
	s += "}"
	return s
}

type DoublePoint struct {
	X float64
	Y float64
}

func NewDoublePoint(x, y float64) *DoublePoint {
	dp := new(DoublePoint)
	dp.X = x
	dp.Y = y
	return dp
}

func CopyDoublePoint(dp *DoublePoint) *DoublePoint {
	dp2 := new(DoublePoint)
	dp2.X, dp2.Y = dp.X, dp.Y
	return dp2
}

func (ip *IntPoint) ToDoublePoint() *DoublePoint {
	dp := new(DoublePoint)
	dp.X = float64(ip.X)
	dp.Y = float64(ip.Y)
	return dp
}

//------------------------------------------------------------------------------
// PolyTree & PolyNode classes
//------------------------------------------------------------------------------

type PolyTree struct {
	PolyNode
	m_AllPolys []*PolyNode
}

func NewPolyTree() *PolyTree {
	pt := new(PolyTree)
	pt.m_AllPolys = make([]*PolyNode, 0)
	return pt
}

func (tree *PolyTree) toPolyNode() *PolyNode {
	node := new(PolyNode)
	node.m_Parent = tree.m_Parent
	node.m_polygon = tree.m_polygon
	node.m_Index = tree.m_Index
	node.m_jointype = tree.m_jointype
	node.m_endtype = tree.m_endtype
	node.m_Childs = tree.m_Childs
	node.IsOpen = tree.IsOpen
	return node
}

func (pt *PolyTree) Clear() {
	pt.m_AllPolys = make([]*PolyNode, 0)
	pt.m_Childs = make([]*PolyNode, 0)
}

func (pt *PolyTree) GetFirst() *PolyNode {
	if len(pt.m_Childs) > 0 {
		return pt.m_Childs[0]
	} else {
		return nil
	}
}

func (pt *PolyTree) Total() int {
	return len(pt.m_AllPolys)
}

type PolyNode struct {
	m_Parent   *PolyNode
	m_polygon  Path
	m_Index    int
	m_jointype JoinType
	m_endtype  EndType
	m_Childs   []*PolyNode
	IsOpen     bool
}

func NewPolyNode() *PolyNode {
	pn := new(PolyNode)
	pn.m_polygon = NewPath()
	pn.m_Childs = make([]*PolyNode, 0)
	return pn
}

func (pn *PolyNode) IsHoleNode() bool {
	result := true
	node := pn.m_Parent
	for node != nil {
		result = !result
		node = node.m_Parent
	}
	return result
}

func (pn *PolyNode) ChildCount() int {
	return len(pn.m_Childs)
}

func (pn *PolyNode) Contour() Path {
	return pn.m_polygon
}

func (pn *PolyNode) AddChild(Child *PolyNode) {
	cnt := pn.ChildCount()
	pn.m_Childs = append(pn.m_Childs, Child)
	Child.m_Parent = pn
	Child.m_Index = cnt
}

func (pn *PolyNode) GetNext() *PolyNode {
	if len(pn.m_Childs) > 0 {
		return pn.m_Childs[0]
	} else {
		return pn.GetNextSiblingUp()
	}
}

func (pn *PolyNode) GetNextSiblingUp() *PolyNode {
	if pn.m_Parent == nil {
		return nil
	} else if pn.m_Index == len(pn.m_Parent.m_Childs)-1 {
		return pn.m_Parent.GetNextSiblingUp()
	} else {
		return pn.m_Parent.m_Childs[pn.m_Index+1]
	}
}

func (pn *PolyNode) Childs() []*PolyNode {
	return pn.m_Childs
}

func (pn *PolyNode) Parent() *PolyNode {
	return pn.m_Parent
}

func (pn *PolyNode) IsHole() bool {
	return pn.IsHoleNode()
}

//------------------------------------------------------------------------------

func Int128Mul(lhs, rhs CInt) *big.Int {
	a := big.NewInt(int64(lhs))
	b := big.NewInt(int64(rhs))
	c := new(big.Int)
	return c.Mul(a, b)
}

//------------------------------------------------------------------------------

// The == operator, which
// will also compare Z values if they exist, is usually used instead of this one.
// This may be a problem when using Z values.
func (a *IntPoint) Equals(b *IntPoint) bool {
	return a.X == b.X && a.Y == b.Y
}

func (a *IntPoint) NotEqual(b *IntPoint) bool {
	return a.X != b.X || a.Y != b.Y
}

type IntRect struct {
	left, top, right, bottom CInt
}

func NewIntRect(l, t, r, b CInt) *IntRect {
	this := new(IntRect)
	this.left = l
	this.top = t
	this.right = r
	this.bottom = b
	return this
}

func (ir *IntRect) Copy() IntRect {
	return IntRect{
		left:   ir.left,
		top:    ir.top,
		right:  ir.right,
		bottom: ir.bottom,
	}
}

type ClipType int

const (
	CtIntersection ClipType = iota
	CtUnion
	CtDifference
	CtXor
)

type PolyType int

const (
	PtSubject PolyType = iota
	PtClip
)

//By far the most widely used winding rules for polygon filling are
//EvenOdd & NonZero (GDI, GDI+, XLib, OpenGL, Cairo, AGG, Quartz, SVG, Gr32)
//Others rules include Positive, Negative and ABS_GTR_EQ_TWO (only in OpenGL)
//see http://glprogramming.com/red/chapter11.html
type PolyFillType int

const (
	PftEvenOdd PolyFillType = iota
	PftNonZero
	PftPositive
	PftNegative
)

type JoinType int

const (
	JtSquare JoinType = iota
	JtRound
	JtMiter
)

type EndType int

const (
	EtClosedPolygon EndType = iota
	EtClosedLine
	EtOpenButt
	EtOpenSquare
	EtOpenRound
)

type EdgeSide int

const (
	EsLeft EdgeSide = iota
	EsRight
)

type Direction int

const (
	DRightToLeft Direction = iota
	DLeftToRight
)

type TEdge struct {
	Bot, Curr, Top, Delta            IntPoint
	Dx                               float64
	PolyTyp                          PolyType
	Side                             EdgeSide
	WindDelta                        int //1 or -1 depending on winding direction
	WindCnt                          int
	WindCnt2                         int //winding count of the opposite polytype
	OutIdx                           int
	Next, Prev, NextInLML, NextInAEL *TEdge
	PrevInAEL, NextInSEL, PrevInSEL  *TEdge
}

func (e *TEdge) String() string {
	return fmt.Sprintf("Bot: %v, Curr: %v, Top: %v, Delta: %v, Dx: %v",
		e.Bot, e.Curr, e.Top, e.Delta, e.Dx)
}

//func (e *TEdge) Copy() *TEdge {
//	o := new(TEdge)
//	o.Bot = e.Bot.Copy()
//	o.Curr = e.Curr.Copy()
//	o.Top = e.Top.Copy()
//	o.Delta = e.Delta.Copy()
//	o.Dx = e.Dx
//	o.PolyTyp = e.PolyTyp
//	o.Side = e.Side
//	o.WindDelta, o.WindCnt, o.WindCnt2, o.OutIdx =
//		e.WindDelta, e.WindCnt, e.WindCnt2, e.OutIdx
//	return o
//}

func (e *TEdge) printEdges() string {
	E := e
	s := ""
	for {
		s += fmt.Sprintf("%v\n", E)
		E = E.Next
		if E == nil || E == e {
			break
		}
	}
	return s
}

type IntersectNode struct {
	Edge1, Edge2 *TEdge
	Pt           *IntPoint
}

func (in *IntersectNode) String() string {
	return fmt.Sprintf("Edge1: %v, Edge2: %v, Pt: %v",
		in.Edge1, in.Edge2, in.Pt)
}

type IntersectNodeList []*IntersectNode

func (i IntersectNodeList) Len() int           { return len(i) }
func (i IntersectNodeList) Less(a, b int) bool { return i[a].Pt.Y < i[b].Pt.Y }
func (i IntersectNodeList) Swap(a, b int)      { i[a], i[b] = i[b], i[a] }

//  class MyIntersectNodeSort : IComparer<IntersectNode>
//  {
//    int Compare(IntersectNode node1, IntersectNode node2)
//    {
//      return (int)(node2.Pt.Y - node1.Pt.Y);
//    }
//  }

type LocalMinima struct {
	Y                     CInt
	LeftBound, RightBound *TEdge
	Next                  *LocalMinima
}

func (lm *LocalMinima) String() string {
	return fmt.Sprintf("Y: %v, LeftBound: %v, RightBound: %v",
		lm.Y, lm.LeftBound, lm.RightBound)
}

func (lm *LocalMinima) printLML() string {
	s := ""
	LM := lm
	for {
		s += fmt.Sprint(LM)
		if LM.Next == nil {
			break
		}
		LM = LM.Next
		s += "\n"
	}
	return s
}

type Scanbeam struct {
	Y    CInt
	Next *Scanbeam
}

func (s *Scanbeam) String() string {
	return fmt.Sprintf("Y: %v", s.Y)
}

func (s *Scanbeam) printScanbeams() string {
	str := ""
	s2 := s
	for s2 != nil {
		str += " " + fmt.Sprint(s2)
		s2 = s2.Next
	}
	return str
}

type OutRec struct {
	Idx            int
	IsHole, IsOpen bool
	FirstLeft      *OutRec //see comments in clipper.pas
	Pts            *OutPt
	BottomPt       *OutPt
	PolyNode       *PolyNode
}

type OutPt struct {
	Idx        int
	Pt         *IntPoint
	Next, Prev *OutPt
}

func (o *OutPt) String() string {
	return fmt.Sprintf("Idx: %v, Pt: %v",
		o.Idx, o.Pt)
}

type Join struct {
	OutPt1, OutPt2 *OutPt
	OffPt          *IntPoint
}

func (j *Join) String() string {
	return fmt.Sprintf("OutPt1: %v, OutPt2: %v, OffPt: %v",
		j.OutPt1, j.OutPt2, j.OffPt)
}

var horizontal = math.Inf(-1)

const (
	Skip       int     = -2
	Unassigned int     = -1
	tolerance  float64 = 1.0e-20
)

type ClipperBase struct {
	m_MinimaList                   *LocalMinima
	m_CurrentLM                    *LocalMinima
	m_edges                        [][]*TEdge
	m_UseFullRange, m_HasOpenPaths bool
	PreserveCollinear              bool
}

func NewClipperBase() *ClipperBase {
	c := new(ClipperBase)
	c.m_edges = make([][]*TEdge, 0)
	return c
}

func near_zero(val float64) bool {
	return (val > -tolerance) && (val < tolerance)
}

func (c *ClipperBase) Swap(val1, val2 *CInt) {
	*val1, *val2 = *val2, *val1
}

//------------------------------------------------------------------------------

func (c *ClipperBase) IsHorizontal(e *TEdge) bool {
	return e.Delta.Y == 0
}

//------------------------------------------------------------------------------

func (c *ClipperBase) PointIsVertex(pt *IntPoint, pp *OutPt) bool {
	pp2 := pp
	for {
		if pp2.Pt == pt {
			return true
		}
		pp2 = pp2.Next
		if pp2 == pp {
			break
		}
	}
	return false
}

//------------------------------------------------------------------------------

func (c *ClipperBase) PointOnLineSegment(pt,
	linePt1, linePt2 *IntPoint, UseFullRange bool) bool {
	if UseFullRange {
		return ((pt.X == linePt1.X) && (pt.Y == linePt1.Y)) ||
			((pt.X == linePt2.X) && (pt.Y == linePt2.Y)) ||
			(((pt.X > linePt1.X) == (pt.X < linePt2.X)) &&
				((pt.Y > linePt1.Y) == (pt.Y < linePt2.Y)) &&
				Int128Mul((pt.X-linePt1.X), (linePt2.Y-linePt1.Y)).Cmp(
					Int128Mul((linePt2.X-linePt1.X), (pt.Y-linePt1.Y))) == 0)
	} else {
		return ((pt.X == linePt1.X) && (pt.Y == linePt1.Y)) ||
			((pt.X == linePt2.X) && (pt.Y == linePt2.Y)) ||
			(((pt.X > linePt1.X) == (pt.X < linePt2.X)) &&
				((pt.Y > linePt1.Y) == (pt.Y < linePt2.Y)) &&
				((pt.X-linePt1.X)*(linePt2.Y-linePt1.Y) ==
					(linePt2.X-linePt1.X)*(pt.Y-linePt1.Y)))
	}
}

//------------------------------------------------------------------------------

func (c *ClipperBase) PointOnPolygon(pt *IntPoint, pp *OutPt, UseFullRange bool) bool {
	pp2 := pp
	for {
		if c.PointOnLineSegment(pt, pp2.Pt, pp2.Next.Pt, UseFullRange) {
			return true
		}
		pp2 = pp2.Next
		if pp2 == pp {
			break
		}
	}
	return false
}

//------------------------------------------------------------------------------

func (c *ClipperBase) SlopesEqual(e1, e2 *TEdge, UseFullRange bool) bool {
	if UseFullRange {
		return Int128Mul(e1.Delta.Y, e2.Delta.X).Cmp(
			Int128Mul(e1.Delta.X, e2.Delta.Y)) == 0
	} else {
		return (e1.Delta.Y)*(e2.Delta.X) ==
			(e1.Delta.X)*(e2.Delta.Y)
	}
}

//------------------------------------------------------------------------------

func (c *ClipperBase) SlopesEqual3(pt1, pt2,
	pt3 *IntPoint, UseFullRange bool) bool {
	if UseFullRange {
		return Int128Mul(pt1.Y-pt2.Y, pt2.X-pt3.X).Cmp(
			Int128Mul(pt1.X-pt2.X, pt2.Y-pt3.Y)) == 0
	} else {
		return (pt1.Y-pt2.Y)*(pt2.X-pt3.X)-
			(pt1.X-pt2.X)*(pt2.Y-pt3.Y) == 0
	}
}

//------------------------------------------------------------------------------

func (c *ClipperBase) SlopesEqual4(pt1, pt2,
	pt3, pt4 *IntPoint, UseFullRange bool) bool {
	if UseFullRange {
		return Int128Mul(pt1.Y-pt2.Y, pt3.X-pt4.X).Cmp(
			Int128Mul(pt1.X-pt2.X, pt3.Y-pt4.Y)) == 0
	} else {
		return (pt1.Y-pt2.Y)*(pt3.X-pt4.X)-
			(pt1.X-pt2.X)*(pt3.Y-pt4.Y) == 0
	}
}

//------------------------------------------------------------------------------

func (c *ClipperBase) Clear() {
	c.m_edges = make([][]*TEdge, 0)
	c.m_UseFullRange = false
	c.m_HasOpenPaths = false
}

//------------------------------------------------------------------------------

func (c *ClipperBase) DisposeLocalMinimaList() {
	for c.m_MinimaList != nil {
		tmpLm := c.m_MinimaList.Next
		c.m_MinimaList = nil
		c.m_MinimaList = tmpLm
	}
	c.m_CurrentLM = nil
}

//------------------------------------------------------------------------------

func (c *ClipperBase) RangeTest(Pt *IntPoint, useFullRange *bool) {
	if *useFullRange {
		if Pt.X > hiRange || Pt.Y > hiRange || -Pt.X > hiRange || -Pt.Y > hiRange {
			panic(NewClipperException("Coordinate outside allowed range"))
		}
	} else if Pt.X > loRange || Pt.Y > loRange || -Pt.X > loRange || -Pt.Y > loRange {
		*useFullRange = true
		c.RangeTest(Pt, useFullRange)
	}
}

//------------------------------------------------------------------------------

func (c *ClipperBase) InitEdge(e, eNext, ePrev *TEdge, pt *IntPoint) {
	e.Next = eNext
	e.Prev = ePrev
	e.Curr = pt.Copy()
	e.OutIdx = Unassigned
}

//------------------------------------------------------------------------------

func (c *ClipperBase) InitEdge2(e *TEdge, polyType PolyType) {
	if e.Curr.Y >= e.Next.Curr.Y {
		e.Bot = e.Curr
		e.Top = e.Next.Curr
	} else {
		e.Top = e.Curr
		e.Bot = e.Next.Curr
	}
	c.SetDx(e)
	e.PolyTyp = polyType
}

//------------------------------------------------------------------------------

func (c *ClipperBase) FindNextLocMin(E *TEdge) *TEdge {
	var E2 *TEdge
	for {
		for E.Bot != E.Prev.Bot || E.Curr == E.Top {
			E = E.Next
		}
		if E.Dx != horizontal && E.Prev.Dx != horizontal {
			break
		}
		for E.Prev.Dx == horizontal {
			E = E.Prev
		}
		E2 = E
		for E.Dx == horizontal {
			E = E.Next
		}
		if E.Top.Y == E.Prev.Bot.Y {
			continue //ie just an intermediate horz.
		}
		if E2.Prev.Bot.X < E.Bot.X {
			E = E2
		}
		break
	}
	return E
}

//------------------------------------------------------------------------------

func (c *ClipperBase) ProcessBound(E *TEdge, IsClockwise bool) *TEdge {
	EStart := E
	Result := E
	var Horz *TEdge
	var StartX CInt
	if E.Dx == horizontal {
		//first we need to be careful here with open paths because this
		//may not be a true local minima (ie may be following a skip edge).
		//also, watch for adjacent horz edges to start heading left
		//before finishing right ...
		if IsClockwise {
			if E.Prev.Bot.Y == E.Bot.Y {
				StartX = E.Prev.Bot.X
			} else {
				StartX = E.Prev.Top.X
			}
		} else {
			if E.Next.Bot.Y == E.Bot.Y {
				StartX = E.Next.Bot.X
			} else {
				StartX = E.Next.Top.X
			}
		}
		if E.Bot.X != StartX {
			c.ReverseHorizontal(E)
		}
	}
	if Result.OutIdx != Skip {
		if IsClockwise {
			for Result.Top.Y == Result.Next.Bot.Y && Result.Next.OutIdx != Skip {
				Result = Result.Next
			}
			if Result.Dx == horizontal && Result.Next.OutIdx != Skip {
				//nb: at the top of a bound, horizontals are added to the bound
				//only when the preceding edge attaches to the horizontal's left vertex
				//unless a Skip edge is encountered when that becomes the top divide
				Horz = Result
				for Horz.Prev.Dx == horizontal {
					Horz = Horz.Prev
				}
				if Horz.Prev.Top.X == Result.Next.Top.X {
					if !IsClockwise {
						Result = Horz.Prev
					}
				} else if Horz.Prev.Top.X > Result.Next.Top.X {
					Result = Horz.Prev
				}
			}
			for E != Result {
				E.NextInLML = E.Next
				if E.Dx == horizontal && E != EStart && E.Bot.X != E.Prev.Top.X {
					c.ReverseHorizontal(E)
				}
				E = E.Next
			}
			if E.Dx == horizontal && E != EStart && E.Bot.X != E.Prev.Top.X {
				c.ReverseHorizontal(E)
			}
			Result = Result.Next //move to the edge just beyond current bound
		} else {
			for Result.Top.Y == Result.Prev.Bot.Y && Result.Prev.OutIdx != Skip {
				Result = Result.Prev
			}
			if Result.Dx == horizontal && Result.Prev.OutIdx != Skip {
				Horz = Result
				for Horz.Next.Dx == horizontal {
					Horz = Horz.Next
				}
				if Horz.Next.Top.X == Result.Prev.Top.X {
					if !IsClockwise {
						Result = Horz.Next
					}
				} else if Horz.Next.Top.X > Result.Prev.Top.X {
					Result = Horz.Next
				}
			}

			for E != Result {
				E.NextInLML = E.Prev
				if E.Dx == horizontal && E != EStart && E.Bot.X != E.Next.Top.X {
					c.ReverseHorizontal(E)
				}
				E = E.Prev
			}
			if E.Dx == horizontal && E != EStart && E.Bot.X != E.Next.Top.X {
				c.ReverseHorizontal(E)
			}
			Result = Result.Prev //move to the edge just beyond current bound
		}
	}

	if Result.OutIdx == Skip {
		//if edges still remain in the current bound beyond the skip edge then
		//create another LocMin and call ProcessBound once more
		E = Result
		if IsClockwise {
			for E.Top.Y == E.Next.Bot.Y {
				E = E.Next
			}
			//don't include top horizontals when parsing a bound a second time,
			//they will be contained in the opposite bound ...
			for E != Result && E.Dx == horizontal {
				E = E.Prev
			}
		} else {
			for E.Top.Y == E.Prev.Bot.Y {
				E = E.Prev
			}
			for E != Result && E.Dx == horizontal {
				E = E.Next
			}
		}
		if E == Result {
			if IsClockwise {
				Result = E.Next
			} else {
				Result = E.Prev
			}
		} else {
			//there are more edges in the bound beyond result starting with E
			if IsClockwise {
				E = Result.Next
			} else {
				E = Result.Prev
			}
			locMin := new(LocalMinima)
			locMin.Y = E.Bot.Y
			locMin.RightBound = E
			locMin.RightBound.WindDelta = 0
			Result = c.ProcessBound(locMin.RightBound, IsClockwise)
			c.InsertLocalMinima(locMin)
		}
	}
	return Result
}

//------------------------------------------------------------------------------

func (c *ClipperBase) AddPath(pg Path, polyType PolyType, Closed bool) bool {
	if !Closed && polyType == PtClip {
		panic(NewClipperException("AddPath: Open paths must be subject."))
	}

	highI := len(pg) - 1
	if Closed {
		for highI > 0 && (pg[highI] == pg[0]) {
			highI--
		}
	}
	for highI > 0 && (pg[highI] == pg[highI-1]) {
		highI--
	}
	if (Closed && highI < 2) || (!Closed && highI < 1) {
		return false
	}

	//create a new edge array ...
	edges := make([]*TEdge, highI+1)
	for i := 0; i <= highI; i++ {
		edges[i] = new(TEdge)
	}

	IsFlat := true

	//1. Basic (first) edge initialization ...
	edges[1].Curr = pg[1].Copy()
	c.RangeTest(pg[0], &c.m_UseFullRange)
	c.RangeTest(pg[highI], &c.m_UseFullRange)
	c.InitEdge(edges[0], edges[1], edges[highI], pg[0])
	c.InitEdge(edges[highI], edges[0], edges[highI-1], pg[highI])
	for i := highI - 1; i >= 1; i-- {
		c.RangeTest(pg[i], &c.m_UseFullRange)
		c.InitEdge(edges[i], edges[i+1], edges[i-1], pg[i])
	}
	eStart := edges[0]

	//2. Remove duplicate vertices, and (when closed) collinear edges ...
	E := eStart
	eLoopStop := eStart
	for {
		if E.Curr == E.Next.Curr {
			if E == E.Next {
				break
			}
			if E == eStart {
				eStart = E.Next
			}
			E = c.RemoveEdge(E)
			eLoopStop = E
			continue
		}
		if E.Prev == E.Next {
			break //only two vertices
		} else if Closed &&
			c.SlopesEqual3(&E.Prev.Curr, &E.Curr, &E.Next.Curr, c.m_UseFullRange) &&
			(!c.PreserveCollinear ||
				!c.Pt2IsBetweenPt1AndPt3(&E.Prev.Curr, &E.Curr, &E.Next.Curr)) {
			//Collinear edges are allowed for open paths but in closed paths
			//the default is to merge adjacent collinear edges into a single edge.
			//However, if the PreserveCollinear property is enabled, only overlapping
			//collinear edges (ie spikes) will be removed from closed paths.
			if E == eStart {
				eStart = E.Next
			}
			E = c.RemoveEdge(E)
			E = E.Prev
			eLoopStop = E
			continue
		}
		E = E.Next
		if (E == eLoopStop) || (!Closed && E.Next == eStart) {
			break
		}
	}

	if (!Closed && (E == E.Next)) || (Closed && (E.Prev == E.Next)) {
		return false
	}

	if !Closed {
		c.m_HasOpenPaths = true
		eStart.Prev.OutIdx = Skip
	}

	//3. Do second stage of edge initialization ...
	E = eStart
	for {
		c.InitEdge2(E, polyType)
		E = E.Next
		if IsFlat && E.Curr.Y != eStart.Curr.Y {
			IsFlat = false
		}
		if E == eStart {
			break
		}
	}

	//4. Finally, add edge bounds to LocalMinima list ...

	//Totally flat paths must be handled differently when adding them
	//to LocalMinima list to avoid endless loops etc ...
	if IsFlat {
		if Closed {
			return false
		}
		E.Prev.OutIdx = Skip
		if E.Prev.Bot.X < E.Prev.Top.X {
			c.ReverseHorizontal(E.Prev)
		}
		locMin := new(LocalMinima)
		locMin.Next = nil
		locMin.Y = E.Bot.Y
		locMin.LeftBound = nil
		locMin.RightBound = E
		locMin.RightBound.Side = EsRight
		locMin.RightBound.WindDelta = 0
		for E.Next.OutIdx != Skip {
			E.NextInLML = E.Next
			if E.Bot.X != E.Prev.Top.X {
				c.ReverseHorizontal(E)
			}
			E = E.Next
		}
		c.InsertLocalMinima(locMin)
		c.m_edges = append(c.m_edges, edges)
		return true
	}

	c.m_edges = append(c.m_edges, edges)
	var clockwise bool
	var EMin *TEdge
	for {
		E = c.FindNextLocMin(E)
		if E == EMin {
			break
		} else if EMin == nil {
			EMin = E
		}

		//E and E.Prev now share a local minima (left aligned if horizontal).
		//Compare their slopes to find which starts which bound ...
		locMin := new(LocalMinima)
		locMin.Next = nil
		locMin.Y = E.Bot.Y
		if E.Dx < E.Prev.Dx {
			locMin.LeftBound = E.Prev
			locMin.RightBound = E
			clockwise = false //Q.nextInLML = Q.prev
		} else {
			locMin.LeftBound = E
			locMin.RightBound = E.Prev
			clockwise = true //Q.nextInLML = Q.next
		}
		locMin.LeftBound.Side = EsLeft
		locMin.RightBound.Side = EsRight

		if !Closed {
			locMin.LeftBound.WindDelta = 0
		} else if locMin.LeftBound.Next == locMin.RightBound {
			locMin.LeftBound.WindDelta = -1
		} else {
			locMin.LeftBound.WindDelta = 1
		}
		locMin.RightBound.WindDelta = -locMin.LeftBound.WindDelta

		E = c.ProcessBound(locMin.LeftBound, clockwise)
		E2 := c.ProcessBound(locMin.RightBound, !clockwise)

		if locMin.LeftBound.OutIdx == Skip {
			locMin.LeftBound = nil
		} else if locMin.RightBound.OutIdx == Skip {
			locMin.RightBound = nil
		}
		c.InsertLocalMinima(locMin)
		if !clockwise {
			E = E2
		}
	}
	return true
}

//------------------------------------------------------------------------------

func (c *ClipperBase) AddPaths(ppg Paths, polyType PolyType, closed bool) bool {
	result := false
	for i := 0; i < len(ppg); i++ {
		if c.AddPath(ppg[i], polyType, closed) {
			result = true
		}
	}
	return result
}

//------------------------------------------------------------------------------

func (c *ClipperBase) Pt2IsBetweenPt1AndPt3(pt1, pt2, pt3 *IntPoint) bool {
	if (pt1 == pt3) || (pt1 == pt2) || (pt3 == pt2) {
		return false
	} else if pt1.X != pt3.X {
		return (pt2.X > pt1.X) == (pt2.X < pt3.X)
	} else {
		return (pt2.Y > pt1.Y) == (pt2.Y < pt3.Y)
	}
}

//------------------------------------------------------------------------------

func (c *ClipperBase) RemoveEdge(e *TEdge) *TEdge {
	//removes e from float64_linked_list (but without removing from memory)
	e.Prev.Next = e.Next
	e.Next.Prev = e.Prev
	result := e.Next
	e.Prev = nil //flag as removed (see ClipperBase.Clear)
	return result
}

//------------------------------------------------------------------------------

func (c *ClipperBase) SetDx(e *TEdge) {
	e.Delta = IntPoint{(e.Top.X - e.Bot.X), (e.Top.Y - e.Bot.Y)}
	if e.Delta.Y == 0 {
		e.Dx = horizontal
	} else {
		e.Dx = float64(e.Delta.X) / float64(e.Delta.Y)
	}
}

//---------------------------------------------------------------------------

func (c *ClipperBase) InsertLocalMinima(newLm *LocalMinima) {
	if c.m_MinimaList == nil {
		c.m_MinimaList = newLm
	} else if newLm.Y >= c.m_MinimaList.Y {
		newLm.Next = c.m_MinimaList
		c.m_MinimaList = newLm
	} else {
		tmpLm := c.m_MinimaList
		for tmpLm.Next != nil && (newLm.Y < tmpLm.Next.Y) {
			tmpLm = tmpLm.Next
		}
		newLm.Next = tmpLm.Next
		tmpLm.Next = newLm
	}
}

//------------------------------------------------------------------------------

func (c *ClipperBase) PopLocalMinima() {
	if c.m_CurrentLM == nil {
		return
	}
	c.m_CurrentLM = c.m_CurrentLM.Next
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

func (c *ClipperBase) Reset() {
	c.m_CurrentLM = c.m_MinimaList
	if c.m_CurrentLM == nil {
		return //ie nothing to process
	}

	//reset all edges ...
	lm := c.m_MinimaList
	for lm != nil {
		e := lm.LeftBound
		if e != nil {
			e.Curr = e.Bot
			e.Side = EsLeft
			e.OutIdx = Unassigned
		}
		e = lm.RightBound
		if e != nil {
			e.Curr = e.Bot
			e.Side = EsRight
			e.OutIdx = Unassigned
		}
		lm = lm.Next
	}
}

//------------------------------------------------------------------------------

func GetBounds(paths Paths) *IntRect {
	i := 0
	cnt := len(paths)
	for i < cnt && len(paths[i]) == 0 {
		i++
	}
	if i == cnt {
		return &IntRect{0, 0, 0, 0}
	}
	result := new(IntRect)
	result.left = paths[i][0].X
	result.right = result.left
	result.top = paths[i][0].Y
	result.bottom = result.top
	for i < cnt {
		for j := 0; j < len(paths[i]); j++ {
			if paths[i][j].X < result.left {
				result.left = paths[i][j].X
			} else if paths[i][j].X > result.right {
				result.right = paths[i][j].X
			}
			if paths[i][j].Y < result.top {
				result.top = paths[i][j].Y
			} else if paths[i][j].Y > result.bottom {
				result.bottom = paths[i][j].Y
			}
		}
		i++
	}
	return result
}

//InitOptions that can be passed to the constructor ...
type InitOptions int

const (
	IoNone              InitOptions = 0
	IoReverseSolution   InitOptions = 1
	IoStrictlySimple    InitOptions = 2
	IoPreserveCollinear InitOptions = 4
)

type TZFillCallback interface {
	ZFill(bot1, top1, bot2, top2, intersectPt *IntPoint)
}

type Clipper struct {
	ClipperBase

	m_PolyOuts      []*OutRec
	m_ClipType      ClipType
	m_Scanbeam      *Scanbeam
	m_ActiveEdges   *TEdge
	m_SortedEdges   *TEdge
	m_IntersectList []*IntersectNode
	//	m_IntersectNodeComparer         IComparer_IntersectNode
	m_ExecuteLocked                 bool
	m_ClipFillType                  PolyFillType
	m_SubjFillType                  PolyFillType
	m_Joins                         []*Join
	m_GhostJoins                    []*Join
	m_UsingPolyTree                 bool
	ZFillFunction                   TZFillCallback
	ReverseSolution, StrictlySimple bool
}

func NewClipper(initOptions InitOptions) *Clipper {
	c := new(Clipper)
	c.m_edges = make([][]*TEdge, 0)
	c.m_Scanbeam = nil
	c.m_ActiveEdges = nil
	c.m_SortedEdges = nil
	c.m_IntersectList = make([]*IntersectNode, 0)
	//  c.m_IntersectNodeComparer = new MyIntersectNodeSort();
	c.m_ExecuteLocked = false
	c.m_UsingPolyTree = false
	c.m_PolyOuts = make([]*OutRec, 0)
	c.m_Joins = make([]*Join, 0)
	c.m_GhostJoins = make([]*Join, 0)
	c.ReverseSolution = (IoReverseSolution == initOptions)
	c.StrictlySimple = (IoStrictlySimple == initOptions)
	c.PreserveCollinear = (IoPreserveCollinear == initOptions)
	return c
}

//------------------------------------------------------------------------------

func (c *Clipper) DisposeScanbeamList() {
	for c.m_Scanbeam != nil {
		sb2 := c.m_Scanbeam.Next
		c.m_Scanbeam = nil
		c.m_Scanbeam = sb2
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) Reset() {
	c.m_CurrentLM = c.m_MinimaList
	if c.m_CurrentLM != nil {
		//reset all edges ...
		lm := c.m_MinimaList
		for lm != nil {
			e := lm.LeftBound
			if e != nil {
				e.Curr = e.Bot
				e.Side = EsLeft
				e.OutIdx = Unassigned
			}
			e = lm.RightBound
			if e != nil {
				e.Curr = e.Bot
				e.Side = EsRight
				e.OutIdx = Unassigned
			}
			lm = lm.Next
		}
	}
	c.m_Scanbeam = nil
	c.m_ActiveEdges = nil
	c.m_SortedEdges = nil
	lm := c.m_MinimaList
	for lm != nil {
		c.InsertScanbeam(lm.Y)
		lm = lm.Next
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) InsertScanbeam(Y CInt) {
	if c.m_Scanbeam == nil {
		c.m_Scanbeam = new(Scanbeam)
		c.m_Scanbeam.Next = nil
		c.m_Scanbeam.Y = Y
	} else if Y > c.m_Scanbeam.Y {
		newSb := new(Scanbeam)
		newSb.Y = Y
		newSb.Next = c.m_Scanbeam
		c.m_Scanbeam = newSb
	} else {
		sb2 := c.m_Scanbeam
		for sb2.Next != nil && (Y <= sb2.Next.Y) {
			sb2 = sb2.Next
		}
		if Y == sb2.Y {
			return //ie ignores duplicates
		}
		newSb := new(Scanbeam)
		newSb.Y = Y
		newSb.Next = sb2.Next
		sb2.Next = newSb
	}
	// fmt.Printf("slkjfalsjflasjldfkjaslj%v\n", c.m_Scanbeam.printScanbeams())
}

//------------------------------------------------------------------------------

func (c *Clipper) Execute1(clipType ClipType,
	subjFillType, clipFillType PolyFillType) (solution Paths, succeeded bool) {
	if c.m_ExecuteLocked {
		return
	} else {
		defer func() {
			c.DisposeAllPolyPts()
			c.m_ExecuteLocked = false
		}()
	}
	if c.m_HasOpenPaths {
		panic(NewClipperException("Error: PolyTree struct is needed for " +
			"open path clipping."))
	}

	c.m_ExecuteLocked = true
	c.m_SubjFillType = subjFillType
	c.m_ClipFillType = clipFillType
	c.m_ClipType = clipType
	c.m_UsingPolyTree = false

	succeeded = c.ExecuteInternal()
	//build the return polygons ...
	if succeeded {
		solution = c.BuildResult()
	}
	return
}

//------------------------------------------------------------------------------

func (c *Clipper) Execute2(clipType ClipType,
	subjFillType, clipFillType PolyFillType) (polytree *PolyTree, succeeded bool) {
	if c.m_ExecuteLocked {
		return
	} else {
		defer func() {
			c.DisposeAllPolyPts()
			c.m_ExecuteLocked = false
		}()
	}
	c.m_ExecuteLocked = true
	c.m_SubjFillType = subjFillType
	c.m_ClipFillType = clipFillType
	c.m_ClipType = clipType
	c.m_UsingPolyTree = true

	succeeded = c.ExecuteInternal()
	//build the return polygons ...
	if succeeded {
		polytree = NewPolyTree()
		c.BuildResult2(polytree)
	}
	return
}

//------------------------------------------------------------------------------

func (c *Clipper) FixHoleLinkage(outRec *OutRec) {
	//skip if an outermost polygon or
	//already already points to the correct FirstLeft ...
	if outRec.FirstLeft == nil ||
		(outRec.IsHole != outRec.FirstLeft.IsHole &&
			outRec.FirstLeft.Pts != nil) {
		return
	}

	orfl := outRec.FirstLeft
	for orfl != nil && ((orfl.IsHole == outRec.IsHole) || orfl.Pts == nil) {
		orfl = orfl.FirstLeft
	}
	outRec.FirstLeft = orfl
}

//------------------------------------------------------------------------------

func (c *Clipper) ExecuteInternal() bool {
	defer func() {
		c.m_Joins = make([]*Join, 0)
		c.m_GhostJoins = make([]*Join, 0)
	}()

	c.Reset()
	if c.m_CurrentLM == nil {
		return false
	}

	botY := c.PopScanbeam()
	for {
		c.InsertLocalMinimaIntoAEL(botY)
		c.m_GhostJoins = make([]*Join, 0)
		c.ProcessHorizontals(false)
		if c.m_Scanbeam == nil {
			break
		}
		topY := c.PopScanbeam()
		if !c.ProcessIntersections(topY) {
			return false
		}
		c.ProcessEdgesAtTopOfScanbeam(topY)
		botY = topY
		if !(c.m_Scanbeam != nil || c.m_CurrentLM != nil) {
			break
		}
	}

	//fix orientations ...
	for i := 0; i < len(c.m_PolyOuts); i++ {
		outRec := c.m_PolyOuts[i]
		if outRec.Pts == nil || outRec.IsOpen {
			continue
		}
		if (outRec.IsHole != c.ReverseSolution) == (c.area(outRec) > 0) {
			c.ReversePolyPtLinks(outRec.Pts)
		}
	}

	c.JoinCommonEdges()

	for i := 0; i < len(c.m_PolyOuts); i++ {
		outRec := c.m_PolyOuts[i]
		if outRec.Pts != nil && !outRec.IsOpen {
			c.FixupOutPolygon(outRec)
		}
	}

	if c.StrictlySimple {
		c.DoSimplePolygons()
	}

	return true
}

//------------------------------------------------------------------------------

func (c *Clipper) PopScanbeam() CInt {
	Y := c.m_Scanbeam.Y
	c.m_Scanbeam = c.m_Scanbeam.Next
	return Y
}

//------------------------------------------------------------------------------

func (c *Clipper) DisposeAllPolyPts() {
	for i := 0; i < len(c.m_PolyOuts); i++ {
		c.DisposeOutRec(i)
	}
	c.m_PolyOuts = make([]*OutRec, 0)
}

//------------------------------------------------------------------------------

func (c *Clipper) DisposeOutRec(index int) {
	outRec := c.m_PolyOuts[index]
	outRec.Pts = nil
	outRec = nil
	c.m_PolyOuts[index] = nil
}

//------------------------------------------------------------------------------

func (c *Clipper) AddJoin(Op1, Op2 *OutPt, OffPt *IntPoint) {
	j := new(Join)
	j.OutPt1 = Op1
	j.OutPt2 = Op2
	j.OffPt = OffPt
	c.m_Joins = append(c.m_Joins, j)
}

//------------------------------------------------------------------------------

func (c *Clipper) AddGhostJoin(Op *OutPt, OffPt *IntPoint) {
	j := new(Join)
	j.OutPt1 = Op
	j.OffPt = OffPt
	c.m_GhostJoins = append(c.m_GhostJoins, j)
}

//------------------------------------------------------------------------------

func (c *Clipper) InsertLocalMinimaIntoAEL(botY CInt) {
	for c.m_CurrentLM != nil && (c.m_CurrentLM.Y == botY) {
		lb := c.m_CurrentLM.LeftBound
		rb := c.m_CurrentLM.RightBound
		c.PopLocalMinima()

		var Op1 *OutPt
		if lb == nil {
			c.InsertEdgeIntoAEL(rb, nil)
			c.SetWindingCount(rb)
			if c.IsContributing(rb) {
				Op1 = c.AddOutPt(rb, &rb.Bot)
			}
		} else if rb == nil {
			c.InsertEdgeIntoAEL(lb, nil)
			c.SetWindingCount(lb)
			if c.IsContributing(lb) {
				Op1 = c.AddOutPt(lb, &lb.Bot)
			}
			c.InsertScanbeam(lb.Top.Y)
		} else {
			c.InsertEdgeIntoAEL(lb, nil)
			c.InsertEdgeIntoAEL(rb, lb)
			c.SetWindingCount(lb)
			rb.WindCnt = lb.WindCnt
			rb.WindCnt2 = lb.WindCnt2
			if c.IsContributing(lb) {
				Op1 = c.AddLocalMinPoly(lb, rb, &lb.Bot)
			}
			c.InsertScanbeam(lb.Top.Y)
		}

		if rb != nil {
			if c.IsHorizontal(rb) {
				c.AddEdgeToSEL(rb)
			} else {
				c.InsertScanbeam(rb.Top.Y)
			}
		}

		if lb == nil || rb == nil {
			continue
		}

		//if output polygons share an Edge with a horizontal rb, they'll need joining later ...
		if Op1 != nil && c.IsHorizontal(rb) &&
			len(c.m_GhostJoins) > 0 && rb.WindDelta != 0 {
			for i := 0; i < len(c.m_GhostJoins); i++ {
				//if the horizontal Rb and a 'ghost' horizontal overlap, then convert
				//the 'ghost' join to a real join ready for later ...
				j := c.m_GhostJoins[i]
				if c.HorzSegmentsOverlap(j.OutPt1.Pt.X, j.OffPt.X, rb.Bot.X, rb.Top.X) {
					c.AddJoin(j.OutPt1, Op1, j.OffPt)
				}
			}
		}

		if lb.OutIdx >= 0 && lb.PrevInAEL != nil &&
			lb.PrevInAEL.Curr.X == lb.Bot.X &&
			lb.PrevInAEL.OutIdx >= 0 &&
			c.SlopesEqual(lb.PrevInAEL, lb, c.m_UseFullRange) &&
			lb.WindDelta != 0 && lb.PrevInAEL.WindDelta != 0 {
			Op2 := c.AddOutPt(lb.PrevInAEL, &lb.Bot)
			c.AddJoin(Op1, Op2, &lb.Top)
		}

		if lb.NextInAEL != rb {
			if rb.OutIdx >= 0 && rb.PrevInAEL.OutIdx >= 0 &&
				c.SlopesEqual(rb.PrevInAEL, rb, c.m_UseFullRange) &&
				rb.WindDelta != 0 && rb.PrevInAEL.WindDelta != 0 {
				Op2 := c.AddOutPt(rb.PrevInAEL, &rb.Bot)
				c.AddJoin(Op1, Op2, &rb.Top)
			}

			e := lb.NextInAEL
			if e != nil {
				for e != rb {
					//nb: For calculating winding counts etc, IntersectEdges() assumes
					//that param1 will be to the right of param2 ABOVE the intersection ...
					c.IntersectEdges(rb, e, &lb.Curr, true) //order important here
					e = e.NextInAEL
				}
			}
		}
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) InsertEdgeIntoAEL(edge, startEdge *TEdge) {
	if c.m_ActiveEdges == nil {
		edge.PrevInAEL = nil
		edge.NextInAEL = nil
		c.m_ActiveEdges = edge
	} else if startEdge == nil && c.E2InsertsBeforeE1(c.m_ActiveEdges, edge) {
		edge.PrevInAEL = nil
		edge.NextInAEL = c.m_ActiveEdges
		c.m_ActiveEdges.PrevInAEL = edge
		c.m_ActiveEdges = edge
	} else {
		if startEdge == nil {
			startEdge = c.m_ActiveEdges
		}
		for startEdge.NextInAEL != nil &&
			!c.E2InsertsBeforeE1(startEdge.NextInAEL, edge) {
			startEdge = startEdge.NextInAEL
		}
		edge.NextInAEL = startEdge.NextInAEL
		if startEdge.NextInAEL != nil {
			startEdge.NextInAEL.PrevInAEL = edge
		}
		edge.PrevInAEL = startEdge
		startEdge.NextInAEL = edge
	}
}

//----------------------------------------------------------------------

func (c *Clipper) E2InsertsBeforeE1(e1, e2 *TEdge) bool {
	if e2.Curr.X == e1.Curr.X {
		if e2.Top.Y > e1.Top.Y {
			return e2.Top.X < TopX(e1, &e2.Top.Y)
		}
		return e1.Top.X > TopX(e2, &e1.Top.Y)
	}
	return e2.Curr.X < e1.Curr.X
}

//------------------------------------------------------------------------------

func (c *Clipper) IsEvenOddFillType(edge *TEdge) bool {
	if edge.PolyTyp == PtSubject {
		return c.m_SubjFillType == PftEvenOdd
	} else {
		return c.m_ClipFillType == PftEvenOdd
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) IsEvenOddAltFillType(edge *TEdge) bool {
	if edge.PolyTyp == PtSubject {
		return c.m_ClipFillType == PftEvenOdd
	} else {
		return c.m_SubjFillType == PftEvenOdd
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) IsContributing(edge *TEdge) bool {
	var pft, pft2 PolyFillType
	if edge.PolyTyp == PtSubject {
		pft = c.m_SubjFillType
		pft2 = c.m_ClipFillType
	} else {
		pft = c.m_ClipFillType
		pft2 = c.m_SubjFillType
	}

	switch pft {
	case PftEvenOdd:
		//return false if a subj line has been flagged as inside a subj polygon
		if edge.WindDelta == 0 && edge.WindCnt != 1 {
			return false
		}
		break
	case PftNonZero:
		if intAbs(edge.WindCnt) != 1 {
			return false
		}
		break
	case PftPositive:
		if edge.WindCnt != 1 {
			return false
		}
		break
	default: //PolyFillType.PftNegative
		if edge.WindCnt != -1 {
			return false
		}
		break
	}

	switch c.m_ClipType {
	case CtIntersection:
		switch pft2 {
		case PftEvenOdd, PftNonZero:
			return (edge.WindCnt2 != 0)
		case PftPositive:
			return (edge.WindCnt2 > 0)
		default:
			return (edge.WindCnt2 < 0)
		}
	case CtUnion:
		switch pft2 {
		case PftEvenOdd, PftNonZero:
			return (edge.WindCnt2 == 0)
		case PftPositive:
			return (edge.WindCnt2 <= 0)
		default:
			return (edge.WindCnt2 >= 0)
		}
	case CtDifference:
		if edge.PolyTyp == PtSubject {
			switch pft2 {
			case PftEvenOdd, PftNonZero:
				return (edge.WindCnt2 == 0)
			case PftPositive:
				return (edge.WindCnt2 <= 0)
			default:
				return (edge.WindCnt2 >= 0)
			}
		} else {
			switch pft2 {
			case PftEvenOdd, PftNonZero:
				return (edge.WindCnt2 != 0)
			case PftPositive:
				return (edge.WindCnt2 > 0)
			default:
				return (edge.WindCnt2 < 0)
			}
		}
	case CtXor:
		if edge.WindDelta == 0 { //XOr always contributing unless open
			switch pft2 {
			case PftEvenOdd, PftNonZero:
				return (edge.WindCnt2 == 0)
			case PftPositive:
				return (edge.WindCnt2 <= 0)
			default:
				return (edge.WindCnt2 >= 0)
			}
		} else {
			return true
		}
	}
	return true
}

//------------------------------------------------------------------------------

func (c *Clipper) SetWindingCount(edge *TEdge) {
	e := edge.PrevInAEL
	//find the edge of the same polytype that immediately preceeds 'edge' in AEL
	for e != nil && ((e.PolyTyp != edge.PolyTyp) || (e.WindDelta == 0)) {
		e = e.PrevInAEL
	}
	if e == nil {
		if edge.WindDelta == 0 {
			edge.WindCnt = 1
		} else {
			edge.WindCnt = edge.WindDelta
		}
		edge.WindCnt2 = 0
		e = c.m_ActiveEdges //ie get ready to calc WindCnt2
	} else if edge.WindDelta == 0 && c.m_ClipType != CtUnion {
		edge.WindCnt = 1
		edge.WindCnt2 = e.WindCnt2
		e = e.NextInAEL //ie get ready to calc WindCnt2
	} else if c.IsEvenOddFillType(edge) {
		//EvenOdd filling ...
		if edge.WindDelta == 0 {
			//are we inside a subj polygon ...
			Inside := true
			e2 := e.PrevInAEL
			for e2 != nil {
				if e2.PolyTyp == e.PolyTyp && e2.WindDelta != 0 {
					Inside = !Inside
				}
				e2 = e2.PrevInAEL
			}
			if Inside {
				edge.WindCnt = 0
			} else {
				edge.WindCnt = 1
			}
		} else {
			edge.WindCnt = edge.WindDelta
		}
		edge.WindCnt2 = e.WindCnt2
		e = e.NextInAEL //ie get ready to calc WindCnt2
	} else {
		//nonZero, Positive or Negative filling ...
		if e.WindCnt*e.WindDelta < 0 {
			//prev edge is 'decreasing' WindCount (WC) toward zero
			//so we're outside the previous polygon ...
			if intAbs(e.WindCnt) > 1 {
				//outside prev poly but still inside another.
				//when reversing direction of prev poly use the same WC
				if e.WindDelta*edge.WindDelta < 0 {
					edge.WindCnt = e.WindCnt
					//otherwise continue to 'decrease' WC ...
				} else {
					edge.WindCnt = e.WindCnt + edge.WindDelta
				}
			} else {
				//now outside all polys of same polytype so set own WC ...
				if edge.WindDelta == 0 {
					edge.WindCnt = 1
				} else {
					edge.WindCnt = edge.WindDelta
				}
			}
		} else {
			//prev edge is 'increasing' WindCount (WC) away from zero
			//so we're inside the previous polygon ...
			if edge.WindDelta == 0 {
				if e.WindCnt < 0 {
					e.WindCnt = e.WindCnt - 1
				} else {
					e.WindCnt = e.WindCnt + 1
				}
				//if wind direction is reversing prev then use same WC
			} else if e.WindDelta*edge.WindDelta < 0 {
				edge.WindCnt = e.WindCnt
				//otherwise add to WC ...
			} else {
				edge.WindCnt = e.WindCnt + edge.WindDelta
			}
		}
		edge.WindCnt2 = e.WindCnt2
		e = e.NextInAEL //ie get ready to calc WindCnt2
	}

	//update WindCnt2 ...
	if c.IsEvenOddAltFillType(edge) {
		//EvenOdd filling ...
		for e != edge {
			if e.WindDelta != 0 {
				if edge.WindCnt2 == 0 {
					edge.WindCnt2 = 1
				} else {
					edge.WindCnt2 = 0
				}
			}
			e = e.NextInAEL
		}
	} else {
		//nonZero, Positive or Negative filling ...
		for e != edge {
			edge.WindCnt2 += e.WindDelta
			e = e.NextInAEL
		}
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) AddEdgeToSEL(edge *TEdge) {
	//SEL pointers in PEdge are reused to build a list of horizontal edges.
	//However, we don't need to worry about order with horizontal edge processing.
	if c.m_SortedEdges == nil {
		c.m_SortedEdges = edge
		edge.PrevInSEL = nil
		edge.NextInSEL = nil
	} else {
		edge.NextInSEL = c.m_SortedEdges
		edge.PrevInSEL = nil
		c.m_SortedEdges.PrevInSEL = edge
		c.m_SortedEdges = edge
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) CopyAELToSEL() {
	e := c.m_ActiveEdges
	c.m_SortedEdges = e
	for e != nil {
		e.PrevInSEL = e.PrevInAEL
		e.NextInSEL = e.NextInAEL
		e = e.NextInAEL
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) SwapPositionsInAEL(edge1, edge2 *TEdge) {
	//check that one or other edge hasn't already been removed from AEL ...
	if edge1.NextInAEL == edge1.PrevInAEL ||
		edge2.NextInAEL == edge2.PrevInAEL {
		return
	}

	if edge1.NextInAEL == edge2 {
		next := edge2.NextInAEL
		if next != nil {
			next.PrevInAEL = edge1
		}
		prev := edge1.PrevInAEL
		if prev != nil {
			prev.NextInAEL = edge2
		}
		edge2.PrevInAEL = prev
		edge2.NextInAEL = edge1
		edge1.PrevInAEL = edge2
		edge1.NextInAEL = next
	} else if edge2.NextInAEL == edge1 {
		next := edge1.NextInAEL
		if next != nil {
			next.PrevInAEL = edge2
		}
		prev := edge2.PrevInAEL
		if prev != nil {
			prev.NextInAEL = edge1
		}
		edge1.PrevInAEL = prev
		edge1.NextInAEL = edge2
		edge2.PrevInAEL = edge1
		edge2.NextInAEL = next
	} else {
		next := edge1.NextInAEL
		prev := edge1.PrevInAEL
		edge1.NextInAEL = edge2.NextInAEL
		if edge1.NextInAEL != nil {
			edge1.NextInAEL.PrevInAEL = edge1
		}
		edge1.PrevInAEL = edge2.PrevInAEL
		if edge1.PrevInAEL != nil {
			edge1.PrevInAEL.NextInAEL = edge1
		}
		edge2.NextInAEL = next
		if edge2.NextInAEL != nil {
			edge2.NextInAEL.PrevInAEL = edge2
		}
		edge2.PrevInAEL = prev
		if edge2.PrevInAEL != nil {
			edge2.PrevInAEL.NextInAEL = edge2
		}
	}

	if edge1.PrevInAEL == nil {
		c.m_ActiveEdges = edge1
	} else if edge2.PrevInAEL == nil {
		c.m_ActiveEdges = edge2
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) SwapPositionsInSEL(edge1, edge2 *TEdge) {
	if edge1.NextInSEL == nil && edge1.PrevInSEL == nil {
		return
	}
	if edge2.NextInSEL == nil && edge2.PrevInSEL == nil {
		return
	}

	if edge1.NextInSEL == edge2 {
		next := edge2.NextInSEL
		if next != nil {
			next.PrevInSEL = edge1
		}
		prev := edge1.PrevInSEL
		if prev != nil {
			prev.NextInSEL = edge2
		}
		edge2.PrevInSEL = prev
		edge2.NextInSEL = edge1
		edge1.PrevInSEL = edge2
		edge1.NextInSEL = next
	} else if edge2.NextInSEL == edge1 {
		next := edge1.NextInSEL
		if next != nil {
			next.PrevInSEL = edge2
		}
		prev := edge2.PrevInSEL
		if prev != nil {
			prev.NextInSEL = edge1
		}
		edge1.PrevInSEL = prev
		edge1.NextInSEL = edge2
		edge2.PrevInSEL = edge1
		edge2.NextInSEL = next
	} else {
		next := edge1.NextInSEL
		prev := edge1.PrevInSEL
		edge1.NextInSEL = edge2.NextInSEL
		if edge1.NextInSEL != nil {
			edge1.NextInSEL.PrevInSEL = edge1
		}
		edge1.PrevInSEL = edge2.PrevInSEL
		if edge1.PrevInSEL != nil {
			edge1.PrevInSEL.NextInSEL = edge1
		}
		edge2.NextInSEL = next
		if edge2.NextInSEL != nil {
			edge2.NextInSEL.PrevInSEL = edge2
		}
		edge2.PrevInSEL = prev
		if edge2.PrevInSEL != nil {
			edge2.PrevInSEL.NextInSEL = edge2
		}
	}

	if edge1.PrevInSEL == nil {
		c.m_SortedEdges = edge1
	} else if edge2.PrevInSEL == nil {
		c.m_SortedEdges = edge2
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) AddLocalMaxPoly(e1, e2 *TEdge, pt *IntPoint) {
	c.AddOutPt(e1, pt)
	if e2.WindDelta == 0 {
		c.AddOutPt(e2, pt)
	}
	if e1.OutIdx == e2.OutIdx {
		e1.OutIdx = Unassigned
		e2.OutIdx = Unassigned
	} else if e1.OutIdx < e2.OutIdx {
		c.AppendPolygon(e1, e2)
	} else {
		c.AppendPolygon(e2, e1)
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) AddLocalMinPoly(e1, e2 *TEdge, pt *IntPoint) *OutPt {
	var result *OutPt
	var e, prevE *TEdge
	if c.IsHorizontal(e2) || (e1.Dx > e2.Dx) {
		result = c.AddOutPt(e1, pt)
		e2.OutIdx = e1.OutIdx
		e1.Side = EsLeft
		e2.Side = EsRight
		e = e1
		if e.PrevInAEL == e2 {
			prevE = e2.PrevInAEL
		} else {
			prevE = e.PrevInAEL
		}
	} else {
		result = c.AddOutPt(e2, pt)
		e1.OutIdx = e2.OutIdx
		e1.Side = EsRight
		e2.Side = EsLeft
		e = e2
		if e.PrevInAEL == e1 {
			prevE = e1.PrevInAEL
		} else {
			prevE = e.PrevInAEL
		}
	}

	if prevE != nil && prevE.OutIdx >= 0 &&
		(TopX(prevE, &pt.Y) == TopX(e, &pt.Y)) &&
		c.SlopesEqual(e, prevE, c.m_UseFullRange) &&
		(e.WindDelta != 0) && (prevE.WindDelta != 0) {
		outPt := c.AddOutPt(prevE, pt)
		c.AddJoin(result, outPt, &e.Top)
	}
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) CreateOutRec() *OutRec {
	result := new(OutRec)
	result.Idx = Unassigned
	result.IsHole = false
	result.IsOpen = false
	result.FirstLeft = nil
	result.Pts = nil
	result.BottomPt = nil
	result.PolyNode = nil
	c.m_PolyOuts = append(c.m_PolyOuts, result)
	result.Idx = len(c.m_PolyOuts) - 1
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) AddOutPt(e *TEdge, pt *IntPoint) *OutPt {
	ToFront := (e.Side == EsLeft)
	if e.OutIdx < 0 {
		outRec := c.CreateOutRec()
		outRec.IsOpen = (e.WindDelta == 0)
		newOp := new(OutPt)
		outRec.Pts = newOp
		newOp.Idx = outRec.Idx
		p := pt.Copy()
		newOp.Pt = &p
		newOp.Next = newOp
		newOp.Prev = newOp
		if !outRec.IsOpen {
			c.SetHoleState(e, outRec)
		}
		e.OutIdx = outRec.Idx //nb: do this after SetZ !
		return newOp
	} else {
		outRec := c.m_PolyOuts[e.OutIdx]
		//OutRec.Pts is the 'Left-most' point & OutRec.Pts.Prev is the 'Right-most'
		op := outRec.Pts
		if ToFront && *pt == *op.Pt {
			return op
		} else if !ToFront && *pt == *op.Prev.Pt {
			return op.Prev
		}

		newOp := new(OutPt)
		newOp.Idx = outRec.Idx
		p := pt.Copy()
		newOp.Pt = &p
		newOp.Next = op
		newOp.Prev = op.Prev
		newOp.Prev.Next = newOp
		op.Prev = newOp
		if ToFront {
			outRec.Pts = newOp
		}
		return newOp
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) HorzSegmentsOverlap(seg1a, seg1b, seg2a, seg2b CInt) bool {
	if seg1a > seg1b {
		c.Swap(&seg1a, &seg1b)
	}
	if seg2a > seg2b {
		c.Swap(&seg2a, &seg2b)
	}
	return (seg1a < seg2b) && (seg2a < seg1b)
}

//------------------------------------------------------------------------------

func (c *Clipper) SetHoleState(e *TEdge, outRec *OutRec) {
	isHole := false
	e2 := e.PrevInAEL
	for e2 != nil {
		if e2.OutIdx >= 0 && e2.WindDelta != 0 {
			isHole = !isHole
			if outRec.FirstLeft == nil {
				outRec.FirstLeft = c.m_PolyOuts[e2.OutIdx]
			}
		}
		e2 = e2.PrevInAEL
	}
	if isHole {
		outRec.IsHole = true
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) GetDx(pt1, pt2 *IntPoint) float64 {
	if pt1.Y == pt2.Y {
		return horizontal
	} else {
		return float64(pt2.X-pt1.X) / float64(pt2.Y-pt1.Y)
	}
}

//---------------------------------------------------------------------------

func (c *Clipper) FirstIsBottomPt(btmPt1, btmPt2 *OutPt) bool {
	p := btmPt1.Prev
	for (p.Pt == btmPt1.Pt) && (p != btmPt1) {
		p = p.Prev
	}
	dx1p := math.Abs(c.GetDx(btmPt1.Pt, p.Pt))
	p = btmPt1.Next
	for (p.Pt == btmPt1.Pt) && (p != btmPt1) {
		p = p.Next
	}
	dx1n := math.Abs(c.GetDx(btmPt1.Pt, p.Pt))

	p = btmPt2.Prev
	for (p.Pt == btmPt2.Pt) && (p != btmPt2) {
		p = p.Prev
	}
	dx2p := math.Abs(c.GetDx(btmPt2.Pt, p.Pt))
	p = btmPt2.Next
	for (p.Pt == btmPt2.Pt) && (p != btmPt2) {
		p = p.Next
	}
	dx2n := math.Abs(c.GetDx(btmPt2.Pt, p.Pt))
	return (dx1p >= dx2p && dx1p >= dx2n) || (dx1n >= dx2p && dx1n >= dx2n)
}

//------------------------------------------------------------------------------

func (c *Clipper) GetBottomPt(pp *OutPt) *OutPt {
	var dups *OutPt
	p := pp.Next
	for p != pp {
		if p.Pt.Y > pp.Pt.Y {
			pp = p
			dups = nil
		} else if p.Pt.Y == pp.Pt.Y && p.Pt.X <= pp.Pt.X {
			if p.Pt.X < pp.Pt.X {
				dups = nil
				pp = p
			} else {
				if p.Next != pp && p.Prev != pp {
					dups = p
				}
			}
		}
		p = p.Next
	}
	if dups != nil {
		//there appears to be at least 2 vertices at bottomPt so ...
		for dups != p {
			if !c.FirstIsBottomPt(p, dups) {
				pp = dups
			}
			dups = dups.Next
			for dups.Pt != pp.Pt {
				dups = dups.Next
			}
		}
	}
	return pp
}

//------------------------------------------------------------------------------

func (c *Clipper) GetLowermostRec(outRec1, outRec2 *OutRec) *OutRec {
	//work out which polygon fragment has the correct hole state ...
	if outRec1.BottomPt == nil {
		outRec1.BottomPt = c.GetBottomPt(outRec1.Pts)
	}
	if outRec2.BottomPt == nil {
		outRec2.BottomPt = c.GetBottomPt(outRec2.Pts)
	}
	bPt1 := outRec1.BottomPt
	bPt2 := outRec2.BottomPt
	if bPt1.Pt.Y > bPt2.Pt.Y {
		return outRec1
	} else if bPt1.Pt.Y < bPt2.Pt.Y {
		return outRec2
	} else if bPt1.Pt.X < bPt2.Pt.X {
		return outRec1
	} else if bPt1.Pt.X > bPt2.Pt.X {
		return outRec2
	} else if bPt1.Next == bPt1 {
		return outRec2
	} else if bPt2.Next == bPt2 {
		return outRec1
	} else if c.FirstIsBottomPt(bPt1, bPt2) {
		return outRec1
	} else {
		return outRec2
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) Param1RightOfParam2(outRec1, outRec2 *OutRec) bool {
	for {
		outRec1 = outRec1.FirstLeft
		if outRec1 == outRec2 {
			return true
		}
		if outRec1 == nil {
			break
		}
	}
	return false
}

//------------------------------------------------------------------------------

func (c *Clipper) GetOutRec(idx int) *OutRec {
	outrec := c.m_PolyOuts[idx]
	for outrec != c.m_PolyOuts[outrec.Idx] {
		outrec = c.m_PolyOuts[outrec.Idx]
	}
	return outrec
}

//------------------------------------------------------------------------------

func (c *Clipper) AppendPolygon(e1, e2 *TEdge) {
	//get the start and ends of both output polygons ...
	outRec1 := c.m_PolyOuts[e1.OutIdx]
	outRec2 := c.m_PolyOuts[e2.OutIdx]

	var holeStateRec *OutRec
	if c.Param1RightOfParam2(outRec1, outRec2) {
		holeStateRec = outRec2
	} else if c.Param1RightOfParam2(outRec2, outRec1) {
		holeStateRec = outRec1
	} else {
		holeStateRec = c.GetLowermostRec(outRec1, outRec2)
	}

	p1_lft := outRec1.Pts
	p1_rt := p1_lft.Prev
	p2_lft := outRec2.Pts
	p2_rt := p2_lft.Prev

	var side EdgeSide
	//join e2 poly onto e1 poly and delete pointers to e2 ...
	if e1.Side == EsLeft {
		if e2.Side == EsLeft {
			//z y x a b c
			c.ReversePolyPtLinks(p2_lft)
			p2_lft.Next = p1_lft
			p1_lft.Prev = p2_lft
			p1_rt.Next = p2_rt
			p2_rt.Prev = p1_rt
			outRec1.Pts = p2_rt
		} else {
			//x y z a b c
			p2_rt.Next = p1_lft
			p1_lft.Prev = p2_rt
			p2_lft.Prev = p1_rt
			p1_rt.Next = p2_lft
			outRec1.Pts = p2_lft
		}
		side = EsLeft
	} else {
		if e2.Side == EsRight {
			//a b c z y x
			c.ReversePolyPtLinks(p2_lft)
			p1_rt.Next = p2_rt
			p2_rt.Prev = p1_rt
			p2_lft.Next = p1_lft
			p1_lft.Prev = p2_lft
		} else {
			//a b c x y z
			p1_rt.Next = p2_lft
			p2_lft.Prev = p1_rt
			p1_lft.Prev = p2_rt
			p2_rt.Next = p1_lft
		}
		side = EsRight
	}

	outRec1.BottomPt = nil
	if holeStateRec == outRec2 {
		if outRec2.FirstLeft != outRec1 {
			outRec1.FirstLeft = outRec2.FirstLeft
		}
		outRec1.IsHole = outRec2.IsHole
	}
	outRec2.Pts = nil
	outRec2.BottomPt = nil

	outRec2.FirstLeft = outRec1

	OKIdx := e1.OutIdx
	ObsoleteIdx := e2.OutIdx

	e1.OutIdx = Unassigned //nb: safe because we only get here via AddLocalMaxPoly
	e2.OutIdx = Unassigned

	e := c.m_ActiveEdges
	for e != nil {
		if e.OutIdx == ObsoleteIdx {
			e.OutIdx = OKIdx
			e.Side = side
			break
		}
		e = e.NextInAEL
	}
	outRec2.Idx = outRec1.Idx
}

//------------------------------------------------------------------------------

func (c *Clipper) ReversePolyPtLinks(pp *OutPt) {
	if pp == nil {
		return
	}
	var pp1, pp2 *OutPt
	pp1 = pp
	for {
		pp2 = pp1.Next
		pp1.Next = pp1.Prev
		pp1.Prev = pp2
		pp1 = pp2
		if pp1 == pp {
			break
		}
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) SwapSides(edge1, edge2 *TEdge) {
	side := edge1.Side
	edge1.Side = edge2.Side
	edge2.Side = side
}

//------------------------------------------------------------------------------

func (c *Clipper) SwapPolyIndexes(edge1, edge2 *TEdge) {
	outIdx := edge1.OutIdx
	edge1.OutIdx = edge2.OutIdx
	edge2.OutIdx = outIdx
}

//------------------------------------------------------------------------------

// default protect=false
func (c *Clipper) IntersectEdges(e1, e2 *TEdge, pt *IntPoint, protect bool) {
	//e1 will be to the left of e2 BELOW the intersection. Therefore e1 is before
	//e2 in AEL except when e1 is being inserted at the intersection point ...

	e1stops := !protect && e1.NextInLML == nil &&
		e1.Top.X == pt.X && e1.Top.Y == pt.Y
	e2stops := !protect && e2.NextInLML == nil &&
		e2.Top.X == pt.X && e2.Top.Y == pt.Y
	e1Contributing := (e1.OutIdx >= 0)
	e2Contributing := (e2.OutIdx >= 0)

	c.SetZ(pt, e1, e2)

	//#if use_lines
	//if either edge is on an OPEN path ...
	if e1.WindDelta == 0 || e2.WindDelta == 0 {
		//ignore subject-subject open path intersections UNLESS they
		//are both open paths, AND they are both 'contributing maximas' ...
		if e1.WindDelta == 0 && e2.WindDelta == 0 {
			if (e1stops || e2stops) && e1Contributing && e2Contributing {
				c.AddLocalMaxPoly(e1, e2, pt)
			}
			//if intersecting a subj line with a subj poly ...
		} else if e1.PolyTyp == e2.PolyTyp &&
			e1.WindDelta != e2.WindDelta && c.m_ClipType == CtUnion {
			if e1.WindDelta == 0 {
				if e2Contributing {
					c.AddOutPt(e1, pt)
					if e1Contributing {
						e1.OutIdx = Unassigned
					}
				}
			} else {
				if e1Contributing {
					c.AddOutPt(e2, pt)
					if e2Contributing {
						e2.OutIdx = Unassigned
					}
				}
			}
		} else if e1.PolyTyp != e2.PolyTyp {
			if (e1.WindDelta == 0) && intAbs(e2.WindCnt) == 1 &&
				(c.m_ClipType != CtUnion || e2.WindCnt2 == 0) {
				c.AddOutPt(e1, pt)
				if e1Contributing {
					e1.OutIdx = Unassigned
				}
			} else if (e2.WindDelta == 0) && (intAbs(e1.WindCnt) == 1) &&
				(c.m_ClipType != CtUnion || e1.WindCnt2 == 0) {
				c.AddOutPt(e2, pt)
				if e2Contributing {
					e2.OutIdx = Unassigned
				}
			}
		}

		if e1stops {
			if e1.OutIdx < 0 {
				c.DeleteFromAEL(e1)
			} else {
				panic(NewClipperException("Error intersecting polylines"))
			}
		}
		if e2stops {
			if e2.OutIdx < 0 {
				c.DeleteFromAEL(e2)
			} else {
				panic(NewClipperException("Error intersecting polylines"))
			}
		}
		return
	}
	//#endif

	//update winding counts...
	//assumes that e1 will be to the Right of e2 ABOVE the intersection
	if e1.PolyTyp == e2.PolyTyp {
		if c.IsEvenOddFillType(e1) {
			oldE1WindCnt := e1.WindCnt
			e1.WindCnt = e2.WindCnt
			e2.WindCnt = oldE1WindCnt
		} else {
			if e1.WindCnt+e2.WindDelta == 0 {
				e1.WindCnt = -e1.WindCnt
			} else {
				e1.WindCnt += e2.WindDelta
			}
			if e2.WindCnt-e1.WindDelta == 0 {
				e2.WindCnt = -e2.WindCnt
			} else {
				e2.WindCnt -= e1.WindDelta
			}
		}
	} else {
		if !c.IsEvenOddFillType(e2) {
			e1.WindCnt2 += e2.WindDelta
		} else {
			if e1.WindCnt2 == 0 {
				e1.WindCnt2 = 1
			} else {
				e1.WindCnt2 = 0
			}
		}
		if !c.IsEvenOddFillType(e1) {
			e2.WindCnt2 -= e1.WindDelta
		} else {
			if e2.WindCnt2 == 0 {
				e2.WindCnt2 = 1
			} else {
				e2.WindCnt2 = 0
			}
		}
	}

	var e1FillType, e2FillType, e1FillType2, e2FillType2 PolyFillType
	if e1.PolyTyp == PtSubject {
		e1FillType = c.m_SubjFillType
		e1FillType2 = c.m_ClipFillType
	} else {
		e1FillType = c.m_ClipFillType
		e1FillType2 = c.m_SubjFillType
	}
	if e2.PolyTyp == PtSubject {
		e2FillType = c.m_SubjFillType
		e2FillType2 = c.m_ClipFillType
	} else {
		e2FillType = c.m_ClipFillType
		e2FillType2 = c.m_SubjFillType
	}

	var e1Wc, e2Wc int
	switch e1FillType {
	case PftPositive:
		e1Wc = e1.WindCnt
		break
	case PftNegative:
		e1Wc = -e1.WindCnt
		break
	default:
		e1Wc = intAbs(e1.WindCnt)
		break
	}
	switch e2FillType {
	case PftPositive:
		e2Wc = e2.WindCnt
		break
	case PftNegative:
		e2Wc = -e2.WindCnt
		break
	default:
		e2Wc = intAbs(e2.WindCnt)
		break
	}

	if e1Contributing && e2Contributing {
		if e1stops || e2stops ||
			(e1Wc != 0 && e1Wc != 1) || (e2Wc != 0 && e2Wc != 1) ||
			(e1.PolyTyp != e2.PolyTyp && c.m_ClipType != CtXor) {
			c.AddLocalMaxPoly(e1, e2, pt)
		} else {
			c.AddOutPt(e1, pt)
			c.AddOutPt(e2, pt)
			c.SwapSides(e1, e2)
			c.SwapPolyIndexes(e1, e2)
		}
	} else if e1Contributing {
		if e2Wc == 0 || e2Wc == 1 {
			c.AddOutPt(e1, pt)
			c.SwapSides(e1, e2)
			c.SwapPolyIndexes(e1, e2)
		}

	} else if e2Contributing {
		if e1Wc == 0 || e1Wc == 1 {
			c.AddOutPt(e2, pt)
			c.SwapSides(e1, e2)
			c.SwapPolyIndexes(e1, e2)
		}
	} else if (e1Wc == 0 || e1Wc == 1) &&
		(e2Wc == 0 || e2Wc == 1) && !e1stops && !e2stops {
		//neither edge is currently contributing ...
		var e1Wc2, e2Wc2 int
		switch e1FillType2 {
		case PftPositive:
			e1Wc2 = e1.WindCnt2
			break
		case PftNegative:
			e1Wc2 = -e1.WindCnt2
			break
		default:
			e1Wc2 = intAbs(e1.WindCnt2)
			break
		}
		switch e2FillType2 {
		case PftPositive:
			e2Wc2 = e2.WindCnt2
			break
		case PftNegative:
			e2Wc2 = -e2.WindCnt2
			break
		default:
			e2Wc2 = intAbs(e2.WindCnt2)
			break
		}

		if e1.PolyTyp != e2.PolyTyp {
			c.AddLocalMinPoly(e1, e2, pt)
		} else if e1Wc == 1 && e2Wc == 1 {
			switch c.m_ClipType {
			case CtIntersection:
				if e1Wc2 > 0 && e2Wc2 > 0 {
					c.AddLocalMinPoly(e1, e2, pt)
					break
				}
			case CtUnion:
				if e1Wc2 <= 0 && e2Wc2 <= 0 {
					c.AddLocalMinPoly(e1, e2, pt)
					break
				}
			case CtDifference:
				if ((e1.PolyTyp == PtClip) && (e1Wc2 > 0) && (e2Wc2 > 0)) ||
					((e1.PolyTyp == PtSubject) && (e1Wc2 <= 0) && (e2Wc2 <= 0)) {
					c.AddLocalMinPoly(e1, e2, pt)
					break
				}
			case CtXor:
				c.AddLocalMinPoly(e1, e2, pt)
				break
			}
		} else {
			c.SwapSides(e1, e2)
		}
	}

	if (e1stops != e2stops) &&
		((e1stops && (e1.OutIdx >= 0)) || (e2stops && (e2.OutIdx >= 0))) {
		c.SwapSides(e1, e2)
		c.SwapPolyIndexes(e1, e2)
	}

	//finally, delete any non-contributing maxima edges  ...
	if e1stops {
		c.DeleteFromAEL(e1)
	}
	if e2stops {
		c.DeleteFromAEL(e2)
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) DeleteFromAEL(e *TEdge) {
	AelPrev := e.PrevInAEL
	AelNext := e.NextInAEL
	if AelPrev == nil && AelNext == nil && (e != c.m_ActiveEdges) {
		return //already deleted
	}
	if AelPrev != nil {
		AelPrev.NextInAEL = AelNext
	} else {
		c.m_ActiveEdges = AelNext
	}
	if AelNext != nil {
		AelNext.PrevInAEL = AelPrev
	}
	e.NextInAEL = nil
	e.PrevInAEL = nil
}

//------------------------------------------------------------------------------

func (c *Clipper) DeleteFromSEL(e *TEdge) {
	SelPrev := e.PrevInSEL
	SelNext := e.NextInSEL
	if SelPrev == nil && SelNext == nil && (e != c.m_SortedEdges) {
		return //already deleted
	}
	if SelPrev != nil {
		SelPrev.NextInSEL = SelNext
	} else {
		c.m_SortedEdges = SelNext
	}
	if SelNext != nil {
		SelNext.PrevInSEL = SelPrev
	}
	e.NextInSEL = nil
	e.PrevInSEL = nil
}

//------------------------------------------------------------------------------

func (c *Clipper) UpdateEdgeIntoAEL(e *TEdge) *TEdge {
	if e.NextInLML == nil {
		panic(NewClipperException("UpdateEdgeIntoAEL: invalid call"))
	}
	AelPrev := e.PrevInAEL
	AelNext := e.NextInAEL
	e.NextInLML.OutIdx = e.OutIdx
	if AelPrev != nil {
		AelPrev.NextInAEL = e.NextInLML
	} else {
		c.m_ActiveEdges = e.NextInLML
	}
	if AelNext != nil {
		AelNext.PrevInAEL = e.NextInLML
	}
	e.NextInLML.Side = e.Side
	e.NextInLML.WindDelta = e.WindDelta
	e.NextInLML.WindCnt = e.WindCnt
	e.NextInLML.WindCnt2 = e.WindCnt2
	e = e.NextInLML
	e.Curr = e.Bot
	e.PrevInAEL = AelPrev
	e.NextInAEL = AelNext
	if !c.IsHorizontal(e) {
		c.InsertScanbeam(e.Top.Y)
	}
	return e
}

//------------------------------------------------------------------------------

func (c *Clipper) ProcessHorizontals(isTopOfScanbeam bool) {
	horzEdge := c.m_SortedEdges
	for horzEdge != nil {
		c.DeleteFromSEL(horzEdge)
		c.ProcessHorizontal(horzEdge, isTopOfScanbeam)
		horzEdge = c.m_SortedEdges
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) GetHorzDirection(HorzEdge *TEdge, Dir *Direction, Left, Right *CInt) {
	if HorzEdge.Bot.X < HorzEdge.Top.X {
		*Left = HorzEdge.Bot.X
		*Right = HorzEdge.Top.X
		*Dir = DLeftToRight
	} else {
		*Left = HorzEdge.Top.X
		*Right = HorzEdge.Bot.X
		*Dir = DRightToLeft
	}
}

//------------------------------------------------------------------------

func (c *Clipper) ProcessHorizontal(horzEdge *TEdge, isTopOfScanbeam bool) {
	var dir Direction
	var horzLeft, horzRight CInt

	c.GetHorzDirection(horzEdge, &dir, &horzLeft, &horzRight)

	eLastHorz := horzEdge
	var eMaxPair *TEdge
	for eLastHorz.NextInLML != nil && c.IsHorizontal(eLastHorz.NextInLML) {
		eLastHorz = eLastHorz.NextInLML
	}
	if eLastHorz.NextInLML == nil {
		eMaxPair = c.GetMaximaPair(eLastHorz)
	}

	for {
		IsLastHorz := (horzEdge == eLastHorz)
		e := c.GetNextInAEL(horzEdge, dir)
		for e != nil {
			//Break if we've got to the end of an intermediate horizontal edge ...
			//nb: Smaller Dx's are to the right of larger Dx's ABOVE the horizontal.
			if e.Curr.X == horzEdge.Top.X && horzEdge.NextInLML != nil &&
				e.Dx < horzEdge.NextInLML.Dx {
				break
			}

			eNext := c.GetNextInAEL(e, dir) //saves eNext for later

			if (dir == DLeftToRight && e.Curr.X <= horzRight) ||
				(dir == DRightToLeft && e.Curr.X >= horzLeft) {
				//so far we're still in range of the horizontal Edge  but make sure
				//we're at the last of consec. horizontals when matching with eMaxPair
				if e == eMaxPair && IsLastHorz {
					if horzEdge.OutIdx >= 0 {
						op1 := c.AddOutPt(horzEdge, &horzEdge.Top)
						eNextHorz := c.m_SortedEdges
						for eNextHorz != nil {
							if eNextHorz.OutIdx >= 0 &&
								c.HorzSegmentsOverlap(horzEdge.Bot.X,
									horzEdge.Top.X, eNextHorz.Bot.X, eNextHorz.Top.X) {
								op2 := c.AddOutPt(eNextHorz, &eNextHorz.Bot)
								c.AddJoin(op2, op1, &eNextHorz.Top)
							}
							eNextHorz = eNextHorz.NextInSEL
						}
						c.AddGhostJoin(op1, &horzEdge.Bot)
						c.AddLocalMaxPoly(horzEdge, eMaxPair, &horzEdge.Top)
					}
					c.DeleteFromAEL(horzEdge)
					c.DeleteFromAEL(eMaxPair)
					return
				} else if dir == DLeftToRight {
					Pt := &IntPoint{e.Curr.X, horzEdge.Curr.Y}
					c.IntersectEdges(horzEdge, e, Pt, true)
				} else {
					Pt := &IntPoint{e.Curr.X, horzEdge.Curr.Y}
					c.IntersectEdges(e, horzEdge, Pt, true)
				}
				c.SwapPositionsInAEL(horzEdge, e)
			} else if (dir == DLeftToRight && e.Curr.X >= horzRight) ||
				(dir == DRightToLeft && e.Curr.X <= horzLeft) {
				break
			}
			e = eNext
		}

		if horzEdge.NextInLML != nil && c.IsHorizontal(horzEdge.NextInLML) {
			horzEdge = c.UpdateEdgeIntoAEL(horzEdge)
			if horzEdge.OutIdx >= 0 {
				c.AddOutPt(horzEdge, &horzEdge.Bot)
			}
			c.GetHorzDirection(horzEdge, &dir, &horzLeft, &horzRight)
		} else {
			break
		}
	}

	if horzEdge.NextInLML != nil {
		if horzEdge.OutIdx >= 0 {
			op1 := c.AddOutPt(horzEdge, &horzEdge.Top)
			if isTopOfScanbeam {
				c.AddGhostJoin(op1, &horzEdge.Bot)
			}

			horzEdge = c.UpdateEdgeIntoAEL(horzEdge)
			if horzEdge.WindDelta == 0 {
				return
			}
			//nb: HorzEdge is no longer horizontal here
			ePrev := horzEdge.PrevInAEL
			eNext := horzEdge.NextInAEL
			if ePrev != nil && ePrev.Curr.X == horzEdge.Bot.X &&
				ePrev.Curr.Y == horzEdge.Bot.Y && ePrev.WindDelta != 0 &&
				(ePrev.OutIdx >= 0 && ePrev.Curr.Y > ePrev.Top.Y &&
					c.SlopesEqual(horzEdge, ePrev, c.m_UseFullRange)) {
				op2 := c.AddOutPt(ePrev, &horzEdge.Bot)
				c.AddJoin(op1, op2, &horzEdge.Top)
			} else if eNext != nil && eNext.Curr.X == horzEdge.Bot.X &&
				eNext.Curr.Y == horzEdge.Bot.Y && eNext.WindDelta != 0 &&
				eNext.OutIdx >= 0 && eNext.Curr.Y > eNext.Top.Y &&
				c.SlopesEqual(horzEdge, eNext, c.m_UseFullRange) {
				op2 := c.AddOutPt(eNext, &horzEdge.Bot)
				c.AddJoin(op1, op2, &horzEdge.Top)
			}
		} else {
			horzEdge = c.UpdateEdgeIntoAEL(horzEdge)
		}
	} else {
		if horzEdge.OutIdx >= 0 {
			c.AddOutPt(horzEdge, &horzEdge.Top)
		}
		c.DeleteFromAEL(horzEdge)
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) GetNextInAEL(e *TEdge, direction Direction) *TEdge {
	if direction == DLeftToRight {
		return e.NextInAEL
	} else {
		return e.PrevInAEL
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) IsMinima(e *TEdge) bool {
	return e != nil && (e.Prev.NextInLML != e) && (e.Next.NextInLML != e)
}

//------------------------------------------------------------------------------

func (c *Clipper) IsMaxima(e *TEdge, Y CInt) bool {
	return (e != nil && e.Top.Y == Y && e.NextInLML == nil)
}

//------------------------------------------------------------------------------

func (c *Clipper) IsIntermediate(e *TEdge, Y CInt) bool {
	return (e.Top.Y == Y && e.NextInLML != nil)
}

//------------------------------------------------------------------------------

func (c *Clipper) GetMaximaPair(e *TEdge) *TEdge {
	var result *TEdge
	if e.Next.Top == e.Top && e.Next.NextInLML == nil {
		result = e.Next
	} else if e.Prev.Top == e.Top && e.Prev.NextInLML == nil {
		result = e.Prev
	}
	if result != nil && (result.OutIdx == Skip ||
		(result.NextInAEL == result.PrevInAEL && !c.IsHorizontal(result))) {
		return nil
	}
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) ProcessIntersections(topY CInt) bool {
	if c.m_ActiveEdges == nil {
		return true
	}
	//	defer func() {
	//		if r := recover(); r != nil {
	//			c.m_SortedEdges = nil
	//			c.m_IntersectList = make([]*IntersectNode, 0)
	//			panic(NewClipperException("ProcessIntersections error: " +
	//				r.(error).Error()))
	//		}
	//	}()
	c.BuildIntersectList(topY)
	if len(c.m_IntersectList) == 0 {
		return true
	}
	if len(c.m_IntersectList) == 1 || c.FixupIntersectionOrder() {
		c.ProcessIntersectList()
	} else {
		return false
	}

	c.m_SortedEdges = nil
	return true
}

//------------------------------------------------------------------------------

func (c *Clipper) BuildIntersectList(topY CInt) {
	if c.m_ActiveEdges == nil {
		return
	}

	//prepare for sorting ...
	e := c.m_ActiveEdges
	c.m_SortedEdges = e
	for e != nil {
		e.PrevInSEL = e.PrevInAEL
		e.NextInSEL = e.NextInAEL
		e.Curr.X = TopX(e, &topY)
		e = e.NextInAEL
	}

	//bubblesort ...
	isModified := true
	for isModified && c.m_SortedEdges != nil {
		isModified = false
		e = c.m_SortedEdges
		for e.NextInSEL != nil {
			eNext := e.NextInSEL

			if e.Curr.X > eNext.Curr.X {
				newNode := new(IntersectNode)
				newNode.Edge1 = e
				newNode.Edge2 = eNext
				newNode.Pt = c.IntersectPoint(e, eNext)
				c.m_IntersectList = append(c.m_IntersectList, newNode)

				c.SwapPositionsInSEL(e, eNext)
				isModified = true
			} else {
				e = eNext
			}
		}
		if e.PrevInSEL != nil {
			e.PrevInSEL.NextInSEL = nil
		} else {
			break
		}
	}
	// fmt.Println("IntersectList", c.m_IntersectList)
	c.m_SortedEdges = nil
}

//------------------------------------------------------------------------------

func (c *Clipper) EdgesAdjacent(inode *IntersectNode) bool {
	return (inode.Edge1.NextInSEL == inode.Edge2) ||
		(inode.Edge1.PrevInSEL == inode.Edge2)
}

//------------------------------------------------------------------------------

func (c *Clipper) IntersectNodeSort(node1, node2 *IntersectNode) int {
	//the following typecast is safe because the differences in Pt.Y will
	//be limited to the height of the scanbeam.
	return int(node2.Pt.Y - node1.Pt.Y)
}

//------------------------------------------------------------------------------

func (c *Clipper) FixupIntersectionOrder() bool {
	//pre-condition: intersections are sorted bottom-most first.
	//Now it's crucial that intersections are made only between adjacent edges,
	//so to ensure this the order of intersections may need adjusting ...
	sort.Sort(IntersectNodeList(c.m_IntersectList))

	c.CopyAELToSEL()
	cnt := len(c.m_IntersectList)
	for i := 0; i < cnt; i++ {
		if !c.EdgesAdjacent(c.m_IntersectList[i]) {
			j := i + 1
			for j < cnt && !c.EdgesAdjacent(c.m_IntersectList[j]) {
				j++
			}
			if j == cnt {
				return false
			}

			tmp := c.m_IntersectList[i]
			c.m_IntersectList[i] = c.m_IntersectList[j]
			c.m_IntersectList[j] = tmp

		}
		c.SwapPositionsInSEL(c.m_IntersectList[i].Edge1, c.m_IntersectList[i].Edge2)
	}
	return true
}

//------------------------------------------------------------------------------

func (c *Clipper) ProcessIntersectList() {
	for i := 0; i < len(c.m_IntersectList); i++ {
		iNode := c.m_IntersectList[i]
		{
			c.IntersectEdges(iNode.Edge1, iNode.Edge2, iNode.Pt, true)
			c.SwapPositionsInAEL(iNode.Edge1, iNode.Edge2)
		}
	}
	c.m_IntersectList = make([]*IntersectNode, 0)
}

//------------------------------------------------------------------------------

func Round(value float64) CInt {
	if value < 0 {
		return CInt(value - 0.5)
	} else {
		return CInt(value + 0.5)
	}
}

//------------------------------------------------------------------------------

func TopX(edge *TEdge, currentY *CInt) CInt {
	if *currentY == edge.Top.Y {
		return edge.Top.X
	}
	// fmt.Println("TopX", edge, "currentY", *currentY)
	return edge.Bot.X + Round(edge.Dx*float64(*currentY-edge.Bot.Y))
}

//------------------------------------------------------------------------------

func (c *Clipper) IntersectPoint(edge1, edge2 *TEdge) (ip *IntPoint) {
	ip = new(IntPoint)
	var b1, b2 float64
	//nb: with very large coordinate values, it's possible for SlopesEqual() to
	//return false but for the edge.Dx value be equal due to float64 precision rounding.
	if edge1.Dx == edge2.Dx {
		ip.Y = edge1.Curr.Y
		ip.X = TopX(edge1, &ip.Y)
		return
	}

	if edge1.Delta.X == 0 {
		ip.X = edge1.Bot.X
		if c.IsHorizontal(edge2) {
			ip.Y = edge2.Bot.Y
		} else {
			b2 = float64(edge2.Bot.Y) - (float64(edge2.Bot.X) / edge2.Dx)
			ip.Y = Round(float64(ip.X)/edge2.Dx + b2)
		}
	} else if edge2.Delta.X == 0 {
		ip.X = edge2.Bot.X
		if c.IsHorizontal(edge1) {
			ip.Y = edge1.Bot.Y
		} else {
			b1 = float64(edge1.Bot.Y) - (float64(edge1.Bot.X) / edge1.Dx)
			ip.Y = Round(float64(ip.X)/edge1.Dx + b1)
		}
	} else {
		b1 = float64(edge1.Bot.X) - float64(edge1.Bot.Y)*edge1.Dx
		b2 = float64(edge2.Bot.X) - float64(edge2.Bot.Y)*edge2.Dx
		q := (b2 - b1) / (edge1.Dx - edge2.Dx)
		ip.Y = Round(q)
		if math.Abs(edge1.Dx) < math.Abs(edge2.Dx) {
			ip.X = Round(edge1.Dx*q + b1)
		} else {
			ip.X = Round(edge2.Dx*q + b2)
		}
	}

	if ip.Y < edge1.Top.Y || ip.Y < edge2.Top.Y {
		if edge1.Top.Y > edge2.Top.Y {
			ip.Y = edge1.Top.Y
		} else {
			ip.Y = edge2.Top.Y
		}
		if math.Abs(edge1.Dx) < math.Abs(edge2.Dx) {
			ip.X = TopX(edge1, &ip.Y)
		} else {
			ip.X = TopX(edge2, &ip.Y)
		}
	}
	//finally, don't allow 'ip' to be BELOW curr.Y (ie bottom of scanbeam) ...
	if ip.Y > edge1.Curr.Y {
		ip.Y = edge1.Curr.Y
		//better to use the more vertical edge to derive X ...
		if math.Abs(edge1.Dx) > math.Abs(edge2.Dx) {
			ip.X = TopX(edge2, &ip.Y)
		} else {
			ip.X = TopX(edge1, &ip.Y)
		}
	}
	return
}

//------------------------------------------------------------------------------

func (c *Clipper) ProcessEdgesAtTopOfScanbeam(topY CInt) {
	e := c.m_ActiveEdges
	for e != nil {
		//1. process maxima, treating them as if they're 'bent' horizontal edges,
		//   but exclude maxima with horizontal edges. nb: e can't be a horizontal.
		IsMaximaEdge := c.IsMaxima(e, topY)

		if IsMaximaEdge {
			eMaxPair := c.GetMaximaPair(e)
			IsMaximaEdge = (eMaxPair == nil || !c.IsHorizontal(eMaxPair))
		}

		if IsMaximaEdge {
			ePrev := e.PrevInAEL
			c.DoMaxima(e)
			if ePrev == nil {
				e = c.m_ActiveEdges
			} else {
				e = ePrev.NextInAEL
			}
		} else {
			//2. promote horizontal edges, otherwise update Curr.X and Curr.Y ...
			if c.IsIntermediate(e, topY) && c.IsHorizontal(e.NextInLML) {
				e = c.UpdateEdgeIntoAEL(e)
				if e.OutIdx >= 0 {
					c.AddOutPt(e, &e.Bot)
				}
				c.AddEdgeToSEL(e)
			} else {
				e.Curr.X = TopX(e, &topY)
				e.Curr.Y = topY
			}

			if c.StrictlySimple {
				ePrev := e.PrevInAEL
				if (e.OutIdx >= 0) && (e.WindDelta != 0) && ePrev != nil &&
					(ePrev.OutIdx >= 0) && (ePrev.Curr.X == e.Curr.X) &&
					(ePrev.WindDelta != 0) {
					ip := e.Curr.Copy()
					c.SetZ(&ip, ePrev, e)
					op := c.AddOutPt(ePrev, &ip)
					op2 := c.AddOutPt(e, &ip)
					c.AddJoin(op, op2, &ip) //StrictlySimple (type-3) join
				}
			}

			e = e.NextInAEL
		}
	}

	//3. Process horizontals at the Top of the scanbeam ...
	c.ProcessHorizontals(true)

	//4. Promote intermediate vertices ...
	e = c.m_ActiveEdges
	for e != nil {
		if c.IsIntermediate(e, topY) {
			var op *OutPt
			if e.OutIdx >= 0 {
				op = c.AddOutPt(e, &e.Top)
			}
			e = c.UpdateEdgeIntoAEL(e)

			//if output polygons share an edge, they'll need joining later ...
			ePrev := e.PrevInAEL
			eNext := e.NextInAEL
			if ePrev != nil && ePrev.Curr.X == e.Bot.X &&
				ePrev.Curr.Y == e.Bot.Y && op != nil &&
				ePrev.OutIdx >= 0 && ePrev.Curr.Y > ePrev.Top.Y &&
				c.SlopesEqual(e, ePrev, c.m_UseFullRange) &&
				(e.WindDelta != 0) && (ePrev.WindDelta != 0) {
				op2 := c.AddOutPt(ePrev, &e.Bot)
				c.AddJoin(op, op2, &e.Top)
			} else if eNext != nil && eNext.Curr.X == e.Bot.X &&
				eNext.Curr.Y == e.Bot.Y && op != nil &&
				eNext.OutIdx >= 0 && eNext.Curr.Y > eNext.Top.Y &&
				c.SlopesEqual(e, eNext, c.m_UseFullRange) &&
				(e.WindDelta != 0) && (eNext.WindDelta != 0) {
				op2 := c.AddOutPt(eNext, &e.Bot)
				c.AddJoin(op, op2, &e.Top)
			}
		}
		e = e.NextInAEL
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) DoMaxima(e *TEdge) {
	eMaxPair := c.GetMaximaPair(e)
	if eMaxPair == nil {
		if e.OutIdx >= 0 {
			c.AddOutPt(e, &e.Top)
		}
		c.DeleteFromAEL(e)
		return
	}

	eNext := e.NextInAEL
	for eNext != nil && eNext != eMaxPair {
		c.IntersectEdges(e, eNext, &e.Top, true)
		c.SwapPositionsInAEL(e, eNext)
		eNext = e.NextInAEL
	}

	if e.OutIdx == Unassigned && eMaxPair.OutIdx == Unassigned {
		c.DeleteFromAEL(e)
		c.DeleteFromAEL(eMaxPair)
	} else if e.OutIdx >= 0 && eMaxPair.OutIdx >= 0 {
		if e.OutIdx >= 0 {
			c.AddLocalMaxPoly(e, eMaxPair, &e.Top)
		}
		c.DeleteFromAEL(e)
		c.DeleteFromAEL(eMaxPair)
		//#if use_lines
	} else if e.WindDelta == 0 {
		if e.OutIdx >= 0 {
			c.AddOutPt(e, &e.Top)
			e.OutIdx = Unassigned
		}
		c.DeleteFromAEL(e)

		if eMaxPair.OutIdx >= 0 {
			c.AddOutPt(eMaxPair, &e.Top)
			eMaxPair.OutIdx = Unassigned
		}
		c.DeleteFromAEL(eMaxPair)
		//#endif
	} else {
		panic(NewClipperException("DoMaxima error"))
	}
}

//------------------------------------------------------------------------------

func reversePath(s Path) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}

func (c *Clipper) ReversePaths(polys Paths) {
	for _, poly := range polys {
		reversePath(poly)
	}
}

//------------------------------------------------------------------------------

func Orientation(poly Path) bool {
	return Area(poly) >= 0
}

//------------------------------------------------------------------------------

func (c *Clipper) PointCount(pts *OutPt) int {
	if pts == nil {
		return 0
	}
	result := 0
	p := pts
	for {
		result++
		p = p.Next
		if p == pts {
			break
		}
	}
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) BuildResult() Paths {
	polyg := Paths(make([]Path, 0, len(c.m_PolyOuts)))
	for i := 0; i < len(c.m_PolyOuts); i++ {
		outRec := c.m_PolyOuts[i]
		if outRec.Pts == nil {
			continue
		}
		p := outRec.Pts.Prev
		cnt := c.PointCount(p)
		if cnt < 2 {
			continue
		}
		pg := Path(make([]*IntPoint, 0, cnt))
		for j := 0; j < cnt; j++ {
			pg = append(pg, p.Pt)
			p = p.Prev
		}
		polyg = append(polyg, pg)
	}
	return polyg
}

//------------------------------------------------------------------------------

func (c *Clipper) BuildResult2(polytree *PolyTree) {
	polytree.Clear()

	//add each output polygon/contour to polytree ...
	polytree.m_AllPolys = make([]*PolyNode, len(c.m_PolyOuts))
	for i := 0; i < len(c.m_PolyOuts); i++ {
		outRec := c.m_PolyOuts[i]
		cnt := c.PointCount(outRec.Pts)
		if (outRec.IsOpen && cnt < 2) ||
			(!outRec.IsOpen && cnt < 3) {
			continue
		}
		c.FixHoleLinkage(outRec)
		pn := new(PolyNode)
		polytree.m_AllPolys[i] = pn
		outRec.PolyNode = pn
		pn.m_polygon = make([]*IntPoint, cnt)
		op := outRec.Pts.Prev
		for j := 0; j < cnt; j++ {
			pn.m_polygon[j] = op.Pt
			op = op.Prev
		}
	}

	//fixup PolyNode links etc ...
	polytree.m_Childs = make([]*PolyNode, 0, len(c.m_PolyOuts))
	for i := 0; i < len(c.m_PolyOuts); i++ {
		outRec := c.m_PolyOuts[i]
		if outRec.PolyNode == nil {
			continue
		} else if outRec.IsOpen {
			outRec.PolyNode.IsOpen = true
			polytree.AddChild(outRec.PolyNode)
		} else if outRec.FirstLeft != nil &&
			outRec.FirstLeft.PolyNode != nil {
			outRec.FirstLeft.PolyNode.AddChild(outRec.PolyNode)
		} else {
			polytree.AddChild(outRec.PolyNode)
		}
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) FixupOutPolygon(outRec *OutRec) {
	//FixupOutPolygon() - removes duplicate points and simplifies consecutive
	//parallel edges by removing the middle vertex.
	var lastOK *OutPt
	outRec.BottomPt = nil
	pp := outRec.Pts
	for {
		if pp.Prev == pp || pp.Prev == pp.Next {
			outRec.Pts = nil
			return
		}
		//test for duplicate points and collinear edges ...
		if (pp.Pt == pp.Next.Pt) || (pp.Pt == pp.Prev.Pt) ||
			(c.SlopesEqual3(pp.Prev.Pt, pp.Pt, pp.Next.Pt, c.m_UseFullRange) &&
				(!c.PreserveCollinear || !c.Pt2IsBetweenPt1AndPt3(pp.Prev.Pt, pp.Pt, pp.Next.Pt))) {
			lastOK = nil
			pp.Prev.Next = pp.Next
			pp.Next.Prev = pp.Prev
			pp = pp.Prev
		} else if pp == lastOK {
			break
		} else {
			if lastOK == nil {
				lastOK = pp
			}
			pp = pp.Next
		}
	}
	outRec.Pts = pp
}

//------------------------------------------------------------------------------

func (c *Clipper) DupOutPt(outPt *OutPt, InsertAfter bool) *OutPt {
	result := new(OutPt)
	result.Pt = outPt.Pt
	result.Idx = outPt.Idx
	if InsertAfter {
		result.Next = outPt.Next
		result.Prev = outPt
		outPt.Next.Prev = result
		outPt.Next = result
	} else {
		result.Prev = outPt.Prev
		result.Next = outPt
		outPt.Prev.Next = result
		outPt.Prev = result
	}
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) GetOverlap(a1, a2, b1, b2 CInt, Left, Right *CInt) bool {
	if a1 < a2 {
		if b1 < b2 {
			*Left = max(a1, b1)
			*Right = min(a2, b2)
		} else {
			*Left = max(a1, b2)
			*Right = min(a2, b1)
		}
	} else {
		if b1 < b2 {
			*Left = max(a2, b1)
			*Right = min(a1, b2)
		} else {
			*Left = max(a2, b2)
			*Right = min(a1, b1)
		}
	}
	return *Left < *Right
}

//------------------------------------------------------------------------------

func (c *Clipper) JoinHorz(op1, op1b, op2, op2b *OutPt,
	Pt *IntPoint, DiscardLeft bool) bool {
	var Dir1 Direction
	if op1.Pt.X > op1b.Pt.X {
		Dir1 = DRightToLeft
	} else {
		Dir1 = DLeftToRight
	}
	var Dir2 Direction
	if op2.Pt.X > op2b.Pt.X {
		Dir2 = DRightToLeft
	} else {
		Dir2 = DLeftToRight
	}
	if Dir1 == Dir2 {
		return false
	}

	//When DiscardLeft, we want Op1b to be on the Left of Op1, otherwise we
	//want Op1b to be on the Right. (And likewise with Op2 and Op2b.)
	//So, to facilitate this while inserting Op1b and Op2b ...
	//when DiscardLeft, make sure we're AT or RIGHT of Pt before adding Op1b,
	//otherwise make sure we're AT or LEFT of Pt. (Likewise with Op2b.)
	if Dir1 == DLeftToRight {
		for op1.Next.Pt.X <= Pt.X &&
			op1.Next.Pt.X >= op1.Pt.X && op1.Next.Pt.Y == Pt.Y {
			op1 = op1.Next
		}
		if DiscardLeft && (op1.Pt.X != Pt.X) {
			op1 = op1.Next
		}
		op1b = c.DupOutPt(op1, !DiscardLeft)
		if op1b.Pt != Pt {
			op1 = op1b
			op1.Pt = Pt
			op1b = c.DupOutPt(op1, !DiscardLeft)
		}
	} else {
		for op1.Next.Pt.X >= Pt.X &&
			op1.Next.Pt.X <= op1.Pt.X && op1.Next.Pt.Y == Pt.Y {
			op1 = op1.Next
		}
		if !DiscardLeft && (op1.Pt.X != Pt.X) {
			op1 = op1.Next
		}
		op1b = c.DupOutPt(op1, DiscardLeft)
		if op1b.Pt != Pt {
			op1 = op1b
			op1.Pt = Pt
			op1b = c.DupOutPt(op1, DiscardLeft)
		}
	}

	if Dir2 == DLeftToRight {
		for op2.Next.Pt.X <= Pt.X &&
			op2.Next.Pt.X >= op2.Pt.X && op2.Next.Pt.Y == Pt.Y {
			op2 = op2.Next
		}
		if DiscardLeft && (op2.Pt.X != Pt.X) {
			op2 = op2.Next
		}
		op2b = c.DupOutPt(op2, !DiscardLeft)
		if op2b.Pt != Pt {
			op2 = op2b
			op2.Pt = Pt
			op2b = c.DupOutPt(op2, !DiscardLeft)
		}
	} else {
		for op2.Next.Pt.X >= Pt.X &&
			op2.Next.Pt.X <= op2.Pt.X && op2.Next.Pt.Y == Pt.Y {
			op2 = op2.Next
		}
		if !DiscardLeft && (op2.Pt.X != Pt.X) {
			op2 = op2.Next
		}
		op2b = c.DupOutPt(op2, DiscardLeft)
		if op2b.Pt != Pt {
			op2 = op2b
			op2.Pt = Pt
			op2b = c.DupOutPt(op2, DiscardLeft)
		}
	}

	if (Dir1 == DLeftToRight) == DiscardLeft {
		op1.Prev = op2
		op2.Next = op1
		op1b.Next = op2b
		op2b.Prev = op1b
	} else {
		op1.Next = op2
		op2.Prev = op1
		op1b.Prev = op2b
		op2b.Next = op1b
	}
	return true
}

//------------------------------------------------------------------------------

func (c *Clipper) JoinPoints(j *Join, outRec1, outRec2 *OutRec) bool {
	op1 := j.OutPt1
	op2 := j.OutPt2
	var op1b, op2b *OutPt

	//There are 3 kinds of joins for output polygons ...
	//1. Horizontal joins where Join.OutPt1 & Join.OutPt2 are a vertices anywhere
	//along (horizontal) collinear edges (& Join.OffPt is on the same horizontal).
	//2. Non-horizontal joins where Join.OutPt1 & Join.OutPt2 are at the same
	//location at the Bottom of the overlapping segment (& Join.OffPt is above).
	//3. StrictlySimple joins where edges touch but are not collinear and where
	//Join.OutPt1, Join.OutPt2 & Join.OffPt all share the same point.
	isHorizontal := (j.OutPt1.Pt.Y == j.OffPt.Y)

	if isHorizontal && (j.OffPt == j.OutPt1.Pt) && (j.OffPt == j.OutPt2.Pt) {
		//Strictly Simple join ...
		if outRec1 != outRec2 {
			return false
		}
		op1b = j.OutPt1.Next
		for op1b != op1 && (op1b.Pt == j.OffPt) {
			op1b = op1b.Next
		}
		reverse1 := (op1b.Pt.Y > j.OffPt.Y)
		op2b = j.OutPt2.Next
		for op2b != op2 && (op2b.Pt == j.OffPt) {
			op2b = op2b.Next
		}
		reverse2 := (op2b.Pt.Y > j.OffPt.Y)
		if reverse1 == reverse2 {
			return false
		}
		if reverse1 {
			op1b = c.DupOutPt(op1, false)
			op2b = c.DupOutPt(op2, true)
			op1.Prev = op2
			op2.Next = op1
			op1b.Next = op2b
			op2b.Prev = op1b
			j.OutPt1 = op1
			j.OutPt2 = op1b
			return true
		} else {
			op1b = c.DupOutPt(op1, true)
			op2b = c.DupOutPt(op2, false)
			op1.Next = op2
			op2.Prev = op1
			op1b.Prev = op2b
			op2b.Next = op1b
			j.OutPt1 = op1
			j.OutPt2 = op1b
			return true
		}
	} else if isHorizontal {
		//treat horizontal joins differently to non-horizontal joins since with
		//them we're not yet sure where the overlapping is. OutPt1.Pt & OutPt2.Pt
		//may be anywhere along the horizontal edge.
		op1b = op1
		for op1.Prev.Pt.Y == op1.Pt.Y && op1.Prev != op1b && op1.Prev != op2 {
			op1 = op1.Prev
		}
		for op1b.Next.Pt.Y == op1b.Pt.Y && op1b.Next != op1 && op1b.Next != op2 {
			op1b = op1b.Next
		}
		if op1b.Next == op1 || op1b.Next == op2 {
			return false //a flat 'polygon'
		}

		op2b = op2
		for op2.Prev.Pt.Y == op2.Pt.Y && op2.Prev != op2b && op2.Prev != op1b {
			op2 = op2.Prev
		}
		for op2b.Next.Pt.Y == op2b.Pt.Y && op2b.Next != op2 && op2b.Next != op1 {
			op2b = op2b.Next
		}
		if op2b.Next == op2 || op2b.Next == op1 {
			return false //a flat 'polygon'
		}

		var Left, Right CInt
		//Op1 -. Op1b & Op2 -. Op2b are the extremites of the horizontal edges
		if !c.GetOverlap(op1.Pt.X, op1b.Pt.X, op2.Pt.X, op2b.Pt.X, &Left, &Right) {
			return false
		}

		//DiscardLeftSide: when overlapping edges are joined, a spike will created
		//which needs to be cleaned up. However, we don't want Op1 or Op2 caught up
		//on the discard Side as either may still be needed for other joins ...
		var Pt *IntPoint
		var DiscardLeftSide bool
		if op1.Pt.X >= Left && op1.Pt.X <= Right {
			Pt = op1.Pt
			DiscardLeftSide = (op1.Pt.X > op1b.Pt.X)
		} else if op2.Pt.X >= Left && op2.Pt.X <= Right {
			Pt = op2.Pt
			DiscardLeftSide = (op2.Pt.X > op2b.Pt.X)
		} else if op1b.Pt.X >= Left && op1b.Pt.X <= Right {
			Pt = op1b.Pt
			DiscardLeftSide = op1b.Pt.X > op1.Pt.X
		} else {
			Pt = op2b.Pt
			DiscardLeftSide = (op2b.Pt.X > op2.Pt.X)
		}
		j.OutPt1 = op1
		j.OutPt2 = op2
		return c.JoinHorz(op1, op1b, op2, op2b, Pt, DiscardLeftSide)
	} else {
		//nb: For non-horizontal joins ...
		//    1. Jr.OutPt1.Pt.Y == Jr.OutPt2.Pt.Y
		//    2. Jr.OutPt1.Pt > Jr.OffPt.Y

		//make sure the polygons are correctly oriented ...
		op1b = op1.Next
		for (op1b.Pt == op1.Pt) && (op1b != op1) {
			op1b = op1b.Next
		}
		Reverse1 := ((op1b.Pt.Y > op1.Pt.Y) ||
			!c.SlopesEqual3(op1.Pt, op1b.Pt, j.OffPt, c.m_UseFullRange))
		if Reverse1 {
			op1b = op1.Prev
			for (op1b.Pt == op1.Pt) && (op1b != op1) {
				op1b = op1b.Prev
			}
			if (op1b.Pt.Y > op1.Pt.Y) ||
				!c.SlopesEqual3(op1.Pt, op1b.Pt, j.OffPt, c.m_UseFullRange) {
				return false
			}
		}
		op2b = op2.Next
		for (op2b.Pt == op2.Pt) && (op2b != op2) {
			op2b = op2b.Next
		}
		Reverse2 := ((op2b.Pt.Y > op2.Pt.Y) ||
			!c.SlopesEqual3(op2.Pt, op2b.Pt, j.OffPt, c.m_UseFullRange))
		if Reverse2 {
			op2b = op2.Prev
			for (op2b.Pt == op2.Pt) && (op2b != op2) {
				op2b = op2b.Prev
			}
			if (op2b.Pt.Y > op2.Pt.Y) ||
				!c.SlopesEqual3(op2.Pt, op2b.Pt, j.OffPt, c.m_UseFullRange) {
				return false
			}
		}

		if (op1b == op1) || (op2b == op2) || (op1b == op2b) ||
			((outRec1 == outRec2) && (Reverse1 == Reverse2)) {
			return false
		}

		if Reverse1 {
			op1b = c.DupOutPt(op1, false)
			op2b = c.DupOutPt(op2, true)
			op1.Prev = op2
			op2.Next = op1
			op1b.Next = op2b
			op2b.Prev = op1b
			j.OutPt1 = op1
			j.OutPt2 = op1b
			return true
		} else {
			op1b = c.DupOutPt(op1, true)
			op2b = c.DupOutPt(op2, false)
			op1.Next = op2
			op2.Prev = op1
			op1b.Prev = op2b
			op2b.Next = op1b
			j.OutPt1 = op1
			j.OutPt2 = op1b
			return true
		}
	}
}

//----------------------------------------------------------------------

//returns 0 if false, +1 if true, -1 if pt ON polygon boundary
//See "The Point in Polygon Problem for Arbitrary Polygons" by Hormann & Agathos
//http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.88.5498&rep=rep1&type=pdf
func PointInPolygon(pt *IntPoint, path Path) int {
	result := 0
	cnt := len(path)
	if cnt < 3 {
		return 0
	}
	ip := path[0]
	for i := 1; i <= cnt; i++ {
		var ipNext *IntPoint
		if i == cnt {
			ipNext = path[0]
		} else {
			ipNext = path[i]
		}
		if ipNext.Y == pt.Y {
			if (ipNext.X == pt.X) || (ip.Y == pt.Y &&
				((ipNext.X > pt.X) == (ip.X < pt.X))) {
				return -1
			}
		}
		if (ip.Y < pt.Y) != (ipNext.Y < pt.Y) {
			if ip.X >= pt.X {
				if ipNext.X > pt.X {
					result = 1 - result
				} else {
					d := float64(ip.X-pt.X)*float64(ipNext.Y-pt.Y) -
						float64(ipNext.X-pt.X)*float64(ip.Y-pt.Y)
					if d == 0 {
						return -1
					} else if (d > 0) == (ipNext.Y > ip.Y) {
						result = 1 - result
					}
				}
			} else {
				if ipNext.X > pt.X {
					d := float64(ip.X-pt.X)*float64(ipNext.Y-pt.Y) -
						float64(ipNext.X-pt.X)*float64(ip.Y-pt.Y)
					if d == 0 {
						return -1
					} else if (d > 0) == (ipNext.Y > ip.Y) {
						result = 1 - result
					}
				}
			}
		}
		ip = ipNext
	}
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) pointInPolygon(pt *IntPoint, op *OutPt) int {
	//returns 0 if false, +1 if true, -1 if pt ON polygon boundary
	//See "The Point in Polygon Problem for Arbitrary Polygons" by Hormann & Agathos
	//http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.88.5498&rep=rep1&type=pdf
	result := 0
	startOp := op
	ptx := pt.X
	pty := pt.Y
	poly0x := op.Pt.X
	poly0y := op.Pt.Y
	for {
		op = op.Next
		poly1x := op.Pt.X
		poly1y := op.Pt.Y

		if poly1y == pty {
			if (poly1x == ptx) || (poly0y == pty &&
				((poly1x > ptx) == (poly0x < ptx))) {
				return -1
			}
		}
		if (poly0y < pty) != (poly1y < pty) {
			if poly0x >= ptx {
				if poly1x > ptx {
					result = 1 - result
				} else {
					d := float64(poly0x-ptx)*float64(poly1y-pty) -
						float64(poly1x-ptx)*float64(poly0y-pty)
					if d == 0 {
						return -1
					}
					if (d > 0) == (poly1y > poly0y) {
						result = 1 - result
					}
				}
			} else {
				if poly1x > ptx {
					d := float64(poly0x-ptx)*float64(poly1y-pty) -
						float64(poly1x-ptx)*float64(poly0y-pty)
					if d == 0 {
						return -1
					}
					if (d > 0) == (poly1y > poly0y) {
						result = 1 - result
					}
				}
			}
		}
		poly0x = poly1x
		poly0y = poly1y
		if startOp == op {
			break
		}
	}
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) Poly2ContainsPoly1(outPt1, outPt2 *OutPt) bool {
	op := outPt1
	for {
		res := c.pointInPolygon(op.Pt, outPt2)
		if res >= 0 {
			return res != 0
		}
		op = op.Next
		if op == outPt1 {
			break
		}
	}
	return true
}

//----------------------------------------------------------------------

func (c *Clipper) FixupFirstLefts1(OldOutRec, NewOutRec *OutRec) {
	for i := 0; i < len(c.m_PolyOuts); i++ {
		outRec := c.m_PolyOuts[i]
		if outRec.Pts != nil && outRec.FirstLeft == OldOutRec {
			if c.Poly2ContainsPoly1(outRec.Pts, NewOutRec.Pts) {
				outRec.FirstLeft = NewOutRec
			}
		}
	}
}

//----------------------------------------------------------------------

func (c *Clipper) FixupFirstLefts2(OldOutRec, NewOutRec *OutRec) {
	for _, outRec := range c.m_PolyOuts {
		if outRec.FirstLeft == OldOutRec {
			outRec.FirstLeft = NewOutRec
		}
	}
}

//----------------------------------------------------------------------

func (c *Clipper) ParseFirstLeft(FirstLeft *OutRec) *OutRec {
	for FirstLeft != nil && FirstLeft.Pts == nil {
		FirstLeft = FirstLeft.FirstLeft
	}
	return FirstLeft
}

//------------------------------------------------------------------------------

func (c *Clipper) JoinCommonEdges() {
	for i := 0; i < len(c.m_Joins); i++ {
		join := c.m_Joins[i]

		outRec1 := c.GetOutRec(join.OutPt1.Idx)
		outRec2 := c.GetOutRec(join.OutPt2.Idx)

		if outRec1.Pts == nil || outRec2.Pts == nil {
			continue
		}

		//get the polygon fragment with the correct hole state (FirstLeft)
		//before calling JoinPoints() ...
		var holeStateRec *OutRec
		if outRec1 == outRec2 {
			holeStateRec = outRec1
		} else if c.Param1RightOfParam2(outRec1, outRec2) {
			holeStateRec = outRec2
		} else if c.Param1RightOfParam2(outRec2, outRec1) {
			holeStateRec = outRec1
		} else {
			holeStateRec = c.GetLowermostRec(outRec1, outRec2)
		}

		if !c.JoinPoints(join, outRec1, outRec2) {
			continue
		}

		if outRec1 == outRec2 {
			//instead of joining two polygons, we've just created a new one by
			//splitting one polygon into two.
			outRec1.Pts = join.OutPt1
			outRec1.BottomPt = nil
			outRec2 = c.CreateOutRec()
			outRec2.Pts = join.OutPt2

			//update all OutRec2.Pts Idx's ...
			c.UpdateOutPtIdxs(outRec2)

			//We now need to check every OutRec.FirstLeft pointer. If it points
			//to OutRec1 it may need to point to OutRec2 instead ...
			if c.m_UsingPolyTree {
				for j := 0; j < len(c.m_PolyOuts)-1; j++ {
					oRec := c.m_PolyOuts[j]
					if oRec.Pts == nil || c.ParseFirstLeft(oRec.FirstLeft) != outRec1 ||
						oRec.IsHole == outRec1.IsHole {
						continue
					}
					if c.Poly2ContainsPoly1(oRec.Pts, join.OutPt2) {
						oRec.FirstLeft = outRec2
					}
				}
			}

			if c.Poly2ContainsPoly1(outRec2.Pts, outRec1.Pts) {
				//outRec2 is contained by outRec1 ...
				outRec2.IsHole = !outRec1.IsHole
				outRec2.FirstLeft = outRec1

				//fixup FirstLeft pointers that may need reassigning to OutRec1
				if c.m_UsingPolyTree {
					c.FixupFirstLefts2(outRec2, outRec1)
				}

				if (outRec2.IsHole != c.ReverseSolution) == (c.area(outRec2) > 0) {
					c.ReversePolyPtLinks(outRec2.Pts)
				}

			} else if c.Poly2ContainsPoly1(outRec1.Pts, outRec2.Pts) {
				//outRec1 is contained by outRec2 ...
				outRec2.IsHole = outRec1.IsHole
				outRec1.IsHole = !outRec2.IsHole
				outRec2.FirstLeft = outRec1.FirstLeft
				outRec1.FirstLeft = outRec2

				//fixup FirstLeft pointers that may need reassigning to OutRec1
				if c.m_UsingPolyTree {
					c.FixupFirstLefts2(outRec1, outRec2)
				}

				if (outRec1.IsHole != c.ReverseSolution) == (c.area(outRec1) > 0) {
					c.ReversePolyPtLinks(outRec1.Pts)
				}
			} else {
				//the 2 polygons are completely separate ...
				outRec2.IsHole = outRec1.IsHole
				outRec2.FirstLeft = outRec1.FirstLeft

				//fixup FirstLeft pointers that may need reassigning to OutRec2
				if c.m_UsingPolyTree {
					c.FixupFirstLefts1(outRec1, outRec2)
				}
			}

		} else {
			//joined 2 polygons together ...

			outRec2.Pts = nil
			outRec2.BottomPt = nil
			outRec2.Idx = outRec1.Idx

			outRec1.IsHole = holeStateRec.IsHole
			if holeStateRec == outRec2 {
				outRec1.FirstLeft = outRec2.FirstLeft
			}
			outRec2.FirstLeft = outRec1

			//fixup FirstLeft pointers that may need reassigning to OutRec1
			if c.m_UsingPolyTree {
				c.FixupFirstLefts2(outRec2, outRec1)
			}
		}
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) UpdateOutPtIdxs(outrec *OutRec) {
	op := outrec.Pts
	for {
		op.Idx = outrec.Idx
		op = op.Prev
		if op == outrec.Pts {
			break
		}
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) DoSimplePolygons() {
	for i := 0; i < len(c.m_PolyOuts); i++ {
		outrec := c.m_PolyOuts[i]
		op := outrec.Pts
		if op == nil {
			continue
		}
		for { //for each Pt in Polygon until duplicate found do ...
			op2 := op.Next
			for op2 != outrec.Pts {
				if op.Pt.Equals(op2.Pt) && op2.Next != op && op2.Prev != op {
					//split the polygon into two ...
					op3 := op.Prev
					op4 := op2.Prev
					op.Prev = op4
					op4.Next = op
					op2.Prev = op3
					op3.Next = op2

					outrec.Pts = op
					outrec2 := c.CreateOutRec()
					outrec2.Pts = op2
					c.UpdateOutPtIdxs(outrec2)
					if c.Poly2ContainsPoly1(outrec2.Pts, outrec.Pts) {
						//OutRec2 is contained by OutRec1 ...
						outrec2.IsHole = !outrec.IsHole
						outrec2.FirstLeft = outrec
					} else {
						if c.Poly2ContainsPoly1(outrec.Pts, outrec2.Pts) {
							//OutRec1 is contained by OutRec2 ...
							outrec2.IsHole = outrec.IsHole
							outrec.IsHole = !outrec2.IsHole
							outrec2.FirstLeft = outrec.FirstLeft
							outrec.FirstLeft = outrec2
						} else {
							//the 2 polygons are separate ...
							outrec2.IsHole = outrec.IsHole
							outrec2.FirstLeft = outrec.FirstLeft
						}
					}
					op2 = op //ie get ready for the next iteration
				}
				op2 = op2.Next
			}
			op = op.Next
			if op == outrec.Pts {
				break
			}
		}
	}
}

//------------------------------------------------------------------------------

func Area(poly Path) float64 {
	cnt := len(poly)
	if cnt < 3 {
		return 0
	}
	a := float64(0)
	j := cnt - 1
	for i := 0; i < cnt; i++ {
		a += (float64(poly[j].X) + float64(poly[i].X)) *
			(float64(poly[j].Y - poly[i].Y))
		j = i
	}
	return -a * 0.5
}

func AreaCombined(polygons Paths) float64 {
	a := 0.
	for _, polygon := range polygons {
		a += Area(polygon)
	}
	return a
}

//------------------------------------------------------------------------------

func (c *Clipper) area(outRec *OutRec) float64 {
	op := outRec.Pts
	if op == nil {
		return 0
	}
	a := float64(0)
	for {
		a = a + float64(op.Prev.Pt.X+op.Pt.X)*float64(op.Prev.Pt.Y-op.Pt.Y)
		op = op.Next
		if op == outRec.Pts {
			break
		}
	}
	return a * 0.5
}

//------------------------------------------------------------------------------
// SimplifyPolygon functions ...
// Convert self-intersecting polygons into simple polygons
//------------------------------------------------------------------------------
// default fillType = PftEvenOdd
func (c *Clipper) SimplifyPolygon(poly Path,
	fillType PolyFillType) Paths {
	var result Paths
	c2 := NewClipper(IoNone)
	c2.StrictlySimple = true
	c2.AddPath(poly, PtSubject, true)
	result, _ = c2.Execute1(CtUnion, fillType, fillType)
	return result
}

//------------------------------------------------------------------------------

// default fillType = PftEvenOdd
func (c *Clipper) SimplifyPolygons(polys Paths,
	fillType PolyFillType) Paths {
	var result Paths
	c2 := NewClipper(IoNone)
	c2.StrictlySimple = true
	c2.AddPaths(polys, PtSubject, true)
	result, _ = c2.Execute1(CtUnion, fillType, fillType)
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) DistanceSqrd(pt1, pt2 *IntPoint) float64 {
	dx := float64(pt1.X - pt2.X)
	dy := float64(pt1.Y - pt2.Y)
	return (dx*dx + dy*dy)
}

//------------------------------------------------------------------------------

func (c *Clipper) DistanceFromLineSqrd(pt, ln1, ln2 *IntPoint) float64 {
	//The equation of a line in general form (Ax + By + C = 0)
	//given 2 points (x¹,y¹) & (x²,y²) is ...
	//(y¹ - y²)x + (x² - x¹)y + (y² - y¹)x¹ - (x² - x¹)y¹ = 0
	//A = (y¹ - y²); B = (x² - x¹); C = (y² - y¹)x¹ - (x² - x¹)y¹
	//perpendicular distance of point (x³,y³) = (Ax³ + By³ + C)/Sqrt(A² + B²)
	//see http://en.wikipedia.org/wiki/Perpendicular_distance
	A := float64(ln1.Y - ln2.Y)
	B := float64(ln2.X - ln1.X)
	C := A*float64(ln1.X) + B*float64(ln1.Y)
	C = A*float64(pt.X) + B*float64(pt.Y) - C
	return (C * C) / (A*A + B*B)
}

//---------------------------------------------------------------------------

func (c *Clipper) SlopesNearCollinear(pt1,
	pt2, pt3 *IntPoint, distSqrd float64) bool {
	return c.DistanceFromLineSqrd(pt2, pt1, pt3) < distSqrd
}

//------------------------------------------------------------------------------

func (c *Clipper) PointsAreClose(pt1, pt2 *IntPoint, distSqrd float64) bool {
	dx := float64(pt1.X - pt2.X)
	dy := float64(pt1.Y - pt2.Y)
	return ((dx*dx)+(dy*dy) <= distSqrd)
}

//------------------------------------------------------------------------------

func (c *Clipper) ExcludeOp(op *OutPt) *OutPt {
	result := op.Prev
	result.Next = op.Next
	op.Next.Prev = result
	result.Idx = 0
	return result
}

//------------------------------------------------------------------------------

//distance = proximity in units/pixels below which vertices will be stripped.
//Default ~= sqrt(2) so when adjacent vertices or semi-adjacent vertices have
//both x & y coords within 1 unit, then the second vertex will be stripped.
// default distance=1.415
func (c *Clipper) CleanPolygon(path Path, distance float64) Path {

	cnt := len(path)

	if cnt == 0 {
		return NewPath()
	}

	outPts := make([]*OutPt, cnt)
	for i := 0; i < cnt; i++ {
		outPts[i] = new(OutPt)
	}

	for i := 0; i < cnt; i++ {
		outPts[i].Pt = path[i]
		outPts[i].Next = outPts[(i+1)%cnt]
		outPts[i].Next.Prev = outPts[i]
		outPts[i].Idx = 0
	}

	distSqrd := distance * distance
	op := outPts[0]
	for op.Idx == 0 && op.Next != op.Prev {
		if c.PointsAreClose(op.Pt, op.Prev.Pt, distSqrd) {
			op = c.ExcludeOp(op)
			cnt--
		} else if c.PointsAreClose(op.Prev.Pt, op.Next.Pt, distSqrd) {
			c.ExcludeOp(op.Next)
			op = c.ExcludeOp(op)
			cnt -= 2
		} else if c.SlopesNearCollinear(op.Prev.Pt, op.Pt, op.Next.Pt, distSqrd) {
			op = c.ExcludeOp(op)
			cnt--
		} else {
			op.Idx = 1
			op = op.Next
		}
	}

	if cnt < 3 {
		cnt = 0
	}
	result := Path(make([]*IntPoint, cnt))
	for i := 0; i < cnt; i++ {
		result[i] = op.Pt
		op = op.Next
	}
	outPts = nil
	return result
}

//------------------------------------------------------------------------------
// default distance = 1.415
func (c *Clipper) CleanPolygons(polys Paths,
	distance float64) Paths {
	result := Paths(make([]Path, len(polys)))
	for i := 0; i < len(polys); i++ {
		result[i] = c.CleanPolygon(polys[i], distance)
	}
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) Minkowski(pattern, path Path, IsSum, IsClosed bool) Paths {
	var delta int
	if IsClosed {
		delta = 1
	} else {
		delta = 0
	}
	polyCnt := len(pattern)
	pathCnt := len(path)
	result := Paths(make([]Path, pathCnt))
	if IsSum {
		for i := 0; i < pathCnt; i++ {
			result[i] = make([]*IntPoint, polyCnt)
			for j, ip := range pattern {
				result[i][j] = &IntPoint{path[i].X + ip.X, path[i].Y + ip.Y}
			}
		}
	} else {
		for i := 0; i < pathCnt; i++ {
			result[i] = make([]*IntPoint, polyCnt)
			for j, ip := range pattern {
				result[i][j] = &IntPoint{path[i].X - ip.X, path[i].Y - ip.Y}
			}
		}
	}

	quads := Paths(make([]Path, 0, (pathCnt+delta)*(polyCnt+1)))
	for i := 0; i < pathCnt-1+delta; i++ {
		for j := 0; j < polyCnt; j++ {
			quad := Path(make([]*IntPoint, 4))
			quad[0] = result[i%pathCnt][j%polyCnt]
			quad[1] = result[(i+1)%pathCnt][j%polyCnt]
			quad[2] = result[(i+1)%pathCnt][(j+1)%polyCnt]
			quad[3] = result[i%pathCnt][(j+1)%polyCnt]
			if !Orientation(quad) {
				reversePath(quad)
			}
			quads = append(quads, quad)
		}
	}
	return quads
}

//------------------------------------------------------------------------------

func (c *Clipper) MinkowskiSum(pattern, path Path, pathIsClosed bool) Paths {
	paths := c.Minkowski(pattern, path, true, pathIsClosed)
	c2 := NewClipper(IoNone)
	c2.AddPaths(paths, PtSubject, true)
	paths, _ = c2.Execute1(CtUnion, PftNonZero, PftNonZero)
	return paths
}

//------------------------------------------------------------------------------

func (c *Clipper) TranslatePath(path Path, delta *IntPoint) Path {
	outPath := Path(make([]*IntPoint, len(path)))
	for i := 0; i < len(path); i++ {
		outPath[i] = &IntPoint{path[i].X + delta.X, path[i].Y + delta.Y}
	}
	return outPath
}

//------------------------------------------------------------------------------

func (c *Clipper) MinkowskiSumAll(pattern Path, paths Paths, pathIsClosed bool) Paths {
	var solution Paths
	c2 := NewClipper(IoNone)
	for i := 0; i < len(paths); i++ {
		tmp := c.Minkowski(pattern, paths[i], true, pathIsClosed)
		c.AddPaths(tmp, PtSubject, true)
		if pathIsClosed {
			path := c.TranslatePath(paths[i], pattern[0])
			c.AddPath(path, PtClip, true)
		}
	}
	solution, _ = c2.Execute1(CtUnion, PftNonZero, PftNonZero)
	return solution
}

//------------------------------------------------------------------------------

func (c *Clipper) MinkowskiDiff(poly1, poly2 Path) Paths {
	paths := c.Minkowski(poly1, poly2, false, true)
	c2 := NewClipper(IoNone)
	c2.AddPaths(paths, PtSubject, true)
	paths, _ = c2.Execute1(CtUnion, PftNonZero, PftNonZero)
	return paths
}

//------------------------------------------------------------------------------

type NodeType int

const (
	ntAny NodeType = iota
	ntOpen
	ntClosed
)

func (c *Clipper) PolyTreeToPaths(polytree *PolyTree) Paths {

	result := Paths(make([]Path, 0, polytree.Total()))
	c.AddPolyNodeToPaths(polytree.toPolyNode(), ntAny, result)
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) AddPolyNodeToPaths(polynode *PolyNode, nt NodeType, paths Paths) {
	match := true
	switch nt {
	case ntOpen:
		return
	case ntClosed:
		match = !polynode.IsOpen
		break
	default:
		break
	}

	if len(polynode.m_polygon) > 0 && match {
		paths = append(paths, polynode.m_polygon)
	}
	for _, pn := range polynode.m_Childs {
		c.AddPolyNodeToPaths(pn, nt, paths)
	}
}

//------------------------------------------------------------------------------

func (c *Clipper) OpenPathsFromPolyTree(polytree *PolyTree) Paths {
	result := Paths(make([]Path, 0, polytree.ChildCount()))
	for i := 0; i < polytree.ChildCount(); i++ {
		if polytree.m_Childs[i].IsOpen {
			result = append(result, polytree.m_Childs[i].m_polygon)
		}
	}
	return result
}

//------------------------------------------------------------------------------

func (c *Clipper) ClosedPathsFromPolyTree(polytree *PolyTree) Paths {
	result := Paths(make([]Path, 0, polytree.Total()))
	c.AddPolyNodeToPaths(polytree.toPolyNode(), ntClosed, result)
	return result
}

//------------------------------------------------------------------------------

const two_pi = math.Pi * 2
const def_arc_tolerance = 0.25

type ClipperOffset struct {
	m_destPolys                   Paths
	m_srcPoly                     Path
	m_destPoly                    Path
	m_normals                     []*DoublePoint
	m_delta, m_sinA, m_sin, m_cos float64
	m_miterLim, m_StepsPerRad     float64

	m_lowest    *IntPoint
	m_polyNodes *PolyNode

	ArcTolerance float64
	MiterLimit   float64
}

func NewClipperOffset() *ClipperOffset {
	co := new(ClipperOffset)
	co.m_normals = make([]*DoublePoint, 0)
	co.m_polyNodes = NewPolyNode()
	co.MiterLimit = 2.0
	co.ArcTolerance = def_arc_tolerance
	co.m_lowest = NewIntPoint(-1, 0)
	return co
}

func (co *ClipperOffset) Clear() {
	co.m_polyNodes.m_Childs = make([]*PolyNode, 0)
	co.m_lowest.X = -1
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) Round(value float64) CInt {
	if value < 0 {
		return CInt(value - 0.5)
	} else {
		return CInt(value + 0.5)
	}
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) AddPath(path Path, joinType JoinType, endType EndType) {
	highI := len(path) - 1
	if highI < 0 {
		return
	}
	newNode := NewPolyNode()
	newNode.m_jointype = joinType
	newNode.m_endtype = endType

	//strip duplicate points from path and also get index to the lowest point ...
	if endType == EtClosedLine || endType == EtClosedPolygon {
		for highI > 0 && path[0] == path[highI] {
			highI--
		}
	}
	newNode.m_polygon = Path(make([]*IntPoint, 0, highI+1))
	newNode.m_polygon = append(newNode.m_polygon, path[0])
	j := 0
	k := 0
	for i := 1; i <= highI; i++ {
		if newNode.m_polygon[j] != path[i] {
			j++
			newNode.m_polygon = append(newNode.m_polygon, path[i])
			if path[i].Y > newNode.m_polygon[k].Y ||
				(path[i].Y == newNode.m_polygon[k].Y &&
					path[i].X < newNode.m_polygon[k].X) {
				k = j
			}
		}
	}
	if endType == EtClosedPolygon && j < 2 {
		return
	}

	co.m_polyNodes.AddChild(newNode)

	//if this path's lowest pt is lower than all the others then update m_lowest
	if endType != EtClosedPolygon {
		return
	}
	if co.m_lowest.X < 0 {
		co.m_lowest = &IntPoint{CInt(co.m_polyNodes.ChildCount() - 1), CInt(k)}
	} else {
		ip := co.m_polyNodes.m_Childs[int(co.m_lowest.X)].m_polygon[int(co.m_lowest.Y)]
		if newNode.m_polygon[k].Y > ip.Y ||
			(newNode.m_polygon[k].Y == ip.Y &&
				newNode.m_polygon[k].X < ip.X) {
			co.m_lowest = &IntPoint{CInt(co.m_polyNodes.ChildCount() - 1), CInt(k)}
		}
	}
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) AddPaths(paths Paths, joinType JoinType, endType EndType) {
	for _, p := range paths {
		co.AddPath(p, joinType, endType)
	}
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) FixOrientations() {
	//fixup orientations of all closed paths if the orientation of the
	//closed path with the lowermost vertex is wrong ...
	if co.m_lowest.X >= 0 &&
		!Orientation(co.m_polyNodes.m_Childs[int(co.m_lowest.X)].m_polygon) {
		for i := 0; i < co.m_polyNodes.ChildCount(); i++ {
			node := co.m_polyNodes.m_Childs[i]
			if node.m_endtype == EtClosedPolygon ||
				(node.m_endtype == EtClosedLine &&
					Orientation(node.m_polygon)) {
				reversePath(node.m_polygon)
			}
		}
	} else {
		for i := 0; i < co.m_polyNodes.ChildCount(); i++ {
			node := co.m_polyNodes.m_Childs[i]
			if node.m_endtype == EtClosedLine &&
				!Orientation(node.m_polygon) {
				reversePath(node.m_polygon)
			}
		}
	}
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) GetUnitNormal(pt1, pt2 *IntPoint) *DoublePoint {
	dx := float64(pt2.X - pt1.X)
	dy := float64(pt2.Y - pt1.Y)
	if (dx == 0) && (dy == 0) {
		return new(DoublePoint)
	}

	f := 1 * 1.0 / math.Sqrt(dx*dx+dy*dy)
	dx *= f
	dy *= f
	return &DoublePoint{dy, -dx}
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) DoOffset(delta float64) {
	co.m_destPolys = Paths(make([]Path, 0, co.m_polyNodes.ChildCount()*2))
	co.m_delta = delta

	//if Zero offset, just copy any CLOSED polygons to m_p and return ...
	if near_zero(delta) {
		for i := 0; i < co.m_polyNodes.ChildCount(); i++ {
			node := co.m_polyNodes.m_Childs[i]
			if node.m_endtype == EtClosedPolygon {
				co.m_destPolys = append(co.m_destPolys, node.m_polygon)
			}
		}
		return
	}

	//see offset_triginometry3.svg in the documentation folder ...
	if co.MiterLimit > 2 {
		co.m_miterLim = 2 / (co.MiterLimit * co.MiterLimit)
	} else {
		co.m_miterLim = 0.5
	}

	var y float64
	if co.ArcTolerance <= 0.0 {
		y = def_arc_tolerance
	} else if co.ArcTolerance > math.Abs(delta)*def_arc_tolerance {
		y = math.Abs(delta) * def_arc_tolerance
	} else {
		y = co.ArcTolerance
	}
	//see offset_triginometry2.svg in the documentation folder ...
	steps := math.Pi / math.Acos(1-y/math.Abs(delta))
	co.m_sin = math.Sin(two_pi / steps)
	co.m_cos = math.Cos(two_pi / steps)
	co.m_StepsPerRad = steps / two_pi
	if delta < 0.0 {
		co.m_sin = -co.m_sin
	}

	for i := 0; i < co.m_polyNodes.ChildCount(); i++ {
		node := co.m_polyNodes.m_Childs[i]
		co.m_srcPoly = node.m_polygon

		length := len(co.m_srcPoly)

		if length == 0 || (delta <= 0 && (length < 3 ||
			node.m_endtype != EtClosedPolygon)) {
			continue
		}

		co.m_destPoly = Path(make([]*IntPoint, 0, int(steps)+4))

		if length == 1 {
			if node.m_jointype == JtRound {
				X := 1.0
				Y := 0.0
				for j := 1; j <= int(steps); j++ {
					co.m_destPoly = append(co.m_destPoly, &IntPoint{
						co.Round(float64(co.m_srcPoly[0].X) + X*delta),
						co.Round(float64(co.m_srcPoly[0].Y) + Y*delta)})
					X2 := X
					X = X*co.m_cos - co.m_sin*Y
					Y = X2*co.m_sin + Y*co.m_cos
				}
			} else {
				X := -1.0
				Y := -1.0
				for j := 0; j < 4; j++ {
					co.m_destPoly = append(co.m_destPoly, &IntPoint{
						co.Round(float64(co.m_srcPoly[0].X) + X*delta),
						co.Round(float64(co.m_srcPoly[0].Y) + Y*delta)})
					if X < 0 {
						X = 1
					} else if Y < 0 {
						Y = 1
					} else {
						X = -1
					}
				}
			}
			co.m_destPolys = append(co.m_destPolys, co.m_destPoly)
			continue
		}

		//build m_normals ...
		co.m_normals = make([]*DoublePoint, 0, length)
		for j := 0; j < length-1; j++ {
			co.m_normals = append(co.m_normals,
				co.GetUnitNormal(co.m_srcPoly[j], co.m_srcPoly[j+1]))
		}
		if node.m_endtype == EtClosedLine ||
			node.m_endtype == EtClosedPolygon {
			co.m_normals = append(co.m_normals,
				co.GetUnitNormal(co.m_srcPoly[length-1], co.m_srcPoly[0]))
		} else {
			co.m_normals = append(co.m_normals,
				CopyDoublePoint(co.m_normals[length-2]))
		}

		if node.m_endtype == EtClosedPolygon {
			k := length - 1
			for j := 0; j < length; j++ {
				co.OffsetPoint(j, &k, node.m_jointype)
			}
			co.m_destPolys = append(co.m_destPolys, co.m_destPoly)
		} else if node.m_endtype == EtClosedLine {
			k := length - 1
			for j := 0; j < length; j++ {
				co.OffsetPoint(j, &k, node.m_jointype)
			}
			co.m_destPolys = append(co.m_destPolys, co.m_destPoly)
			co.m_destPoly = Path(make([]*IntPoint, 0))
			//re-build m_normals ...
			n := co.m_normals[length-1]
			for j := length - 1; j > 0; j-- {
				co.m_normals[j] = &DoublePoint{-co.m_normals[j-1].X, -co.m_normals[j-1].Y}
			}
			co.m_normals[0] = &DoublePoint{-n.X, -n.Y}
			k = 0
			for j := length - 1; j >= 0; j-- {
				co.OffsetPoint(j, &k, node.m_jointype)
			}
			co.m_destPolys = append(co.m_destPolys, co.m_destPoly)
		} else {
			k := 0
			for j := 1; j < length-1; j++ {
				co.OffsetPoint(j, &k, node.m_jointype)
			}

			var pt1 *IntPoint
			if node.m_endtype == EtOpenButt {
				j := length - 1
				pt1 = &IntPoint{CInt(co.Round(float64(co.m_srcPoly[j].X) +
					co.m_normals[j].X*delta)),
					CInt(co.Round(float64(co.m_srcPoly[j].Y) +
						co.m_normals[j].Y*delta))}
				co.m_destPoly = append(co.m_destPoly, pt1)
				pt1 = &IntPoint{CInt(co.Round(float64(co.m_srcPoly[j].X) -
					co.m_normals[j].X*delta)),
					CInt(co.Round(float64(co.m_srcPoly[j].Y) -
						co.m_normals[j].Y*delta))}
				co.m_destPoly = append(co.m_destPoly, pt1)
			} else {
				j := length - 1
				k := length - 2
				co.m_sinA = 0
				co.m_normals[j] = &DoublePoint{-co.m_normals[j].X,
					-co.m_normals[j].Y}
				if node.m_endtype == EtOpenSquare {
					co.DoSquare(j, k)
				} else {
					co.DoRound(j, k)
				}
			}

			//re-build m_normals ...
			for j := length - 1; j > 0; j-- {
				co.m_normals[j] = &DoublePoint{-co.m_normals[j-1].X, -co.m_normals[j-1].Y}
			}

			co.m_normals[0] = &DoublePoint{-co.m_normals[1].X, -co.m_normals[1].Y}

			k = length - 1
			for j := k - 1; j > 0; j-- {
				co.OffsetPoint(j, &k, node.m_jointype)
			}

			if node.m_endtype == EtOpenButt {
				pt1 = &IntPoint{CInt(co.Round(float64(co.m_srcPoly[0].X) -
					co.m_normals[0].X*delta)),
					CInt(co.Round(float64(co.m_srcPoly[0].Y) -
						co.m_normals[0].Y*delta))}
				co.m_destPoly = append(co.m_destPoly, pt1)
				pt1 = &IntPoint{CInt(co.Round(float64(co.m_srcPoly[0].X) +
					co.m_normals[0].X*delta)),
					CInt(co.Round(float64(co.m_srcPoly[0].Y) +
						co.m_normals[0].Y*delta))}
				co.m_destPoly = append(co.m_destPoly, pt1)
			} else {
				k = 1
				co.m_sinA = 0
				if node.m_endtype == EtOpenSquare {
					co.DoSquare(0, 1)
				} else {
					co.DoRound(0, 1)
				}
			}
			co.m_destPolys = append(co.m_destPolys, co.m_destPoly)
		}
	}
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) Execute(delta float64) (solution Paths) {
	solution = make([]Path, 0)
	co.FixOrientations()
	co.DoOffset(delta)
	//now clean up 'corners' ...
	clpr := NewClipper(IoNone)
	clpr.AddPaths(co.m_destPolys, PtSubject, true)
	if delta > 0 {
		solution, _ = clpr.Execute1(CtUnion,
			PftPositive, PftPositive)
	} else {
		r := GetBounds(co.m_destPolys)
		outer := Path(make([]*IntPoint, 4))

		outer[0] = &IntPoint{r.left - 10, r.bottom + 10}
		outer[1] = &IntPoint{r.right + 10, r.bottom + 10}
		outer[2] = &IntPoint{r.right + 10, r.top - 10}
		outer[3] = &IntPoint{r.left - 10, r.top - 10}

		clpr.AddPath(outer, PtSubject, true)
		clpr.ReverseSolution = true
		solution, _ = clpr.Execute1(CtUnion, PftNegative, PftNegative)
		if len(solution) > 0 {
			solution = solution[1:]
		}
	}
	return
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) Execute2(delta float64) (solution *PolyTree) {
	solution = NewPolyTree()
	co.FixOrientations()
	co.DoOffset(delta)

	//now clean up 'corners' ...
	clpr := NewClipper(IoNone)
	clpr.AddPaths(co.m_destPolys, PtSubject, true)
	if delta > 0 {
		solution, _ = clpr.Execute2(CtUnion,
			PftPositive, PftPositive)
	} else {
		r := GetBounds(co.m_destPolys)
		outer := Path(make([]*IntPoint, 4))

		outer[0] = &IntPoint{r.left - 10, r.bottom + 10}
		outer[1] = &IntPoint{r.right + 10, r.bottom + 10}
		outer[2] = &IntPoint{r.right + 10, r.top - 10}
		outer[3] = &IntPoint{r.left - 10, r.top - 10}

		clpr.AddPath(outer, PtSubject, true)
		clpr.ReverseSolution = true
		solution, _ = clpr.Execute2(CtUnion, PftNegative, PftNegative)
		//remove the outer PolyNode rectangle ...
		if solution.ChildCount() == 1 && solution.m_Childs[0].ChildCount() > 0 {
			outerNode := solution.m_Childs[0]
			solution.m_Childs[0] = outerNode.m_Childs[0]
			for i := 1; i < outerNode.ChildCount(); i++ {
				solution.AddChild(outerNode.m_Childs[i])
			}
		} else {
			solution.Clear()
		}
	}
	return
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) OffsetPoint(j int, k *int, jointype JoinType) {
	//cross product ...
	co.m_sinA = (co.m_normals[*k].X*co.m_normals[j].Y - co.m_normals[j].X*co.m_normals[*k].Y)

	if math.Abs(co.m_sinA*co.m_delta) < 1.0 {
		//dot product ...
		cosA := (co.m_normals[*k].X*co.m_normals[j].X +
			co.m_normals[j].Y*co.m_normals[*k].Y)
		if cosA > 0 { // angle ==> 0 degrees {
			co.m_destPoly = append(co.m_destPoly, &IntPoint{co.Round(
				float64(co.m_srcPoly[j].X) + co.m_normals[*k].X*co.m_delta),
				co.Round(float64(co.m_srcPoly[j].Y) +
					co.m_normals[*k].Y*co.m_delta)})
			return
		}
		//else angle ==> 180 degrees
	} else if co.m_sinA > 1.0 {
		co.m_sinA = 1.0
	} else if co.m_sinA < -1.0 {
		co.m_sinA = -1.0
	}

	if co.m_sinA*co.m_delta < 0 {
		co.m_destPoly = append(co.m_destPoly, &IntPoint{co.Round(
			float64(co.m_srcPoly[j].X) + co.m_normals[*k].X*co.m_delta),
			co.Round(float64(co.m_srcPoly[j].Y) + co.m_normals[*k].Y*co.m_delta)})
		co.m_destPoly = append(co.m_destPoly, co.m_srcPoly[j])
		co.m_destPoly = append(co.m_destPoly, &IntPoint{co.Round(
			float64(co.m_srcPoly[j].X) + co.m_normals[j].X*co.m_delta),
			co.Round(float64(co.m_srcPoly[j].Y) + co.m_normals[j].Y*co.m_delta)})
	} else {
		switch jointype {
		case JtMiter:
			{
				r := 1 + (co.m_normals[j].X*co.m_normals[*k].X +
					co.m_normals[j].Y*co.m_normals[*k].Y)
				if r >= co.m_miterLim {
					co.DoMiter(j, *k, r)
				} else {
					co.DoSquare(j, *k)
				}
				break
			}
		case JtSquare:
			co.DoSquare(j, *k)
			break
		case JtRound:
			co.DoRound(j, *k)
			break
		}
	}
	*k = j
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) DoSquare(j, k int) {
	dx := math.Tan(math.Atan2(co.m_sinA,
		co.m_normals[k].X*co.m_normals[j].X+
			co.m_normals[k].Y*co.m_normals[j].Y) / 4)
	co.m_destPoly = append(co.m_destPoly, &IntPoint{
		co.Round(float64(co.m_srcPoly[j].X) +
			co.m_delta*(co.m_normals[k].X-co.m_normals[k].Y*dx)),
		co.Round(float64(co.m_srcPoly[j].Y) +
			co.m_delta*(co.m_normals[k].Y+co.m_normals[k].X*dx))})
	co.m_destPoly = append(co.m_destPoly, &IntPoint{
		co.Round(float64(co.m_srcPoly[j].X) +
			co.m_delta*(co.m_normals[j].X+co.m_normals[j].Y*dx)),
		co.Round(float64(co.m_srcPoly[j].Y) +
			co.m_delta*(co.m_normals[j].Y-co.m_normals[j].X*dx))})
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) DoMiter(j, k int, r float64) {
	q := co.m_delta / r
	co.m_destPoly = append(co.m_destPoly, &IntPoint{co.Round(
		float64(co.m_srcPoly[j].X) + (co.m_normals[k].X+co.m_normals[j].X)*q),
		co.Round(float64(co.m_srcPoly[j].Y) +
			(co.m_normals[k].Y+co.m_normals[j].Y)*q)})
}

//------------------------------------------------------------------------------

func (co *ClipperOffset) DoRound(j, k int) {
	a := math.Atan2(co.m_sinA,
		co.m_normals[k].X*co.m_normals[j].X+co.m_normals[k].Y*co.m_normals[j].Y)
	steps := int(co.Round(co.m_StepsPerRad * math.Abs(a)))

	X := co.m_normals[k].X
	Y := co.m_normals[k].Y
	var X2 float64
	for i := 0; i < steps; i++ {
		co.m_destPoly = append(co.m_destPoly, &IntPoint{
			co.Round(float64(co.m_srcPoly[j].X) + X*co.m_delta),
			co.Round(float64(co.m_srcPoly[j].Y) + Y*co.m_delta)})
		X2 = X
		X = X*co.m_cos - co.m_sin*Y
		Y = X2*co.m_sin + Y*co.m_cos
	}
	co.m_destPoly = append(co.m_destPoly, &IntPoint{
		co.Round(float64(co.m_srcPoly[j].X) + co.m_normals[j].X*co.m_delta),
		co.Round(float64(co.m_srcPoly[j].Y) + co.m_normals[j].Y*co.m_delta)})
}

//------------------------------------------------------------------------------

type ClipperException struct {
	s string
}

func NewClipperException(s string) *ClipperException {
	e := new(ClipperException)
	e.s = s
	return e
}

func (e ClipperException) Error() string {
	return e.s
}

func intAbs(i int) int {
	if i > 0 {
		return i
	} else {
		return i * -1
	}
}
func min(a, b CInt) CInt {
	if a < b {
		return a
	} else {
		return b
	}
}
func max(a, b CInt) CInt {
	if a > b {
		return a
	} else {
		return b
	}
}
