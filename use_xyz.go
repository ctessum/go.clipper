// +build use_xyz

package clipper

import "fmt"

type IntPoint struct {
	X CInt
	Y CInt
	Z CInt
}

func (p *IntPoint) String() string {
	return fmt.Sprintf("{%v, %v}", p.X, p.Y)
}

func NewIntPoint(X, Y, Z CInt) *IntPoint {
	ip := new(IntPoint)
	ip.X, ip.Y, ip.Z = X, Y, Z
	return ip
}
func NewIntPointFromFloat(x, y, z float64) *IntPoint {
	ip := new(IntPoint)
	ip.X, ip.Y = CInt(x), CInt(y), CInt(z)
	return ip
}

func (ip *IntPoint) Copy() IntPoint {
	return IntPoint{
		X: ip.X,
		Y: ip.Y,
		Z: ip.Z,
	}
}

func (c *ClipperBase) ReverseHorizontal(e *TEdge) {
	//swap horizontal edges' top and bottom x's so they follow the natural
	//progression of the bounds - ie so their xbots will align with the
	//adjoining lower edge. [Helpful in the ProcessHorizontal() method.]
	Swap(&e.Top.X, &e.Bot.X)
	Swap(&e.Top.Z, &e.Bot.Z)
}

func (c *Clipper) SetZ(pt *IntPoint, e1, e2 *TEdge) {
	if pt.Z != 0 || ZFillFunction == nil {
		return
	} else if pt == e1.Bot {
		pt.Z = e1.Bot.Z
	} else if pt == e1.Top {
		pt.Z = e1.Top.Z
	} else if pt == e2.Bot {
		pt.Z = e2.Bot.Z
	} else if pt == e2.Top {
		pt.Z = e2.Top.Z
	} else {
		c.ZFillFunction(e1.Bot, e1.Top, e2.Bot, e2.Top, &pt)
	}
}

//------------------------------------------------------------------------------
