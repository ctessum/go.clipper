// +build !use_xyz

package clipper

import "fmt"

type IntPoint struct {
	X CInt
	Y CInt
}

func (p IntPoint) String() string {
	return fmt.Sprintf("{%v, %v}", p.X, p.Y)
}

func NewIntPoint(X, Y CInt) *IntPoint {
	ip := new(IntPoint)
	ip.X, ip.Y = X, Y
	return ip
}
func NewIntPointFromFloat(x, y float64) *IntPoint {
	ip := new(IntPoint)
	ip.X, ip.Y = CInt(x), CInt(y)
	return ip
}

func (ip *IntPoint) Copy() IntPoint {
	return IntPoint{
		X: ip.X,
		Y: ip.Y,
	}
}

func (c *ClipperBase) ReverseHorizontal(e *TEdge) {
	//swap horizontal edges' top and bottom x's so they follow the natural
	//progression of the bounds - ie so their xbots will align with the
	//adjoining lower edge. [Helpful in the ProcessHorizontal() method.]
	c.Swap(&e.Top.X, &e.Bot.X)
}

func (c *Clipper) SetZ(pt *IntPoint, e1, e2 *TEdge) {
	return
}
