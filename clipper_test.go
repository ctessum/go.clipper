package clipper


import (
	"math"
	"fmt"
	"os"
)
 
//===============================================================================
//===============================================================================

func LoadFile1(lines []string) [][]*Point {
    // File type 1: first line is total polygons count and subsequent lines 
    // contain the polygon vertex count followed by its coords 
        polygons = make([][]*Point,0)
        poly = make([]*Point,0)
		for _,l :=range lines {
			vals := re.split(' |, |,', l.strip())
            if len(vals) < 2 {  
                if (len(poly)  > 2) {
                    polygons = append(polygons, poly)
				}
                poly = make([]*Point,0)
			} else { 
                poly = append(poly &Point{int(vals[0]), int(vals[1])})
			}
		}
        if (len(poly)  > 2) {
            polygons = append(polygons, poly)
		}
        return polygons
}
//===============================================================================

func LoadFile2(lines []string) [][]*Point {
    // File type 2: vertex coords on consecutive lines for each polygon 
    // where each polygon is separated by an empty line 
        polygons = make([][]*Point,0)
		for l := range lines {
            l = l.strip()
            if (l == "") { 
                if (len(poly)  > 2){
                    polygons = append(polygons,poly)
				}
        poly = make([]*Point,0)
			} else {
				vals := re.split(' |, |,', l)
                poly = append(poly,&Point{int(vals[0]), int(vals[1])})
			}
		}
        if (len(poly)  > 2) {
            polygons = append(polygons, poly)
		}
        return polygons
}
//===============================================================================

func LoadFile(filename string) [][]*Point {
	f, err := os.Open(filename)
	if err != nil {
		return nil 
	}
	defer f.Close()
	lines := f.readlines()
        // pick file type from format of first line ...
        if len(lines) == 0 { return nil
	}else if not ',' in lines[0] { return LoadFile1(lines)
}else { return LoadFile2(lines)}
}
//===============================================================================
    
func SaveToFile(filename string, polys [][]*Point, scale float64) {
	invScale := 1.0 / scale
	f,err := os.Create(filename)
	if err!=nil{panic(err)}
	defer f.Close()
            if invScale == 1 {
				for _,poly :=range polys {
					for _,pt :=range poly{
                        f.write("{0}, {1}\n".format(pt.x, pt.y))
					}
                    f.write("\n")
				}
				} else {
					for _,poly :=range polys {
						for _,pt :=range poly{
                        f.write("{0:.4f}, {1:.4f}\n".format(pt.x * invScale, pt.y * invScale))
					}
                    f.write("\n")
				}
			}
		}
//===============================================================================
    
func RandomPoly(maxWidth, maxHeight, vertCnt int) []*Point {
    result = make([]*Point,vertCnt)
	for i:=0; i<vertCnt; i++ {
        result[i] = &Point{rand.Intn(maxWidth), rand.Intn(maxHeight)}
	}
    return result
}

//===============================================================================
// SVGBuilder
//===============================================================================
class SVGBuilder(object):
    
    func HtmlColor(self, val):
        return "//{0:06x}".format(val & 0xFFFFFF)
    
    func AlphaClr(self, val):
        return "{0:.2f}".format(float(val >> 24)/255)

    class StyleInfo(object):
        func __init__(self):      
            self.fillType = PolyFillType.EvenOdd
            self.brushClr = 0
            self.penClr = 0
            self.penWidth = 0.8
            self.showCoords = False
    
    class StyleInfoPlus(StyleInfo):
        
        func __init__(self):   
            SVGBuilder.StyleInfo.__init__(self)   
            self.polygons = []
            self.textlines = []
        
    func __init__(self):        
        self.GlobalStyle = SVGBuilder.StyleInfo()
        self.PolyInfoList = []
        self.PathHeader = " <path d=\""
        self.PathFooter = "\"\n style=\"fill:{0}; fill-opacity:{1}; fill-rule:{2}; stroke:{3}; stroke-opacity:{4}; stroke-width:{5:.2f};\" filter=\"url(//Gamma)\"/>\n\n"
        self.Header = """<?xml version=\"1.0\" standalone=\"no\"?> 
<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" 
\"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\"> 
\n<svg width=\"{0}px\" height=\"{1}px\" viewBox=\"0 0 {0} {1}\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">
  <funcs>
    <filter id="Gamma">
      <feComponentTransfer>
        <feFuncR type="gamma" amplitude="1" exponent="0.3" offset="0" />
        <feFuncG type="gamma" amplitude="1" exponent="0.3" offset="0" />
        <feFuncB type="gamma" amplitude="1" exponent="0.3" offset="0" />
      </feComponentTransfer>
    </filter>
  </funcs>\n\n"""
    
    func AddPolygon(self, poly, brushColor, penColor):
        if poly is None or len(poly) == 0: return
        pi = self.StyleInfoPlus()
        pi.penWidth = self.GlobalStyle.penWidth
        pi.fillType = self.GlobalStyle.fillType
        pi.showCoords = self.GlobalStyle.showCoords
        pi.brushClr = brushColor
        pi.penClr = penColor        
        pi.polygons.append(poly)
        self.PolyInfoList.append(pi)
    
    func AddPolygons(self, polys, brushColor, penColor):
        if polys is None or len(polys) == 0: return
        pi = self.StyleInfoPlus()
        pi.penWidth = self.GlobalStyle.penWidth
        pi.fillType = self.GlobalStyle.fillType
        pi.showCoords = self.GlobalStyle.showCoords
        pi.brushClr = brushColor
        pi.penClr = penColor        
        pi.polygons = polys
        self.PolyInfoList.append(pi)
    
    func SaveToFile(self, filename, invScale = 1.0, margin = 10):
        if len(self.PolyInfoList) == 0: return False
        if invScale == 0: invScale = 1.0
        if margin < 0: margin = 0
        pi = self.PolyInfoList[0]
        // get bounding rect ...
        left = right = pi.polygons[0][0].x
        top = bottom = pi.polygons[0][0].y
        for pi in self.PolyInfoList:
            for p in pi.polygons:
                for ip in p:
                    if ip.x < left: left = ip.x
                    if ip.x > right: right = ip.x
                    if ip.y < top: top = ip.y
                    if ip.y > bottom: bottom = ip.y
        left *= invScale
        top *= invScale
        right *= invScale
        bottom *= invScale
        offsetX = -left + margin      
        offsetY = -top + margin      
                    
        f = open(filename, 'w')
        m2 = margin * 2
        f.write(self.Header.format(right - left + m2, bottom - top + m2))
        for pi in self.PolyInfoList:
            f.write(self.PathHeader)
            for p in pi.polygons:
                cnt = len(p)
                if cnt < 3: continue
                f.write(" M {0:.2f} {1:.2f}".format(p[0].x * invScale + offsetX, p[0].y * invScale + offsetY))
                for i in range(1,cnt):
                    f.write(" L {0:.2f} {1:.2f}".format(p[i].x * invScale + offsetX, p[i].y * invScale + offsetY))
                f.write(" z")
            fillRule = "evenodd"
            if pi.fillType != PolyFillType.EvenOdd: fillRule = "nonzero"
            f.write(self.PathFooter.format(self.HtmlColor(pi.brushClr), 
                self.AlphaClr(pi.brushClr), fillRule, 
                self.HtmlColor(pi.penClr), self.AlphaClr(pi.penClr),  pi.penWidth))
            
            if (pi.showCoords):
                f.write("<g font-family=\"Verdana\" font-size=\"11\" fill=\"black\">\n\n")
                for p in pi.polygons:
                    cnt = len(p)
                    if cnt < 3: continue
                    for pt in p:
                        x = pt.x * invScale + offsetX
                        y = pt.y * invScale + offsetY
                        f.write("<text x=\"{0}\" y=\"{1}\">{2},{3}</text>\n".format(x, y, pt.x, pt.y))
                    f.write("\n")
                f.write("</g>\n")
    
        f.write("</svg>\n")
        f.close()
        return True

//===============================================================================
// Main entry ...
//===============================================================================
         
scaleExp := 0
scale := math.Pow(10, scaleExp)
invScale := 1.0 / scale

subj, clip := make([][]*Point,0), ([][]*Point,0) 
//load saved subject and clip polygons ...    
//subj = LoadFile('./subj.txt')
//clip = LoadFile('./clip.txt')

// Generate random subject and clip polygons ...
subj = append(subj,RandomPoly(640 * scale, 480 * scale, 100))
clip = append(clip,RandomPoly(640 * scale, 480 * scale, 100))
//SaveToFile('./subj2.txt', subj, scale)
//SaveToFile('./clip2.txt', clip, scale)

// Load the polygons into Clipper and execute the boolean clip op ...
c = NewClipper()
solution = make([][]*Point)
pft := EvenOdd

c.AddPolygons(subj, Subject)
c.AddPolygons(clip, Clip)
result := c.Execute(Intersection, solution, pft, pft)
if !result {
	panic("Problem!")
}

SaveToFile('./solution2.txt', solution, scale)

// Create an SVG file to display what's happened ...
svgBuilder := newSVGBuilder() 
//svgBuilder.GlobalStyle.showCoords = True
svgBuilder.GlobalStyle.fillType = pft
svgBuilder.AddPolygons(subj, 0x402020FF, 0x802020FF)
//svgBuilder.GlobalStyle.showCoords = False
svgBuilder.AddPolygons(clip, 0x40FFFF20, 0x80FF2020)
svgBuilder.GlobalStyle.penWidth = 0.6
svgBuilder.AddPolygons(solution, 0x60138013, 0xFF003300)

holes := make([][]*Point)
for _,poly := range solution {
    if Area(poly) < 0{ holes = append(holes, poly)}
}
svgBuilder.AddPolygons(holes, 0x0, 0xFFFF0000)
svgBuilder.SaveToFile('./test.svg', invScale, 100)

if result { os.startfile('test.svg')} // call(('open', 'test.svg')) // print("finished") // 
else { fmt.Println("failed")}
