package main

import (
	"fmt"
	"github.com/jbeda/geom"
	"io"
	"math"
	"os"
	"strings"
)

const (
	DEFAULT_STYLE        = "stroke-width: 0.002; stroke-linecap: round; fill: none"
	DOUBLE_STROKE_OFFSET = 0.002
	CUT_STYLE            = "stroke: black"
	MARK1_STYLE          = "stroke: green"
	MARK2_STYLE          = "stroke: red"
)

////////////////////////////////////////////////////////////////////////////
// SVG serialization helper
type SVG struct {
	Writer io.Writer
}

func NewSVG(w io.Writer) *SVG {
	return &SVG{w}
}

func (svg *SVG) printf(format string, a ...interface{}) (n int, errno error) {
	return fmt.Fprintf(svg.Writer, format, a...)
}

// BUGBUG: not quoting aware
func extraparams(s []string) string {
	ep := ""
	for i := 0; i < len(s); i++ {
		if strings.Index(s[i], "=") > 0 {
			ep += (s[i]) + " "
		} else if len(s[i]) > 0 {
			ep += fmt.Sprintf("style='%s' ", s[i])
		}
	}
	return ep
}

func onezero(b bool) string {
	if b {
		return "1"
	}
	return "0"
}

func (svg *SVG) Start(viewBox geom.Rect, s ...string) {
	svg.printf(`<?xml version="1.0"?>
<svg version="1.1"
     viewBox="%f %f %f %f"
     xmlns="http://www.w3.org/2000/svg" %s>
`, viewBox.Min.X, viewBox.Min.Y, viewBox.Width(), viewBox.Height(), extraparams(s))
}

func (svg *SVG) End() {
	svg.printf("</svg>")
}

func (svg *SVG) Line(p1 geom.Coord, p2 geom.Coord, s ...string) {
	svg.printf("<line x1='%f' y1='%f' x2='%f' y2='%f' %s/>\n", p1.X, p1.Y, p2.X, p2.Y, extraparams(s))
}

func (svg *SVG) Circle(c geom.Coord, r float64, s ...string) {
	svg.printf("<circle cx='%f' cy='%f' r='%f' %s/>\n", c.X, c.Y, r, extraparams(s))
}

func (svg *SVG) CircularArc(c geom.Coord, v1 geom.Coord, v2 geom.Coord, r float64, s ...string) {
	p1 := v1.Minus(c).Unit().Times(r).Plus(c)
	p2 := v2.Minus(c).Unit().Times(r).Plus(c)
	a := geom.VertexAngle(v1, c, v2)
	largeArc := a > math.Pi
	sweep := a > 0
	svg.printf("<path d='M%f,%f A%f,%f 0 %s,%s %f,%f' %s/>\n",
		p1.X, p1.Y, r, r, onezero(largeArc), onezero(sweep), p2.X, p2.Y, extraparams(s))
}

////////////////////////////////////////////////////////////////////////////
// Penrose stuff

type CutLine struct {
	A, B geom.Coord
}

type MarkArc struct {
	C, P1, P2 geom.Coord
	R         float64
}

type RenderOutput struct {
	cuts  []CutLine
	mark1 []MarkArc
	mark2 []MarkArc
}

func (me *RenderOutput) AddCutLine(p1, p2 geom.Coord) {
	me.cuts = append(me.cuts, CutLine{p1, p2})
}

func (me *RenderOutput) AddMark1Arc(c, v1, v2 geom.Coord, r float64) {
	me.mark1 = append(me.mark1, MarkArc{c, v1, v2, r})
}
func (me *RenderOutput) AddMark2Arc(c, v1, v2 geom.Coord, r float64) {
	me.mark2 = append(me.mark2, MarkArc{c, v1, v2, r})
}

func (me *RenderOutput) MakeSVG(s *SVG) {
	for _, ma := range me.mark1 {
		s.CircularArc(ma.C, ma.P1, ma.P2, ma.R, MARK1_STYLE)
	}
	for _, ma := range me.mark2 {
		s.CircularArc(ma.C, ma.P1, ma.P2, ma.R+DOUBLE_STROKE_OFFSET, MARK2_STYLE)
		s.CircularArc(ma.C, ma.P1, ma.P2, ma.R-DOUBLE_STROKE_OFFSET, MARK2_STYLE)
	}
	for _, cl := range me.cuts {
		s.Line(cl.A, cl.B, CUT_STYLE)
	}
}

// Penrose primitives
type PenrosePrimitive interface {
	Render(ro *RenderOutput)
	Deflate() []PenrosePrimitive
}

const (
	c1 = math.Phi - 1.0
	c2 = 2.0 - math.Phi
)

type halfKite struct {
	*geom.Triangle
}

func (me halfKite) Render(ro *RenderOutput) {
	ro.AddCutLine(me.B, me.C)
	ro.AddCutLine(me.C, me.A)

	rB := me.A.Minus(me.C).Magnitude()
	rA := rB * 1.0 / math.Phi

	ro.AddMark2Arc(me.A, me.B, me.C, rA)
	ro.AddMark1Arc(me.B, me.C, me.A, rB)
}

func (me halfKite) Deflate() []PenrosePrimitive {
	r := make([]PenrosePrimitive, 3)
	d := me.A.Times(c1).Plus(me.B.Times(c2))
	e := me.B.Times(c1).Plus(me.C.Times(c2))
	r[0] = halfKite{&geom.Triangle{d, me.C, me.A}}
	r[1] = halfKite{&geom.Triangle{d, me.C, e}}
	r[2] = halfDart{&geom.Triangle{me.B, e, d}}
	return r
}

type halfDart struct {
	*geom.Triangle
}

func (me halfDart) Render(ro *RenderOutput) {
	ro.AddCutLine(me.B, me.C)
	ro.AddCutLine(me.C, me.A)

	r := me.A.Minus(me.C).Magnitude()
	rA := r / math.Pow(math.Phi, 2)
	rB := r / math.Pow(math.Phi, 3)

	ro.AddMark1Arc(me.A, me.B, me.C, rA)
	ro.AddMark2Arc(me.B, me.C, me.A, rB)
}

func (me halfDart) Deflate() []PenrosePrimitive {
	r := make([]PenrosePrimitive, 2)
	d := me.A.Times(c2).Plus(me.C.Times(c1))
	r[0] = halfDart{&geom.Triangle{me.C, d, me.B}}
	r[1] = halfKite{&geom.Triangle{me.B, me.A, d}}
	return r
}

func degToRads(d float64) float64 {
	return d * math.Pi / 180.0
}

// Starting shape "constants"
var HalfKite = halfKite{
	&geom.Triangle{
		geom.Coord{0, 0},
		geom.Coord{math.Phi * math.Cos(degToRads(72)), math.Phi * math.Sin(degToRads(72))},
		geom.Coord{1, 0},
	},
}

var HalfDart = halfDart{
	&geom.Triangle{
		geom.Coord{1, 0},
		geom.Coord{c1 * math.Cos(degToRads(36)), c1 * math.Sin(degToRads(36))},
		geom.Coord{0, 0},
	},
}

func Sun() []PenrosePrimitive {
	r := make([]PenrosePrimitive, 0, 10)
	for i := 0; i < 5; i++ {
		r = append(r, halfKite{
			&geom.Triangle{
				geom.Coord{math.Cos(degToRads(float64(72 * i))), math.Sin(degToRads(float64(72 * i)))},
				geom.Coord{0, 0},
				geom.Coord{math.Cos(degToRads(float64(36 + 72.0*i))), math.Sin(degToRads(float64(36 + 72.0*i)))},
			},
		})
		r = append(r, halfKite{
			&geom.Triangle{
				geom.Coord{math.Cos(degToRads(float64(72 * i))), math.Sin(degToRads(float64(72 * i)))},
				geom.Coord{0, 0},
				geom.Coord{math.Cos(degToRads(float64(-36 + 72.0*i))), math.Sin(degToRads(float64(-36 + 72.0*i)))},
			},
		})
	}
	return r
}

func BoundsOfPrimitiveSlice(ps []PenrosePrimitive) geom.Rect {
	bounds := ps[0].(geom.Bounded).Bounds()
	for _, p := range ps[1:] {
		bounds.ExpandToContainRect(p.(geom.Bounded).Bounds())
	}
	return bounds
}

////////////////////////////////////////////////////////////////////////////
func main() {
	shapes := Sun()
	// The laser cutter can cut 400x600.  Make our dimensions match that.
	bounds := geom.Rect{geom.Coord{-300, -200}, geom.Coord{300, 200}}
	bounds.Scale(1.0/350.0, 1.0/350.0)
	fmt.Fprintln(os.Stderr, bounds)

	// Deflate the shapes
	for i := 0; i < 5; i++ {
		newShapes := make([]PenrosePrimitive, 0, 3*len(shapes))
		for _, shape := range shapes {
			newShapes = append(newShapes, shape.Deflate()...)
		}
		shapes = newShapes
	}

	// Remove any shapes out of bounds
	newShapes := make([]PenrosePrimitive, 0, len(shapes))
	for _, shape := range shapes {
		if bounds.ContainsRect(shape.(geom.Bounded).Bounds()) {
			newShapes = append(newShapes, shape)
		}
	}
	shapes = newShapes

	// Render to drawing primitives
	ro := &RenderOutput{}
	for _, shape := range shapes {
		shape.Render(ro)
	}

	s := NewSVG(os.Stdout)
	s.Start(bounds, DEFAULT_STYLE)
	ro.MakeSVG(s)
	s.End()
}
