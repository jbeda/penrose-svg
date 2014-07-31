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

// Penrose primitive
type PenrosePrimitive interface {
	Draw(s *SVG)
	Deflate() []PenrosePrimitive
}

const (
	c1 = math.Phi - 1.0
	c2 = 2.0 - math.Phi
)

type halfKite struct {
	*geom.Triangle
}

func (me halfKite) Draw(s *SVG) {
	s.Line(me.B, me.C, CUT_STYLE)
	s.Line(me.C, me.A, CUT_STYLE)

	rB := me.A.Minus(me.C).Magnitude()
	rA := rB * 1.0 / math.Phi

	s.CircularArc(me.A, me.B, me.C, rA+DOUBLE_STROKE_OFFSET, MARK2_STYLE)
	s.CircularArc(me.A, me.B, me.C, rA-DOUBLE_STROKE_OFFSET, MARK2_STYLE)
	s.CircularArc(me.B, me.C, me.A, rB, MARK1_STYLE)
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

func (me halfDart) Draw(s *SVG) {
	s.Line(me.B, me.C, CUT_STYLE)
	s.Line(me.C, me.A, CUT_STYLE)

	r := me.A.Minus(me.C).Magnitude()
	rA := r / math.Pow(math.Phi, 2)
	rB := r / math.Pow(math.Phi, 3)

	s.CircularArc(me.A, me.B, me.C, rA, MARK1_STYLE)
	s.CircularArc(me.B, me.C, me.A, rB+DOUBLE_STROKE_OFFSET, MARK2_STYLE)
	s.CircularArc(me.B, me.C, me.A, rB-DOUBLE_STROKE_OFFSET, MARK2_STYLE)
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
	bounds := BoundsOfPrimitiveSlice(shapes)

	for i := 0; i < 5; i++ {
		newShapes := make([]PenrosePrimitive, 0, 3*len(shapes))
		for _, shape := range shapes {
			newShapes = append(newShapes, shape.Deflate()...)
		}
		shapes = newShapes
	}

	s := NewSVG(os.Stdout)
	s.Start(bounds, DEFAULT_STYLE)
	for _, shape := range shapes {
		shape.Draw(s)
	}
	s.End()
}
