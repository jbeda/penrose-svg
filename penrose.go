package main

import (
	"fmt"
	"github.com/jbeda/geom"
	"github.com/jbeda/geom/qtree"
	"io"
	"math"
	"os"
	"strings"
)

// Arc drawing functions/strategies
type ArcFuncType int

const (
	ARC_FUNC_CIRCULAR ArcFuncType = iota
	ARC_FUNC_BEZIER
)

// Tunable constants for output
const (
	DEFAULT_STYLE                = "stroke-width: 0.002; stroke-linecap: round; fill: none"
	DOUBLE_STROKE_OFFSET         = 0.002
	CUT_STYLE                    = "stroke: black"
	MARK1_STYLE                  = "stroke: green"
	MARK2_STYLE                  = "stroke: red"
	DEFLATE_LEVEL                = 5
	SQUISH_ARC_FUNC              = ARC_FUNC_CIRCULAR
	SQUISH_ARC_FACTOR            = 0.9
	SQUISH_ARC_BEZIER_ROUNDESS_A = 0.15
	SQUISH_ARC_BEZIER_ROUNDESS_B = 0.1
)

// Mathematical constants for generating Penrose decompositions.
const (
	C1 = math.Phi - 1.0
	C2 = 2.0 - math.Phi
)

////////////////////////////////////////////////////////////////////////////
// SVG serialization helper
type SVG struct {
	writer io.Writer
}

func NewSVG(w io.Writer) *SVG {
	return &SVG{w}
}

func (svg *SVG) printf(format string, a ...interface{}) (n int, errno error) {
	return fmt.Fprintf(svg.writer, format, a...)
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
	svg.printf("</svg>\n")
}

func (svg *SVG) Line(p1 geom.Coord, p2 geom.Coord, s ...string) {
	svg.printf("<line x1='%f' y1='%f' x2='%f' y2='%f' %s/>\n", p1.X, p1.Y, p2.X, p2.Y, extraparams(s))
}

func (svg *SVG) Circle(c geom.Coord, r float64, s ...string) {
	svg.printf("<circle cx='%f' cy='%f' r='%f' %s/>\n", c.X, c.Y, r, extraparams(s))
}

func (svg *SVG) CircularArc(p1, p2 geom.Coord, r float64, largeArc, sweep bool, s ...string) {
	svg.printf("<path d='M%f,%f A%f,%f 0 %s,%s %f,%f' %s/>\n",
		p1.X, p1.Y, r, r, onezero(largeArc), onezero(sweep), p2.X, p2.Y, extraparams(s))
}

func (svg *SVG) QuadBezier(p1 geom.Coord, ctrl1 geom.Coord, p2 geom.Coord, s ...string) {
	svg.printf("<path d='M%f,%f Q%f,%f %f,%f' %s/>\n",
		p1.X, p1.Y, ctrl1.X, ctrl1.Y, p2.X, p2.Y, extraparams(s))
}

func (svg *SVG) CubicBezier(p1, ctrl1, ctrl2, p2 geom.Coord, s ...string) {
	svg.printf("<path d='M%f,%f C%f,%f %f,%f %f,%f' %s/>\n",
		p1.X, p1.Y, ctrl1.X, ctrl1.Y, ctrl2.X, ctrl2.Y, p2.X, p2.Y, extraparams(s))
}

////////////////////////////////////////////////////////////////////////////
// Decomposed Penrose Geometry
//
// This stuff deals with pre-styled lines and arcs and is tailored to Penrose
// tiling needs.  This could be a generic retained mode 2D model thing but it
// is more quick and dirty for now.

// Comparing floating point sucks.  This is probably wrong in the general case
// but is good enough for this application.  Don't assume I know what I'm
// doing here.  I'm pulling stuff out of my butt.
const FLOAT_EQUAL_THRESH = 0.000000001

func FloatAlmostEqual(a, b float64) bool {
	return math.Abs(a-b) < FLOAT_EQUAL_THRESH
}

func AlmostEqualsCoord(a, b geom.Coord) bool {
	return FloatAlmostEqual(a.X, b.X) && FloatAlmostEqual(a.Y, b.Y)
}

// +++ CutLine
type CutLine struct {
	A, B geom.Coord
}

func AlmostEqualsCutLines(a, b CutLine) bool {
	return (AlmostEqualsCoord(a.A, b.A) && AlmostEqualsCoord(a.B, b.B)) ||
		(AlmostEqualsCoord(a.A, b.B) && AlmostEqualsCoord(a.B, b.A))
}

func (cl CutLine) Equals(oi interface{}) bool {
	ocl, ok := oi.(CutLine)
	return ok && AlmostEqualsCutLines(cl, ocl)
}

func (cl CutLine) Bounds() geom.Rect {
	r := geom.Rect{cl.A, cl.A}
	r.ExpandToContainCoord(cl.B)
	return r
}

// +++ MarkArc
type MarkArc struct {
	C, P1, P2 geom.Coord
	R         float64
}

func AlmostEqualMarkArcs(a, b MarkArc) bool {
	return AlmostEqualsCoord(a.C, b.C) && FloatAlmostEqual(a.R, b.R) &&
		((AlmostEqualsCoord(a.P1, b.P1) && AlmostEqualsCoord(a.P2, b.P2)) ||
			(AlmostEqualsCoord(a.P1, b.P2) && AlmostEqualsCoord(a.P2, b.P1)))
}

func (ma MarkArc) Equals(oi interface{}) bool {
	oma, ok := oi.(MarkArc)
	return ok && AlmostEqualMarkArcs(ma, oma)
}

// BUGBUG: This isn't the true bounds but is good enough for finding the
// endpoints.  To realy get bounds we'd need to look for arcs that cross an
// axis and add that point to the bounding rect.
func (ma MarkArc) Bounds() geom.Rect {
	p1 := ma.P1.Minus(ma.C).Unit().Times(ma.R).Plus(ma.C)
	p2 := ma.P2.Minus(ma.C).Unit().Times(ma.R).Plus(ma.C)
	r := geom.Rect{p1, p1}
	r.ExpandToContainCoord(p2)
	return r
}

// +++ RenderOutput
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

func (me *RenderOutput) RemoveDuplicates() {
	// First calculate the bounds of all the cut lines
	allBounds := geom.NilRect()
	for _, cl := range me.cuts {
		allBounds.ExpandToContainRect(cl.Bounds())
	}

	fmt.Fprintf(os.Stderr, "Number of cut lines before RemoveDuplicates: %d\n", len(me.cuts))

	// Create our qtree based on that
	qt := qtree.New(qtree.ConfigDefault(), allBounds)

	// Add all cuts if they aren't already there
	for _, cl := range me.cuts {
		qt.FindOrInsert(cl)
	}

	// Extract out all elements
	col := make(map[qtree.Item]bool)
	qt.Enumerate(col)
	newCuts := []CutLine{}
	for item, _ := range col {
		newCuts = append(newCuts, item.(CutLine))
	}
	me.cuts = newCuts
	fmt.Fprintf(os.Stderr, "Number of cut lines after RemoveDuplicates: %d\n", len(me.cuts))
}

func (me *RenderOutput) MakeSVG(s *SVG) {
	squishFunc := SquishedArcCircle
	if SQUISH_ARC_FUNC == ARC_FUNC_BEZIER {
		squishFunc = SquishedArcBezier
	}

	for _, ma := range me.mark1 {
		squishFunc(s, &ma, 0, MARK1_STYLE)
	}
	for _, ma := range me.mark2 {
		squishFunc(s, &ma, +DOUBLE_STROKE_OFFSET, MARK2_STYLE)
		squishFunc(s, &ma, -DOUBLE_STROKE_OFFSET, MARK2_STYLE)
	}
	for _, cl := range me.cuts {
		s.Line(cl.A, cl.B, CUT_STYLE)
	}
}

func SquishedArcCircle(svg *SVG, ma *MarkArc, rOffset float64, s ...string) {
	p1 := ma.P1.Minus(ma.C).Unit().Times(ma.R*SQUISH_ARC_FACTOR + rOffset).Plus(ma.C)
	p2 := ma.P2.Minus(ma.C).Unit().Times(ma.R + rOffset).Plus(ma.C)
	a := geom.VertexAngle(ma.P1, ma.C, ma.P2)
	largeArc := a > math.Pi
	sweep := a > 0
	svg.CircularArc(p1, p2, ma.R+rOffset/2, largeArc, sweep, s...)
}

func SquishedArcBezier(svg *SVG, ma *MarkArc, rOffset float64, s ...string) {
	r := ma.R
	v1 := ma.P1.Minus(ma.C).Unit()
	v2 := ma.P2.Minus(ma.C).Unit()
	p1 := v1.Times(r*SQUISH_ARC_FACTOR + rOffset).Plus(ma.C)
	p2 := v2.Times(r + rOffset).Plus(ma.C)
	a := geom.VertexAngle(ma.P1, ma.C, ma.P2)

	var v1p, v2p geom.Coord
	if a > 0 {
		v1p = geom.Coord{-v1.Y, v1.X}.Unit()
		v2p = geom.Coord{v2.Y, -v2.X}.Unit()
	} else {
		v1p = geom.Coord{v1.Y, -v1.X}.Unit()
		v2p = geom.Coord{-v2.Y, v2.X}.Unit()
	}
	ctrlDist1 := SQUISH_ARC_BEZIER_ROUNDESS_A * math.Abs(a) * math.Pow(C1, DEFLATE_LEVEL)
	ctrl1 := v1p.Times(ctrlDist1).Plus(p1)
	ctrlDist2 := SQUISH_ARC_BEZIER_ROUNDESS_B*math.Abs(a)*math.Pow(C1, DEFLATE_LEVEL) + rOffset
	ctrl2 := v2p.Times(ctrlDist2).Plus(p2)
	// svg.Circle(p1, 0.002, s...)
	// svg.Circle(ctrl1, 0.001, s...)
	svg.CubicBezier(p1, ctrl1, ctrl2, p2, s...)
}

////////////////////////////////////////////////////////////////////////////
// Primary Penrose tile generation
type PenrosePrimitive interface {
	Render(ro *RenderOutput)
	Deflate() []PenrosePrimitive
}

// +++ halfKite
type halfKite struct {
	*geom.Triangle
}

func (me halfKite) Render(ro *RenderOutput) {
	ro.AddCutLine(me.B, me.C)
	ro.AddCutLine(me.C, me.A)

	rB := me.A.Minus(me.C).Magnitude()
	rA := rB * 1.0 / math.Phi

	ro.AddMark2Arc(me.A, me.B, me.C, rA)
	ro.AddMark1Arc(me.B, me.A, me.C, rB)
}

func (me halfKite) Deflate() []PenrosePrimitive {
	r := make([]PenrosePrimitive, 3)
	d := me.A.Times(C1).Plus(me.B.Times(C2))
	e := me.B.Times(C1).Plus(me.C.Times(C2))
	r[0] = halfKite{&geom.Triangle{d, me.C, me.A}}
	r[1] = halfKite{&geom.Triangle{d, me.C, e}}
	r[2] = halfDart{&geom.Triangle{me.B, e, d}}
	return r
}

// +++ halfDart
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
	ro.AddMark2Arc(me.B, me.A, me.C, rB)
}

func (me halfDart) Deflate() []PenrosePrimitive {
	r := make([]PenrosePrimitive, 2)
	d := me.A.Times(C2).Plus(me.C.Times(C1))
	r[0] = halfDart{&geom.Triangle{me.C, d, me.B}}
	r[1] = halfKite{&geom.Triangle{me.B, me.A, d}}
	return r
}

// +++ Starting Shapes
func degToRads(d float64) float64 {
	return d * math.Pi / 180.0
}

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
		geom.Coord{C1 * math.Cos(degToRads(36)), C1 * math.Sin(degToRads(36))},
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

func DeflatePenrosePrimitives(ps []PenrosePrimitive, levels int) []PenrosePrimitive {
	r := ps
	fmt.Fprintf(os.Stderr, "Starting primitive count: %d\n", len(r))
	for i := 0; i < levels; i++ {
		rNext := make([]PenrosePrimitive, 0, 3*len(r))
		for _, shape := range r {
			rNext = append(rNext, shape.Deflate()...)
		}
		r = rNext
		fmt.Fprintf(os.Stderr, "Primitive count after iteration %d: %d\n", i+1, len(r))
	}
	return r
}

func CullShapes(ps []PenrosePrimitive, bounds *geom.Rect) []PenrosePrimitive {
	r := make([]PenrosePrimitive, 0, len(ps))
	for _, shape := range ps {
		if bounds.ContainsRect(shape.(geom.Bounded).Bounds()) {
			r = append(r, shape)
		}
	}
	fmt.Fprintf(os.Stderr, "Primitive count after cull: %d\n", len(r))
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

	// Deflate the shapes
	shapes = DeflatePenrosePrimitives(shapes, DEFLATE_LEVEL)

	// Remove any shapes out of bounds
	shapes = CullShapes(shapes, &bounds)

	// Render to drawing primitives
	ro := &RenderOutput{}
	for _, shape := range shapes {
		shape.Render(ro)
	}
	ro.RemoveDuplicates()

	s := NewSVG(os.Stdout)
	s.Start(bounds, DEFAULT_STYLE)
	ro.MakeSVG(s)
	s.End()
}
