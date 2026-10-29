import Std

/-!+One-loop projected 3221 intersection with fixed minimal-survival coefficients.
A_i are corrected inverse couplings at MZ; x=ln(MI/MZ)/(2*pi),
y=ln(MG/MI)/(2*pi). Equations are multiplied by 30 to clear denominators.
This file proves algebra, not the beta coefficients or any claim about nature.
-/
namespace GreatPumpkin

-- Standard ordered-field laws are explicit parameters, not new axioms.
variable {R : Type} [Lean.Grind.Field R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [Lean.Grind.OrderedRing R]

def Meets (a b c x y u : R) : Prop :=
  30*a = 30*u + 123*x + 24*y ∧
  30*b = 30*u - 95*x - 90*y ∧
  30*c = 30*u - 210*x - 210*y

theorem invariant (a b c x y u : R) (h : Meets a b c x y u) :
    40*a - 78*b + 38*c = 145*x := by
  unfold Meets at h
  grind

theorem solution (a b c : R) :
    Meets a b c ((40*a-78*b+38*c)/145)
      ((-115*a+333*b-218*c)/435) ((35*a+693*b-293*c)/435) := by
  unfold Meets
  grind

theorem unique (a b c x y u x' y' u' : R)
    (h : Meets a b c x y u) (h' : Meets a b c x' y' u') :
    x = x' ∧ y = y' ∧ u = u' := by
  unfold Meets at *
  grind

theorem threshold_obstruction (a b c d1 d2 d3 x y u bound xmin : R)
    (h : Meets (a+d1) (b+d2) (c+d3) x y u)
    (hb : 40*d1 - 78*d2 + 38*d3 ≤ bound)
    (hx : xmin ≤ x)
    (gap : 40*a - 78*b + 38*c + bound < 145*xmin) : False := by
  have hi := invariant (a+d1) (b+d2) (c+d3) x y u h
  grind

theorem common_shift (a b c x y u k : R) (h : Meets a b c x y u) :
    Meets (a+k) (b+k) (c+k) x y (u+k) := by
  unfold Meets at *
  grind

-- Arbitrary independent offsets always manufacture a desired intersection.
theorem flexible_fit (a b c x y u : R) :
    ∃ d1 d2 d3, Meets (a+d1) (b+d2) (c+d3) x y u := by
  refine ⟨u + 41*x/10 + 4*y/5 - a,
    u - 19*x/6 - 3*y - b, u - 7*x - 7*y - c, ?_⟩
  unfold Meets
  grind

theorem positive_scales_iff (a b c : R) :
    (∃ x y u, Meets a b c x y u ∧ 0 < x ∧ 0 < y ∧ 1 < u) ↔
    (0 < (40*a-78*b+38*c)/145 ∧
     0 < (-115*a+333*b-218*c)/435 ∧
     1 < (35*a+693*b-293*c)/435) := by
  constructor
  · rintro ⟨x, y, u, h, hx, hy, hu⟩
    have hs := solution a b c
    have he := unique a b c x y u _ _ _ h hs
    grind
  · rintro ⟨hx, hy, hu⟩
    exact ⟨_, _, _, solution a b c, hx, hy, hu⟩

-- Non-vacuity: the standard library supplies the ordered-field instance.
example (a b c : Rat) :
    Meets a b c ((40*a-78*b+38*c)/145)
      ((-115*a+333*b-218*c)/435) ((35*a+693*b-293*c)/435) := solution a b c

def Within (z radius : R) : Prop := -radius ≤ z ∧ z ≤ radius

-- Coefficients are 2*w dot (P*b_H), in the documented heavy-field order.
-- Their derivation from representations is independently checked in Python.
theorem mass_box_support (z1 z2 z3 z4 z5 z6 z7 z8 z9 z10 z11 z12 z13 z14 z15 z16 radius : R)
    (h1 : Within z1 radius)
    (h2 : Within z2 radius)
    (h3 : Within z3 radius)
    (h4 : Within z4 radius)
    (h5 : Within z5 radius)
    (h6 : Within z6 radius)
    (h7 : Within z7 radius)
    (h8 : Within z8 radius)
    (h9 : Within z9 radius)
    (h10 : Within z10 radius)
    (h11 : Within z11 radius)
    (h12 : Within z12 radius)
    (h13 : Within z13 radius)
    (h14 : Within z14 radius)
    (h15 : Within z15 radius)
    (h16 : Within z16 radius) :
    9*z1 +9*z2 +9*z3 +9*z4 +207*z5 +75*z6 -201*z7 -129*z8 -28*z9 +8*z10 -18*z11 +14*z12 +14*z13 +19*z14 -26*z15 +8*z16 ≤ 783*radius := by
  unfold Within at *
  grind

omit [LE R] [LT R] [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [Lean.Grind.OrderedRing R] in
theorem hypercharge_hidden_direction (r b k : R) :
    3*(r+2*k) + 2*(b-3*k) = 3*r+2*b := by
  grind

#print axioms mass_box_support
#print axioms hypercharge_hidden_direction
#print axioms invariant
#print axioms solution
#print axioms unique
#print axioms threshold_obstruction
#print axioms common_shift
#print axioms flexible_fit
#print axioms positive_scales_iff
end GreatPumpkin
