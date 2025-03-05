#!/opt/Wolfram/WolframScript/bin/wolframscript -format OutputForm

(* Psif := Psi[x1, y1, z1, x2, y2, z2] *)
Psif := Psi[r1, r2, r12]
r[x_, y_, z_] = Sqrt[x^2 + y^2 + z^2]
r1 = r[x1, y1, z1]
r2 = r[x2, y2, z2]
r12 = r[x2-x1, y2-y1, z2-z1]
eqn = - 1/2*Laplacian[Psif, {x1, y1, z1}] - 1/2*Laplacian[Psif, {x2, y2, z2}] - (2/r1 + 2/r2 - 1/r12) * Psif

(* Print[eqn] *)

(* Hylleraas coordinates
 * r1, theta1, psi1, r2, r12, 
 *)

(* Given the PDE "eqn"
 *
 * Change to r1, r2, r12
 *
 * Assumption: solution only depends on r1, r2, r12
 *
 * eqn[x1, y1, z1, x2, y2, z2]
 *
 * sol[r1, r2, r12]
 *
 * D[Psif, x1] = D[r1, x1] D[Psif, r1] + D[r2, x1] D[Psi, r2] + D[r12, x1] D[Psi, r12]
 *
 * r1 = sqrt(x1^2 + y1^2 + z1^2)
 *)

(* Print[eqn] *)

(* Suggested by https://mathematica.stackexchange.com/a/10460 *)
(* uu := u[x, y, z, t] *)
(* uu[x_, y_, z_, t_] := u[x, y, z, t] *)
(* Print[Dt[uu, t]] *)

r1 := Sqrt[x1^2 + y1^2 + z1^2]
r2 := Sqrt[x2^2 + y2^2 + z2^2]
r12 := Sqrt[(x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2]
uu := u[r1, r2, r12]
uu[r1_, r2_, r12_] := u[r1, r2, r12]

(* Now we want to set a bunch of partial derivatives.  We set TagSet (/:);
 * its use is suggested in the Dt documentation
 *
 * I originally had the partial derivatives listed out like this:
 * y1 /: Dt[y1, x1] = 0
 * z1 /: Dt[z1, x1] = 0
 * (you can't assign Dt[y1, x1] globally; you have to TagSet it on y1)
 *
 * But there are so many of them, I don't want to list them all out.
 * gpt-4o suggested the use of Scan, like this:
 *
 * Scan[(# /: Dt[#, x1] = 0) &, {y1, z1, x2, y2, z2}]
 *
 * gpt-4o also suggested the use of If, so I can include x1 in the list:
 * Scan[If[# =!= x1, (# /: Dt[#, x1] = 0)] &, {x1, y1, z1, x2, y2, z2}]
 *
 * gpt-4o suggested Subsets to form pairs and scan over the pairs:
 *
 * pairs = Tuples[vars, 2];
 * Scan[Function[{pair}, If[pair[[1]] =!= pair[[2]], Evaluate[pair[[1]]] /: Dt[pair[[1]], pair[[2]]] = 0]], pairs]
 *
 * Another gpt-4o prompt or two gets us to the final form (@@ is Apply and =!= is UnsameQ):
 *)

vars = {x1, y1, z1, x2, y2, z2};
pairs = Tuples[vars, 2];

TagSetDifferential[a_, b_] :=
  If[a =!= b,
    a /: Dt[a, b] = 0
  ]

Scan[TagSetDifferential @@ # &, pairs]

(* //. is ReplaceRepeated
 * This doesn't work:
 * Print[Dt[uu, x1] //. {r12 -> R12, r1 -> R1, r2 -> R2}]
 *
 * https://mathematica.stackexchange.com/a/84909 suggests using FullForm to diagnose this, like this:
 * Print[expr // FullForm]
 *
 * The reason it doesn't work is because some of the square roots appear to the -(1/2) power, not the (1/2) power
 *
 * So we include need to include negative powers (and third powers)
 *)

rules = {r12 -> R12, r1 -> R1, r2 -> R2, 1/r1 -> 1/R1, 1/r2 -> 1/R2, 1/r12 -> 1/R12, 1/r1^3 -> 1/R1^3, 1/r2^3 -> 1/R2^3, 1/r12^3 -> 1/R12^3}

LHS[f_] := f[[1]]

(* This doesn't do what I want:
 * eqn //. rules // Simplify //. x1^2 -> X1^2
 *
 * For operator precedence table, see https://reference.wolfram.com/language/tutorial/OperatorInputForms.html.en
 * //. is ReplaceRepeated; it has higher precedence than // (right apply)
 * So the expression nests like this:
 * (eqn //. rules) // (Simplify //. x1^2 -> X1^2)
 *
 * The right parenthesized expression just evaluates to Simplify:
 * Print[Simplify // x1^2 -> X1^2]
 * so the net effect is to drop the right hand rule completely.
 *)

eqnsimp = eqn //. rules // Simplify

(* Print[(eqn //. rules // Simplify) //. x1^2 -> X1^2] *)

Print[(eqn //. rules // Simplify)]
(* Print[(eqn //. rules // Simplify) /. {x1^2 -> (X1^2 - r12^2 + R12^2 // Expand)} /. {x1^2 -> X1^2}] *)
Print[(eqn //. rules // Simplify) /. {x1^2 -> (X1^2 - r12^2 + R12^2 // Expand)} /. {x1^2 -> X1^2} /. {x2^2 -> (X2^2 - r2^2 + R2^2 // Expand)} /. {X2^2 -> x2^2}]
(* Print[(eqn //. rules // Simplify) /. {x1^2 -> (X1^2 - r12^2 + R12^2 // Expand)} /. {x1^2 -> X1^2} /. {x2^2 -> (X2^2 - r2^2 + R2^2 // Expand)} /. {X2^2 -> x2^2} /. {x1*x2 -> X1*X2 - (1/2*(r1^2 + r2^2 - r12^2 // Expand)) + 1/2*(R1^2 + R2^2 - R12^2)}] *)
(* Print[(eqn //. rules // Simplify) /. {x1^2 -> (X1^2 - r12^2 + R12^2 // Expand)} /. {x1^2 -> X1^2} /. {x2^2 -> (X2^2 - r2^2 + R2^2 // Expand)} /. {X2^2 -> x2^2} /. {x1*x2 + y1*y2 + z1*z2 -> - 1/2*(R1^2 + R2^2 - R12^2)}] *)
bwb = (eqn //. rules // Simplify) /. {x1^2 -> (X1^2 - r12^2 + R12^2 // Expand)} /. {x1^2 -> X1^2} /. {x2^2 -> (X2^2 - r2^2 + R2^2 // Expand)} /. {X2^2 -> x2^2} //. {x1*x2 -> - 1/2*(R1^2 + R2^2 - R12^2) - y1*y2 - z1*z2}
Print[Distribute /@ (Simplify /@ bwb)]

(* Print[Psif + (R2^2*(-3*R12^3 + Expand[r12^2])*Derivative[0,0,1][Psi][r1, r2, r12])/(R12^2) //. rules //. x1^2 -> X1^2] *)
