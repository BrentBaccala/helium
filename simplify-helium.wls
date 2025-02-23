#!/opt/Wolfram/WolframScript/bin/wolframscript -format OutputForm

Psif := Psi[x1, y1, z1, x2, y2, z2]
r[x_, y_, z_] = Sqrt[x^2 + y^2 + z^2]
eqn = - 1/2*Laplacian[Psif, {x1, y1, z1}] - 1/2*Laplacian[Psif, {x2, y2, z2}] - (2/r[x1, y1, z1] + 2/r[x2, y2, z2] - 1/r[x2-x1, y2-y1, z2-z1]) * Psif == 0

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

(* /: is TagSet; its use is suggested in the Dt documentation *)
(* I originally had them listed out like this: *)
(* y1 /: Dt[y1, x1] = 0 *)
(* z1 /: Dt[z1, x1] = 0 *)
(* gpt-4o suggested the use of Scan, like this: *)
(* Scan[(# /: Dt[#, x1] = 0) &, {y1, z1, x2, y2, z2}] *)
(* gpt-4o also suggested the use of If *)
(* Scan[If[# =!= x1, (# /: Dt[#, x1] = 0)] &, {x1, y1, z1, x2, y2, z2}] *)
(* gpt-4o suggested Subsets to form pairs: *)

vars = {x1, y1, z1, x2, y2, z2};

(* Form all pairs *)
pairs = Tuples[vars, 2];

(* Scan over the pairs *)
Scan[Function[{pair}, If[pair[[1]] =!= pair[[2]], Evaluate[pair[[1]]] /: Dt[pair[[1]], pair[[2]]] = 0]], pairs]

(* Scan[If[# =!= x1, (# /: Dt[#, x1] = 0)] &, {x1, y1, z1, x2, y2, z2}] *)
(* Scan[Function[var, If[var =!= x1, (var /: Dt[var, x1] = 0)]], {x1, y1, z1, x2, y2, z2}] *)
(* Scan[Function[{pair}, If[pair[[1]] =!= x1, Evaluate[pair[[1]]] /: Dt[pair[[1]], x1] = 0]], pairs] *)

(* //. is ReplaceRepeated *)
(* This doesn't work: *)
(* Print[Dt[uu, x1] //. {r12 -> R12, r1 -> R1, r2 -> R2} // FullForm] *)
(* https://mathematica.stackexchange.com/a/84909 suggests using FullForm to diagnose this, like this: *)
(* Print[expr // FullForm] *)
(* The reason it doesn't work is because some of the square roots appear to the -(1/2) power, not the (1/2) power *)

(* So do this: *)
bwb = Dt[uu,x1]
rules = {r12 -> R12, r1 -> R1, r2 -> R2, 1/r1 -> 1/R1, 1/r2 -> 1/R2, 1/r12 -> 1/R12}
Print[bwb //. rules]
