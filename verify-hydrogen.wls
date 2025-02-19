#!/opt/Wolfram/WolframScript/bin/wolframscript

Psif := Psi[x, y, z]
r[x_, y_, z_] = Sqrt[x^2 + y^2 + z^2]
eqn = - 1 / 2 * Laplacian[Psif, {x, y, z}] - 1 / r[x, y, z] * Psif == 0
sol[x_, y_, z_] := BesselJ[0, 2 * Sqrt[x + r[x, y, z]]]
Print[FullSimplify[eqn /. Psi -> ( {x, y, z} |-> sol[x, y, z] ) ]]
