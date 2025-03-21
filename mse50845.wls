#!/opt/Wolfram/WolframScript/bin/wolframscript -format OutputForm
(*
 * From Mathematica Stack Exchange #50845   https://mathematica.stackexchange.com/questions/50845
 *
 * Given the PDE: x^2 u_xx - u_yy + u_y = 0
 *
 * We want to change to s(x,y) = x e^y
 *                      t(x,y) = x e^-y
 *
 * Goal: convert the PDE to use u_ss, u_tt, u_ts, u_t, u_s
 *)

expr = x^2 D[u[x, y], {x, 2}] - D[u[x, y], {y, 2}] + D[u[x, y], y]
(* Print[expr] *)

(* This is how we figured out the expressions for x and y *)
(* Print[Eliminate[s == x Exp[y] && t == x Exp[-y], x]] *)
(* Print[Eliminate[s == x Exp[y] && t == x Exp[-y], y]] *)

(* Eliminate x and solve for y *)
(* Print[Solve[Eliminate[s == x Exp[y] && t == x Exp[-y], x], y]] *)

(* Eliminate y and solve for x *)
(* Print[Solve[Eliminate[s == x Exp[y] && t == x Exp[-y], y], x]] *)

(* This doesn't work; why? *)
(* Print[Solve[s == x Exp[y] && t == x Exp[-y], x]] *)
Print[Solve[x > 0 && y > 0 && s == x Exp[y] && t == x Exp[-y], {x,y}]]

$Assumptions = {s > 0, t > 0}
(* /. is ReplaceAll
 * -> is Rule
 * (body &) is a pure function
 * // is Postfix ( postfix function application: a // b is b[a] )
 * operator precedence table: https://reference.wolfram.com/language/tutorial/OperatorInputForms.html
 * operator precedence: Rule over ReplaceAll
 *                      Postfix isn't on that table, so I guess it's lower than everything
 * Rule binds tighter than ReplaceAll, and ReplaceAll groups a /. b /. c like ((a /. b) /. c)
 * expr /. u -> (body &) /. {} // Simplify
 * (expr /. (u -> (body &))) /. {} // Simplify
 *  Simplify[ (expr /. (u -> (body &))) /. {} ]
 *
 * expr is the PDE, which references u[x,y]
 * First we map u to its function form
 *)

(* This yields u[E^y*x, x/E^y] *)
(* u[x,y] /. u -> (u[# Exp[#2], # Exp[-#2]] &) *)
(* This yields Derivative[0, 1][u][E^y*x, x/E^y]/E^y + E^y*Derivative[1, 0][u][E^y*x, x/E^y] *)
(* D[u[x,y],x] /. u -> (u[# Exp[#2], # Exp[-#2]] &) *)
(* we use the PDE instead of u[x,y] or D[u[x,y],x], then replace x and y with their s/t forms, which converts u[E^y*x, x/E^y] to u[s,t] and simplify *)

bwb = expr /. u -> (u[# Exp[#2], # Exp[-#2]] &) /. {x -> Sqrt[s t],
                                                    y -> Log[Sqrt[s/t]]} // Simplify

Print[bwb]
