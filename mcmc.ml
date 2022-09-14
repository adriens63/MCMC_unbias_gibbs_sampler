open Randomx





(* launch *)
(*
ocamlfind ocamlopt -o mcmc -linkpkg -package owl gamma.ml randomx.ml mcmc.ml
ocamlfind ocamlopt -o mcmc -linkpkg -package owl gamma.cmx randomx.cmx mcmc.ml
ocamlfind ocamlopt -o program -linkpkg -package pkg module1.cmx module2.cmx
*)


(* simple gibbs sampling *)




(* constants *)
let capital_v = 0.00434 and a = -1. and b = 2. ;;
let z = [|0.395; 0.375; 0.355; 0.334; 0.313; 0.313; 0.291; 0.269; 0.247; 0.247; 0.224; 0.224; 0.224; 0.224; 0.224; 0.200; 0.175; 0.148|];;


(* we use another formula for the sum that depends on theta *)

let sum_over_theta = fun (theta: float array) -> 
    let capital_k = Array.length(theta) 
    in
    let capital_k_float = float_of_int(capital_k)
    in 
    let sum_of_squares = ref 0. and sum = ref 0.
    in
    for i = 0 to capital_k - 1 do
        sum_of_squares := !sum_of_squares +. theta.(i) ** 2.;
        sum := !sum +. theta.(i);
    done;
    (!sum_of_squares -. !sum**2. /. capital_k_float) /. 2. , !sum /. capital_k_float;;


exception ShapeMissmatch
let assign_end = fun (arr_1: 'a array) (arr_2: 'a array) ->
    let len_1 = Array.length arr_1 and len_2 = Array.length arr_2 
    in
    if len_2 > len_1 then raise ShapeMissmatch 
    else
    let start = len_1 - len_2
    in
    for i = start to len_1 do
        arr_1.(i) <- arr_2.(i - start);
    done;
    arr_1 ;;


let gibbs_sampler = fun (x: float array) ->
    (* x is in fact x_t_1, and this fun returns x_t *)
    (* extract the useful information from x_t-1, ie theta in our case *)
    let capital_k_int = Array.length x - 2 in 
    let theta = Array.sub x 2 capital_k_int
    in
    let capital_k = float_of_int(capital_k_int)
    in
    let (sum_gamma, mean_mu) = sum_over_theta theta
    in

    (* simulation of A *)
    let capital_a = Simulate.inv_gamma (a +. (capital_k -. 1.) /. 2.) (b +. sum_gamma)
    in
    (* simulation of mu *)
    let mu = Simulate.normal mean_mu (capital_a /. capital_k)
    in
    (* computation of the constants *)
    let inv_v_plus_a = 1. /. (capital_v +. capital_a)
    in
    (* simulation of theta *)
    for i = 0 to capital_k_int do 
        theta.(i) <- Simulate.normal (inv_v_plus_a *. (mu *. capital_v +. z.(i) *. capital_a)) (inv_v_plus_a *. capital_a *. capital_v);
    done;
    
    (* storage of x_t *)
    let x = assign_end x theta in
    x.(0) <- capital_a;
    x.(1) <- mu;
    x;;






    

(* test *)















