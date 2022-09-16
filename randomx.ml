(* #use "topfind"
#require "owl-top"
 *)

open Owl
open Gamma


(* launch on the go *)
(* ocaml str.cma gamma.cmo randomx.ml *)

(* launch
ocamlfind ocamlopt -o randomx -linkpkg -package owl gamma.ml randomx.ml
ocamlfind ocamlopt -o program -linkpkg -package pkg module1.cmx module2.cmx *)





(* input seed *)
let read_ints () =
    let line = read_line () in
    let ints = Str.split (Str.regexp "  *") line in
    List.map int_of_string ints;;

let seeds = read_ints () in let seed = List.hd seeds in Random.init(seed);;




(* simulation of rv of interest *)

module Simulate = struct 

    (** simulation of N(loc, scale) **)

    let normal = fun ~(loc: float) ~(scale: float) -> Stats.gaussian_rvs loc scale;;

    let normal_pdf = fun ~(loc: float) ~(scale: float) -> Stats.gaussian_pdf ~mu:loc ~sigma:scale;;


    (** simulation of InvGamma(loc, scale) **)

    let inv_gamma = fun ~(shape: float) ~(scale: float) -> 1. /. (Stats.gamma_rvs shape scale);;



    let inv_gamma_pdf = fun (x: float) ~(shape: float) ~(scale: float) ->
        ( 1. /. Lanczos.gamma shape ) *. scale ** shape *. x ** (-.shape -. 1.) *. exp (-.scale /. x)
        

end


(* 

let symetric_uniform = fun () -> 
    let u = Random.float 2. -. 1.
    in 
    u;;


let twice_standard_normal = fun () -> 
    let rec aux_normal (u_1: float) (u_2: float) = match u_1, u_2 with
        |u_1, u_2 when u_1**2. +. u_2**2. <= 1. -> u_1, u_2
        |_,_                                  -> aux_normal (symetric_uniform () ) (symetric_uniform ())
    in
    let u_1, u_2 = aux_normal (symetric_uniform ()) (symetric_uniform ())
    in
    let s = u_1**2. +. u_2**2.
    in
    let temp = sqrt( -2. *. log(s) /. s)
    in
    let x = u_1 *. temp and y = u_2 *. temp
    in
    [|x; y|];;


let normal = fun (loc: float) (scale: float) -> scale *. (twice_standard_normal ()).(0) +. loc;;
 *)










(* personal implementation of inv gamma *)
    (** simulation of InvGamma(loc, scale) **)

    (* let inv_gamma_pdf = fun (x: float) (loc: float) (scale: float) ->
        scale ** loc *. x ** (-loc -. 1.) *. exp (-scale /. x) /. Lanczos.gamma loc 

    let inv_gamma = fun (loc: float) (scale: float) -> 
    let partial_pdf = fun (x: float) -> inv_gamma_pdf x loc scale 
    in
    let m = b ** a / a in let x = pareto a in let u = Random.float 1.
    in
    let rec inv_gamma_aux x u = match x u with 
        |x, u when u *. m *. pareto_pdf x a > inv_gamma_pdf x a b ->  inv_gamma_aux (pareto a) (Random.float 1.)
        |x, _                                                     -> x
    in inv_gamma_aux x u;; *)
