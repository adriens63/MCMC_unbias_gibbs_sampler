open Mcmc


(* main *)


(* launch *)
(* ocaml str.cma gamma.cmo randomx.ml *)




let main = fun () -> 

  let capital_k_int = Array.length z
  in
  let n_steps = read_int () 
  and 
  let x_0 = {capital_a = 1.; mu = 1.; theta = Array.make capital_k_int (snd sum_over_theta z)}
  in  
  let mcmc_chain = generate_chain n_steps x_0

  ;;

let () = main ()



