open Mcmc
open Randomx
open Owl_plplot
open Owl
open Format
(* main *)


(* launch *)
(* ocamlfind ocamlopt -o launch_mcmc -linkpkg -package owl,owl-plplot gamma.ml randomx.ml mcmc.ml launch_mcmc.ml *)
(* QT_QPA_PLATFORM=offscreen ./launch_mcmc  *)


(* the first input is the random seed *)
(* the second input is the number of samples *)


let main = fun () -> 

  let capital_k_int = Array.length z
  in
  let numbers = read_ints () 
  in 
  let n_steps = List.hd numbers 
  in 
  let x_0 = {capital_a = 1.; mu = 1.; theta = Array.make capital_k_int (snd (sum_over_theta z))}
  in  
  let mcmc_chain = generate_chain n_steps x_0
  in
  let m = chain_to_matrix mcmc_chain n_steps
  in
  let h = Plot.create ~m:2 ~n:2 "plot_003.png" in
  Plot.set_background_color h 255 255 255;

  (* focus on the subplot at 0,0 *)
  Plot.subplot h 0 0;
  Plot.set_title h "Distribution of A";
  Plot.set_xlabel h "y";
  Plot.set_ylabel h "Frequency";
  Plot.histogram ~h ~bin:50 (Mat.col m 0);

  (* focus on the subplot at 0,1 *)
  Plot.subplot h 0 1;
  Plot.set_title h "Distribution of mu";
  Plot.set_xlabel h "y";
  Plot.set_ylabel h "Frequency";
  Plot.histogram ~h ~bin:50 (Mat.col m 1);

  (* focus on the subplot at 1,0 *)
  Plot.subplot h 1 0;
  Plot.set_title h "Distribution of theta_0";
  Plot.set_xlabel h "y";
  Plot.set_ylabel h "Frequency";
  Plot.histogram ~h ~bin:50 (Mat.col m 2);

  (* focus on the subplot at 1,1 *)
  Plot.subplot h 1 1;
  Plot.set_foreground_color h 0 50 255;
  Plot.set_title h "Sine function";
  Plot.(plot_fun ~h ~spec:[ LineStyle 2 ] Maths.sin 0. 28.);
  Plot.autocorr ~h (Mat.sequential 1 28);

  (* output your final plot *)
  Plot.output h;;



let () = main ()



