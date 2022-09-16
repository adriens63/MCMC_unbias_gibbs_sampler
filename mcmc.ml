open Randomx
open Owl




(* launch *)
(*
ocamlfind ocamlopt -o mcmc -linkpkg -package owl gamma.ml randomx.ml mcmc.ml
ocamlfind ocamlopt -o mcmc -linkpkg -package owl gamma.cmx randomx.cmx mcmc.ml
ocamlfind ocamlopt -o program -linkpkg -package pkg module1.cmx module2.cmx
*)


(* simple gibbs sampling *)
type bayesianVariables = {capital_a: float; mu: float; theta: float array};;
type mcmcChain = bayesianVariables list;;
type sequentialVariable = float list;;


(* constants *)
let capital_v = 0.00434 and a = 1. and b = 2. ;;
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


let gibbs_sampler = fun (x_t_1: bayesianVariables) ->
    (* input is x_t_1, and this fun returns x_t *)
    (* extract the useful information from x_t-1, ie theta in our case *)
    let theta = x_t_1.theta
    in
    let capital_k_int = Array.length theta 
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
    let new_theta = Array.make capital_k_int 0.
    in
    for i = 0 to capital_k_int-1 do 
        new_theta.(i) <- Simulate.normal (inv_v_plus_a *. (mu *. capital_v +. z.(i) *. capital_a)) (inv_v_plus_a *. capital_a *. capital_v);
    done;
    
    (* storage of x_t *)
    let x_t = {capital_a = capital_a; mu = mu; theta = new_theta}
    in x_t;;




let rec generate_chain (n_steps: int) (x_t: bayesianVariables) = match n_steps with
    |n when n <= 0 -> []
    |_             -> x_t::generate_chain (n_steps - 1) (gibbs_sampler x_t)
;;

    
let chain_to_lists = fun (mcmc_chain: mcmcChain) -> 
    let rec acc (mcmc_chain: mcmcChain) (capital_as: sequentialVariable) (mus: sequentialVariable) (thetas_0: sequentialVariable) = match mcmc_chain with
    |[] -> capital_as, mus, thetas_0
    |e::q -> acc q (e.capital_a::capital_as) (e.mu::mus) (e.theta.(0)::thetas_0)
    in let (capital_as, mus, thetas_0) = acc mcmc_chain [] [] []
    in capital_as, mus, thetas_0
;;

let chain_to_matrix = fun (mcmc_chain: mcmcChain) (n_steps: int) ->
    let m = Mat.zeros n_steps 3
    in 
    let rec aux (mcmc_chain: mcmcChain) (i: int) = match mcmc_chain with
    |[] -> ()
    |e::q ->
                    Mat.set m i 0 e.capital_a;
                    Mat.set m i 1 e.mu;
                    Mat.set m i 2 e.theta.(0);
                    aux q (i+1);
                
    in 
    aux mcmc_chain 0;
    m;;






type coupledBayesianVariables = {x_t: bayesianVariables; y_t_1: bayesianVariables};;
type coupledMcmcChain = coupledBayesianVariables list;;

module CoupledGibbsSampler = struct

    let maximal_coupling = fun p_sampler p_pdf q_sampler q_pdf ->
        let x = p_sampler ()
        in
        let w = Stats.uniform_rvs 0. (p_pdf x)
        in 
        if w <= q_pdf x
        then x, x
        else let y_star = q_sampler ()
        in 
        let w_star = Stats.uniform_rvs 0. (q_pdf y_star)
        in
        let rec aux w_star y_star = match w_star, y_star with
        |w_star,y_star when w_star <= p_pdf y_star  -> begin
                                                        let new_y_star = q_sampler ()
                                                        in 
                                                        aux (Stats.uniform_rvs 0. (q_pdf new_y_star)) new_y_star
                                                        end
        |_, y_star                                  -> x, y_star
        in aux w_star y_star;;


    let coupled_normal = fun loc_1 scale_1 loc_2 scale_2 ->
        let p_pdf = Simulate.normal_pdf ~loc:loc_1 ~scale:scale_1
        in
        let p_sampler = fun () -> Simulate.normal ~loc:loc_1 ~scale:scale_1
        in
        let q_pdf = Simulate.normal_pdf ~loc:loc_2 ~scale:scale_2
        in
        let q_sampler = fun () -> Simulate.normal ~loc:loc_2 ~scale:scale_2
        in
        maximal_coupling p_sampler p_pdf q_sampler q_pdf;;


    let coupled_gamma = fun shape_1 scale_1 shape_2 scale_2 ->
        let p_pdf = Simulate.inv_gamma_pdf ~shape:shape_1 ~scale:scale_1
        in
        let p_sampler = fun () -> Simulate.inv_gamma ~shape:shape_1 ~scale:scale_1
        in
        let q_pdf = Simulate.inv_gamma_pdf ~shape:shape_2 ~scale:scale_2
        in
        let q_sampler = fun () -> Simulate.inv_gamma ~shape:shape_2 ~scale:scale_2
        in
        maximal_coupling p_sampler p_pdf q_sampler q_pdf;;

    let coupled_gibbs_sampler = fun (xy: coupledBayesianVariables) ->
        (* unpack *)
        let x_t = xy.x_t
        and y_t_1 = xy.y_t_1
        in
        
        (* for x_t *)
        let theta_x = x_t.theta
        in
        let capital_k_int = Array.length theta_x 
        in  
        let capital_k = float_of_int(capital_k_int)
        in
        let (sum_gamma_x, mean_mu_x) = sum_over_theta theta_x
        in

        (* for y_t_1 *)
        let theta_y = y_t_1.theta
        in
        let (sum_gamma_y, mean_mu_y) = sum_over_theta theta_y
        in


        (* simulation of A *)
        let (capital_a_x, capital_a_y) = coupled_gamma (a +. (capital_k -. 1.) /. 2.) (b +. sum_gamma_x) (a +. (capital_k -. 1.) /. 2.) (b +. sum_gamma_y)
        in
        (* simulation of mu *)
        let (mu_x, mu_y) = coupled_normal mean_mu_x (capital_a_x /. capital_k) mean_mu_x (capital_a_x /. capital_k)
        in
        (* computation of the constants *)
        let inv_v_plus_a_x = 1. /. (capital_v +. capital_a_x)
        in
        let inv_v_plus_a_y = 1. /. (capital_v +. capital_a_y)
        in
        (* simulation of theta *)
        let new_theta_x = Array.make capital_k_int 0.
        and new_theta_y = Array.make capital_k_int 0.
        in
        for i = 0 to capital_k_int-1 do 
            let (temp_x, temp_y) = coupled_normal (inv_v_plus_a_x *. (mu_x *. capital_v +. z.(i) *. capital_a_x)) (inv_v_plus_a_x *. capital_a_x *. capital_v) (inv_v_plus_a_y *. (mu_y *. capital_v +. z.(i) *. capital_a_y)) (inv_v_plus_a_y *. capital_a_y *. capital_v) in
            new_theta_x.(i) <- temp_x;
            new_theta_y.(i) <- temp_y;
        done;
        
        (* storage of x_t *)
        let x_t_plus_1 = {capital_a = capital_a_x; mu = mu_x; theta = new_theta_x} 
        and y_t = {capital_a = capital_a_y; mu = mu_y; theta = new_theta_y} 
        in
        let new_xy = {x_t = x_t_plus_1; y_t_1 = y_t}
        in new_xy;;
    


    let generate_estimator = fun (burnin: int) (m: int) (xy_0: coupledBayesianVariables) (h: bayesianVariables -> float) ->
        
        let rec generate_coupled_chain (tau: int) (xy_t: coupledBayesianVariables) = match tau, xy_t with
        |0, xy_t                                                       -> xy_t::generate_coupled_chain (tau + 1) (coupled_gibbs_sampler xy_t)
        |1, xy_t                                                       -> xy_t::generate_coupled_chain (tau + 1) (coupled_gibbs_sampler xy_t)        
        |tau, xy_t when ((tau < m) && (xy_t.x_t <> xy_t.y_t_1)) -> xy_t::generate_coupled_chain (tau + 1) (coupled_gibbs_sampler xy_t)
        |_, xy_t -> [xy_t]
        in
        let coupled_chain = generate_coupled_chain 0 xy_0
        in
        let rec _sum_differences (coupled_chain: coupledMcmcChain) = match coupled_chain with
            |[] -> 0.
            |e::q -> h(e.x_t) -. h(e.y_t_1) +. _sum_differences q
        in
        let rec generate_list_estimators (acc: int) (coupled_chain: coupledMcmcChain) = match acc, coupled_chain with
            |_, [] -> []
            |acc, e::q when acc < burnin ->  generate_list_estimators (acc+1) q
            |acc, e::q ->  (h(e.x_t) +. _sum_differences q)::generate_list_estimators (acc+1) q
        in
        let estimators = generate_list_estimators 0 coupled_chain
        in
        let mean = fun l ->
            let rec aux (l: float list) (sum: float) (len: float) = match l with
                |[] -> sum /. len
                |e::q -> aux q (sum +. e) (len +. 1.)
            in aux l 0. 0.
        in
        mean estimators;;
    

    let rec generate_estimators (burnin: int) (m: int) (xy_0: coupledBayesianVariables) (n_steps: int) = match n_steps with
        |n_steps when n_steps <= 0 -> []
        |_             -> {capital_a = generate_estimator burnin m xy_0 (fun z -> z.capital_a); mu = generate_estimator burnin m xy_0 (fun z -> z.mu); theta = Array.make 1 (generate_estimator burnin m xy_0 (fun z -> z.theta.(0))) }::generate_estimators burnin m xy_0 (n_steps - 1) 


end
