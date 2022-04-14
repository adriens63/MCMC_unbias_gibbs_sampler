open Gamma

(*1*)


(* constants *)
let capital_v = 0.00434 and a = -1. and b = 2. ;;
let z = [|0.395; 0.375; 0.355; 0.334; 0.313; 0.313; 0.291; 0.269; 0.247; 0.247; 0.224; 0.224;
0.224; 0.224; 0.224; 0.200; 0.175; 0.148|];;

(* simulation d'une N(0.,1.) par Box-Muller amélioré (cf td1) *)

(** simulation d'une U[-1., 1.] **)

let symetric_uniform = fun () -> 
    let u = Random.float 2. -. 1.
    in 
    u;;

(** simulation d'une N(0.,1.) **)

let standard_normal = fun () -> 
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
    x, y;;

let normal = fun (mu: float) (sigma: float) -> sigma *. standard_normal () +. mu;;

(* simulation d'une InvGamma(a., b.) *)

let sgn x = if x < 0. then -1. else 1.

let pareto_cdf = fun (x: float) (a: float) ->
    a /. x**(a +. 1.)

let pareto = fun (a: float) -> 
    let u = Random.float 1. 
    in
    1. / (1. -. u)**(1. /. a)

let inv_gamma_pdf = fun (x: float) (a: float) (b: float) ->
    ( 1. / Lanczos.gamma a ) *. b ** a *. x ** (-a -. 1.) *. exp (-b /. x)

let inv_gamma = fun a b ->
    let m = b ** a / a in let x = pareto a in let u = Random.float 1.
    in
    let rec inv_gamma_aux x u = match x u with 
        |x, u when u *. m *. pareto_pdf x a > inv_gamma_pdf x a b ->  inv_gamma_aux (pareto a) (Random.float 1.)
        |x, _                                                     -> x
    in inv_gamma_aux x u;;




(* gibbs sampling *)

(* we use another formula for the sum that depends on theta *)

let sum_over_theta = fun (theta: float array) -> 
    let capital_k = Array.length(theta) 
    in
    let capital_k_float = float_of_int(capital_k)
    in 
    let sum_of_squares = ref 0 and sum = ref 0
    in
    for i = 0 to capital_k - 1 do
        sum_of_squares := !sum_of_squares + theta.(i)**2
        sum := !sum + theta.(i)
    done;
    (!sum_of_squares -. !sum /. capital_k_float) /. 2. , !sum /. capital_k_float;;

let gibbs_sampler = fun (theta: float array) ->
    let capital_k_int = Array.length(theta) 
    in
    let capital_k = float_of_int(capital_k_int)
    in
    let (sum_gamma, mean_mu) = sum_over_theta theta
    in
    (* simulation of A *)
    let capital_a = inv_gamma (a +. (capital_k -. 1.) /. 2.) (b +. sum_for_gamma theta)
    and
    (* simulation of mu *)
    mu = normal mean_mu (capital_a /. capital_k)
    in
    (* computation of the constants *)
    let inv_v_plus_a = 1. /. (capital_v +. capital_a)
    in
    (* simulation of theta *)
    for i = 0 to capital_k_int do 
        theta.(i) <- normal (inv_v_plus_a *. (mu *. capital_v +. z.(i))) (inv_v_plus_a *. capital_a *. capital_v)
    done;;
    capital_a, mu, theta

(* coupled gibbs sampler with maximal coupling *)

let coupled_gibbs_sampler = fun (theta: float array) ->

    let maximal_coupling = fun sampler_p p sampler_q q ->
        (* sampler_p: function that sample from the law that has distribution p
           p: first distribution, val p : float -> float
           sampler_q: function that sample from the law that has distribution q
           q: second distribution, val p : float -> float
         *)
        let capital_x = sampler_p theta
        in
        let capital_w = Random.float p capital_x
        in 
        if capital_w <= q capital_x then capital_x, capital_x
        else (
            let y_star = sampler_q theta
            in
            let w_star = Random.float q y_star
            in 
            let rec aux_step_2 y_star w_star = match y_star w_star with
                |y_star, w_star when w_star <= p y_star -> begin (* la terminaison est assurée p.s. car le temps d'arret suit une loi géométrique 
                                                                    de parametre : 1− ́int_X( min (p(y), q(y)) dy )  *)
                                                            let y_star = sampler_q theta
                                                            in
                                                            let w_star = Random.float q y_star
                                                            in aux_step_2 y_star w_star
                                                            end
                |y_star,_                              -> capital_x, y_star
            in
            aux_step_2 y_star w_star;;

        );;

 (* in our case, sampler_p and sampler_q will be the gibbs_sampler, that yields the random vector (A, mu, theta_i)
    both p and q will be the distribution of the random vector (A, mu, theta_i)  *)


















