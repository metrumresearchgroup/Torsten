open Core_kernel

let pmx_ode_func = [ ( UnsizedType.AutoDiffable
                     , UnsizedType.UFun
                         ( [ (UnsizedType.AutoDiffable, UnsizedType.UReal)
                           ; (UnsizedType.AutoDiffable, UArray UReal)
                           ; (UnsizedType.AutoDiffable, UArray UReal)
                           ; (DataOnly, UArray UReal); (DataOnly, UArray UInt) ]
                         , ReturnType (UArray UReal) ) ) ]

let pmx_integrate_ode_arg = [ (UnsizedType.AutoDiffable, UnsizedType.UArray UReal) (* y0 *)
                            ; (UnsizedType.AutoDiffable, UReal)        (* t0 *)
                            ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal) (* ts *)
                            ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal) (* theta *)
                            ; (DataOnly, UnsizedType.UArray UReal)     (* x_r *)
                            ; (DataOnly, UnsizedType.UArray UInt) ]    (* x_i *)

let pmx_integrate_ode_tol = [ (UnsizedType.DataOnly, UnsizedType.UReal)            (* rtol *)
                            ; (UnsizedType.DataOnly, UnsizedType.UReal)            (* atol *)
                            ; (UnsizedType.DataOnly, UnsizedType.UReal) ]         (* mxstep *)

let pmx_algebra_sol_tol = pmx_integrate_ode_tol

let pmx_integrate_ode_group_arg = [ (UnsizedType.AutoDiffable, (UnsizedType.UArray (UnsizedType.UArray UReal))) (* y0 *)
                                  ; (UnsizedType.AutoDiffable, UReal)                   (* t0 *)
                                  ; (DataOnly, UnsizedType.UArray UInt)                 (* len *)
                                  ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* ts *)
                                  ; (UnsizedType.AutoDiffable, (UnsizedType.UArray (UnsizedType.UArray UReal))) (* theta *)
                                  ; (DataOnly, (UnsizedType.UArray (UnsizedType.UArray UReal)))     (* x_r *)
                                  ; (DataOnly, (UnsizedType.UArray (UnsizedType.UArray UInt))) ]    (* x_i *)

let rec pmx_ode add_func name args_list = match args_list with
  | [] -> ()
  | head::tail -> add_func(name, UnsizedType.ReturnType (UArray (UArray UReal)), pmx_ode_func@head); pmx_ode add_func name tail

let rec pmx_ode_group add_func name args_list = match args_list with
  | [] -> ()
  | head::tail -> add_func(name, UnsizedType.ReturnType UMatrix, pmx_ode_func@head); pmx_ode_group add_func name tail

(* from rosettacode.org *)
let rec cart_prod_list l = 
    (* We need to do the cross product of our current list and all the others
     * so we define a helper function for that *)
    let rec aux ~acc l1 l2 = match l1, l2 with
    | [], _ | _, [] -> acc
    | h1::t1, h2::t2 -> 
        let acc = (h1::h2)::acc in
        let acc = (aux ~acc t1 l2) in
        aux ~acc [h1] t2
    (* now we can do the actual computation *)
    in match l with
    | [] -> []
    | [l1] -> List.map ~f:(fun x -> [x]) l1
    | l1::tl ->
        let tail_product = cart_prod_list tl in
        aux ~acc:[] l1 tail_product

let rec pmx_solve_cpt add_func name args_list =
  let pmx_event_args = [ (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* time *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* amt *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* rate *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* ii *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt)                 (* evid *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt)                 (* cmt *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt)                 (* addl *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt) ] in            (* ss *)
  match args_list with
  | [] -> ()
  | head::tail -> add_func(name, UnsizedType.ReturnType UMatrix, pmx_event_args @ head); pmx_solve_cpt add_func name tail

let rec pmx_solve add_func name args_list =
  let pmx_event_args = [ (UnsizedType.DataOnly, UnsizedType.UInt)                        (* nCmt *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* time *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* amt *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* rate *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* ii *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt)                 (* evid *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt)                 (* cmt *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt)                 (* addl *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt) ] in            (* ss *)
  match args_list with
  | [] -> ()
  | head::tail -> add_func(name, UnsizedType.ReturnType UMatrix, pmx_ode_func @ pmx_event_args @ head); pmx_solve add_func name tail

let rec pmx_solve_group add_func name args_list =
  let pmx_event_args = [ (UnsizedType.DataOnly, UnsizedType.UInt)                        (* nCmt *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UInt)            (* length *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* time *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* amt *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* rate *)
                       ; (UnsizedType.AutoDiffable, UnsizedType.UArray UReal)            (* ii *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt)                 (* evid *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt)                 (* cmt *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt)                 (* addl *)
                       ; (UnsizedType.DataOnly, UnsizedType.UArray UInt) ] in            (* ss *)
  match args_list with
  | [] -> ()
  | head::tail -> add_func(name, UnsizedType.ReturnType UMatrix, pmx_ode_func @ pmx_event_args @ head); pmx_solve_group add_func name tail

let (pmx_solve_args, pmx_group_args, pmx_solve_cpt_args) =
  let (pmx_param_1d, pmx_param_2d) = ( [(UnsizedType.AutoDiffable, (UnsizedType.UArray UReal))], [(UnsizedType.AutoDiffable, (UnsizedType.UArray (UArray UReal)))] ) in
  let (pmx_r_data, pmx_i_data) = ( [(UnsizedType.DataOnly, (UnsizedType.UArray (UArray UReal)))], [(UnsizedType.DataOnly, (UnsizedType.UArray (UArray UInt)))] ) in
  let pmx_param = pmx_param_1d @ pmx_param_2d in
  let pmx_solve_args_0 = cart_prod_list [ pmx_param; pmx_param; pmx_param; pmx_r_data; pmx_i_data ] @
                           cart_prod_list [ pmx_param; pmx_param; pmx_param; pmx_r_data ] @
                             cart_prod_list [ pmx_param; pmx_param; pmx_param ] @
                               cart_prod_list [ pmx_param; pmx_param ] @
                                 cart_prod_list [ pmx_param ] in
  let pmx_group_args_0 = cart_prod_list [ pmx_param_2d; pmx_param_2d; pmx_param_2d; pmx_r_data; pmx_i_data ] @
                           cart_prod_list [ pmx_param_2d; pmx_param_2d; pmx_param_2d; pmx_r_data ] @
                             cart_prod_list [ pmx_param_2d; pmx_param_2d; pmx_param_2d ] @
                               cart_prod_list [ pmx_param_2d; pmx_param_2d ] @
                                 cart_prod_list [ pmx_param_2d ] in
  let pmx_solve_args_w_ode_tol = List.map ~f:(fun l -> l @ pmx_integrate_ode_tol) pmx_solve_args_0 in 
  let pmx_solve_args_w_all_tol = List.map ~f:(fun l -> l @ pmx_integrate_ode_tol @ pmx_algebra_sol_tol) pmx_solve_args_0 in 
  let pmx_group_args_w_ode_tol = List.map ~f:(fun l -> l @ pmx_integrate_ode_tol) pmx_group_args_0 in 
  let pmx_group_args_w_all_tol = List.map ~f:(fun l -> l @ pmx_integrate_ode_tol @ pmx_algebra_sol_tol) pmx_group_args_0 in 
  (pmx_solve_args_0 @ pmx_solve_args_w_ode_tol @ pmx_solve_args_w_all_tol,
   pmx_group_args_0 @ pmx_group_args_w_ode_tol @ pmx_group_args_w_all_tol,
   cart_prod_list [ pmx_param; pmx_param; pmx_param ] @ cart_prod_list [ pmx_param; pmx_param ] @ cart_prod_list [ pmx_param ] )

(* torsten functions *)
let add_torsten_qualified add_func =

  let sol_names = List.map ~f:(fun sol -> "pmx_integrate_ode_" ^ sol) ["adams"; "bdf"; "rk45"] in
  List.iter ~f:(fun sol -> pmx_ode add_func sol [pmx_integrate_ode_arg @ pmx_integrate_ode_tol ; pmx_integrate_ode_arg ]) sol_names;

  let sol_names = List.map ~f:(fun sol -> "pmx_integrate_ode_group_" ^ sol) ["adams"; "bdf"; "rk45"] in
  List.iter ~f:(fun sol -> pmx_ode_group add_func sol [pmx_integrate_ode_group_arg @ pmx_integrate_ode_tol ; pmx_integrate_ode_group_arg ]) sol_names;

  let sol_names = List.map ~f:(fun sol -> "pmx_solve_" ^ sol) ["adams"; "bdf"; "rk45"] in
  List.iter ~f:(fun sol -> pmx_solve add_func sol pmx_solve_args) sol_names;

  let sol_names = List.map ~f:(fun sol -> "pmx_solve_group_" ^ sol) ["adams"; "bdf"; "rk45"] in
  List.iter ~f:(fun sol -> pmx_solve_group add_func sol pmx_group_args) sol_names;

  let sol_names = List.map ~f:(fun sol -> "pmx_solve_" ^ sol) ["onecpt"; "twocpt"] in
  List.iter ~f:(fun sol -> pmx_solve_cpt add_func sol pmx_solve_cpt_args) sol_names;

(* pmx_solve_linode *)
  add_func
    ( "pmx_solve_linode"
    , ReturnType UMatrix
    , [ (AutoDiffable, UArray UReal)       (* time *)
      ; (AutoDiffable, UArray UReal)       (* amt *)
      ; (AutoDiffable, UArray UReal)       (* rate *)
      ; (AutoDiffable, UArray UReal)       (* ii *)
      ; (DataOnly, UArray UInt)            (* evid *)
      ; (DataOnly, UArray UInt)            (* cmt *)
      ; (DataOnly, UArray UInt)            (* addl *)
      ; (DataOnly, UArray UInt)            (* ss *)
      ; (AutoDiffable, UArray UMatrix)     (* pMatrix *)
      ; (AutoDiffable, (UArray (UArray UReal)))       (* biovar *)
      ; (AutoDiffable, (UArray (UArray UReal))) ] ) ; (* tlag *)
  add_func
    ( "pmx_solve_linode"
    , ReturnType UMatrix
    , [ (AutoDiffable, UArray UReal)       (* time *)
      ; (AutoDiffable, UArray UReal)       (* amt *)
      ; (AutoDiffable, UArray UReal)       (* rate *)
      ; (AutoDiffable, UArray UReal)       (* ii *)
      ; (DataOnly, UArray UInt)            (* evid *)
      ; (DataOnly, UArray UInt)            (* cmt *)
      ; (DataOnly, UArray UInt)            (* addl *)
      ; (DataOnly, UArray UInt)            (* ss *)
      ; (AutoDiffable, UMatrix)            (* pMatrix *)
      ; (AutoDiffable, UArray UReal)       (* biovar *)
      ; (AutoDiffable, UArray UReal) ] ) ; (* tlag *)
  add_func
    ( "pmx_solve_linode"
    , ReturnType UMatrix
    , [ (AutoDiffable, UArray UReal)       (* time *)
      ; (AutoDiffable, UArray UReal)       (* amt *)
      ; (AutoDiffable, UArray UReal)       (* rate *)
      ; (AutoDiffable, UArray UReal)       (* ii *)
      ; (DataOnly, UArray UInt)            (* evid *)
      ; (DataOnly, UArray UInt)            (* cmt *)
      ; (DataOnly, UArray UInt)            (* addl *)
      ; (DataOnly, UArray UInt)            (* ss *)
      ; (AutoDiffable, UArray UMatrix)     (* pMatrix *)
      ; (AutoDiffable, UArray UReal)       (* biovar *)
      ; (AutoDiffable, UArray UReal) ] ) ; (* tlag *)

(* pmx_solve_onecpt_bdf *)
  add_func
    ( "pmx_solve_onecpt_bdf"
    , ReturnType UMatrix
    , [ ( AutoDiffable
        , UFun
            ( [ (AutoDiffable, UReal)
              ; (AutoDiffable, UArray UReal)
              ; (AutoDiffable, UArray UReal)
              ; (AutoDiffable, UArray UReal)
              ; (DataOnly, UArray UReal); (DataOnly, UArray UInt) ]
            , ReturnType (UArray UReal) ) )
      ; (DataOnly, UInt)                         (* nOde *)
      ; (AutoDiffable, UArray UReal)             (* time *)
      ; (AutoDiffable, UArray UReal)             (* amt *)
      ; (AutoDiffable, UArray UReal)             (* rate *)
      ; (AutoDiffable, UArray UReal)             (* ii *)
      ; (DataOnly, UArray UInt)                  (* evid *)
      ; (DataOnly, UArray UInt)                  (* cmt *)
      ; (DataOnly, UArray UInt)                  (* addl *)
      ; (DataOnly, UArray UInt)                  (* ss *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* pMatrix *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* biovar *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* tlag *)
      ; (DataOnly,  UReal)                      (* rtol *)
      ; (DataOnly,  UReal)                      (* atol *)
      ; (DataOnly,  UReal)                       (* mxstep *)
      ; (DataOnly,  UReal)                      (* as_rtol *)
      ; (DataOnly,  UReal)                      (* as_atol *)
      ; (DataOnly,  UReal) ] ) ;                 (* as_mxstep *)
(* pmx_solve_onecpt_rk45 *)
  add_func
    ( "pmx_solve_onecpt_rk45"
    , ReturnType UMatrix
    , [ ( AutoDiffable
        , UFun
            ( [ (AutoDiffable, UReal)
              ; (AutoDiffable, UArray UReal)
              ; (AutoDiffable, UArray UReal)
              ; (AutoDiffable, UArray UReal)
              ; (DataOnly, UArray UReal); (DataOnly, UArray UInt) ]
            , ReturnType (UArray UReal) ) )
      ; (DataOnly, UInt)                         (* nOde *)
      ; (AutoDiffable, UArray UReal)             (* time *)
      ; (AutoDiffable, UArray UReal)             (* amt *)
      ; (AutoDiffable, UArray UReal)             (* rate *)
      ; (AutoDiffable, UArray UReal)             (* ii *)
      ; (DataOnly, UArray UInt)                  (* evid *)
      ; (DataOnly, UArray UInt)                  (* cmt *)
      ; (DataOnly, UArray UInt)                  (* addl *)
      ; (DataOnly, UArray UInt)                  (* ss *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* pMatrix *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* biovar *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* tlag *)
      ; (DataOnly,  UReal)                      (* rtol *)
      ; (DataOnly,  UReal)                      (* atol *)
      ; (DataOnly,  UReal)                       (* mxstep *)
      ; (DataOnly,  UReal)                      (* as_rtol *)
      ; (DataOnly,  UReal)                      (* as_atol *)
      ; (DataOnly,  UReal) ] ) ;                 (* as_mxstep *)
(* pmx_solve_twocpt_bdf *)
  add_func
    ( "pmx_solve_twocpt_bdf"
    , ReturnType UMatrix
    , [ ( AutoDiffable
        , UFun
            ( [ (AutoDiffable, UReal)
              ; (AutoDiffable, UArray UReal)
              ; (AutoDiffable, UArray UReal)
              ; (AutoDiffable, UArray UReal)
              ; (DataOnly, UArray UReal); (DataOnly, UArray UInt) ]
            , ReturnType (UArray UReal) ) )
      ; (DataOnly, UInt)                         (* nOde *)
      ; (AutoDiffable, UArray UReal)             (* time *)
      ; (AutoDiffable, UArray UReal)             (* amt *)
      ; (AutoDiffable, UArray UReal)             (* rate *)
      ; (AutoDiffable, UArray UReal)             (* ii *)
      ; (DataOnly, UArray UInt)                  (* evid *)
      ; (DataOnly, UArray UInt)                  (* cmt *)
      ; (DataOnly, UArray UInt)                  (* addl *)
      ; (DataOnly, UArray UInt)                  (* ss *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* pMatrix *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* biovar *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* tlag *)
      ; (DataOnly,  UReal)                      (* rtol *)
      ; (DataOnly,  UReal)                      (* atol *)
      ; (DataOnly,  UReal)                       (* mxstep *)
      ; (DataOnly,  UReal)                      (* as_rtol *)
      ; (DataOnly,  UReal)                      (* as_atol *)
      ; (DataOnly,  UReal) ] ) ;                 (* as_mxstep *)
(* pmx_solve_twocpt_rk45 *)
  add_func
    ( "pmx_solve_twocpt_rk45"
    , ReturnType UMatrix
    , [ ( AutoDiffable
        , UFun
            ( [ (AutoDiffable, UReal)
              ; (AutoDiffable, UArray UReal)
              ; (AutoDiffable, UArray UReal)
              ; (AutoDiffable, UArray UReal)
              ; (DataOnly, UArray UReal); (DataOnly, UArray UInt) ]
            , ReturnType (UArray UReal) ) )
      ; (DataOnly, UInt)                         (* nOde *)
      ; (AutoDiffable, UArray UReal)             (* time *)
      ; (AutoDiffable, UArray UReal)             (* amt *)
      ; (AutoDiffable, UArray UReal)             (* rate *)
      ; (AutoDiffable, UArray UReal)             (* ii *)
      ; (DataOnly, UArray UInt)                  (* evid *)
      ; (DataOnly, UArray UInt)                  (* cmt *)
      ; (DataOnly, UArray UInt)                  (* addl *)
      ; (DataOnly, UArray UInt)                  (* ss *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* pMatrix *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* biovar *)
      ; (AutoDiffable, (UArray (UArray UReal))) (* tlag *)
      ; (DataOnly,  UReal)                      (* rtol *)
      ; (DataOnly,  UReal)                      (* atol *)
      ; (DataOnly,  UReal)                       (* mxstep *)
      ; (DataOnly,  UReal)                      (* as_rtol *)
      ; (DataOnly,  UReal)                      (* as_atol *)
      ; (DataOnly,  UReal) ] ) ;                 (* as_mxstep *)
(* linear interpolation *)
  add_func
    ( "linear_interpolation"
    , ReturnType UReal
    , [ (AutoDiffable, UReal)              (* x_out *)
      ; (AutoDiffable, UArray UReal)       (* x *)
      ; (AutoDiffable, UArray UReal) ] ) ; (* y *)
  add_func
    ( "linear_interpolation"
    , ReturnType (UArray UReal)
    , [ (AutoDiffable, UArray UReal)       (* x_out *)
      ; (AutoDiffable, UArray UReal)       (* x *)
      ; (AutoDiffable, UArray UReal) ] )   (* y *)

let () =
  add_torsten_qualified Stan_math_signatures.add_qualified

(* end of torsten sig helpers *)

       
