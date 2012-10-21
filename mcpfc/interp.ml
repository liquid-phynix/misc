open Batteries_uni
open Bigarray
open Unix
open Gsl_integration
open Printf

let pi = 2.0 *. asin 1.0

let print_array2d ar =
  for i = 0 to (Array2.dim1 ar - 1) do
    for j = 0 to (Array2.dim2 ar - 1) do
      print_float (Array2.get ar i j);
      print_string " "
    done;
    print_newline () 
  done

let read_array2d fn dim =
  let fd = openfile fn [O_RDONLY] 0 in
  let ar2d = Array2.map_file fd float64 c_layout false dim dim in
  close fd;
  ar2d

let save_array2d fn ar =
  let fd = openfile fn [O_RDWR; O_CREAT; O_TRUNC] 0o640
  and d1 = Array2.dim1 ar
  and d2 = Array2.dim2 ar in
  let ar2d = Array2.map_file fd float64 c_layout true d1 d2 in
  for i = 0 to (d1 - 1) do
    for j = 0 to (d2 - 1) do
      Array2.set ar2d i j (Array2.get ar i j)
    done
  done;
  close fd

let make_interp_func_unif ar =
  let dim = float_of_int (Array2.dim1 ar) in
  fun x y ->
    if x >= 0.0 && x < dim && y >= 0.0 && y <= dim then
      let x = int_of_float x
      and y = int_of_float y in
      Array2.get ar x y
    else 0.0

let make_interp_func_lin ar =
  let dim = float_of_int (Array2.dim1 ar) in
  fun x y ->
    if x >= 0.0 && x < dim && y >= 0.0 && y < dim then
      let xi = int_of_float x
      and yi = int_of_float y in
      let dx = x -. float_of_int xi
      and dy = y -. float_of_int yi in
      let xip1 = xi + 1
      and yip1 = yi + 1 in
      let pA = Array2.get ar xi yi
      and pB = Array2.get ar xip1 yi
      and pC = Array2.get ar xi yip1
      and pD = Array2.get ar xip1 yip1 in
      pA *. (1.0 -. dx) *. (1.0 -. dy) +. pB *. dx *. (1.0 -. dy) +. pC *. (1.0 -. dx) *. dy +. pD *. dx *. dy
    else 0.0
    
let integrate_fun f lower upper =
  let ws = make_ws 10000 in
  let res = qag f lower upper 1e-2 1e-2 GAUSS15 ws
  in res.Gsl_fun.res

let save_points fn points =
  let oc = open_out fn in
  List.iter (fun (x, y) -> fprintf oc "%f\t%f\n" x y) points;
  close_out oc

(* real purpose main *)

let () = match Sys.argv with
  | [|_; kind; in_file; out_file; len; lower; upper; dk|] ->
      let lower = float_of_string lower
      and upper = float_of_string upper
      and dk = float_of_string dk
      and len = int_of_string len in
      let ar = read_array2d in_file len in
      let ipfun = match kind with
        | "unif" -> make_interp_func_unif ar
        | "lin" -> make_interp_func_lin ar
        | _ -> invalid_arg "interpolation must be lin/unif" in
      let at_k k = integrate_fun (fun phi ->  ipfun (k *. cos phi) (k *. sin phi)) 0. (2. *. pi) in
      save_points out_file ([? List : (k, at_k k) | k <- ((lower,dk) --. upper) ?])
  | _ -> invalid_arg "kind(lin/unif), infile, outfile, len, lower, upper, dk"


(* testing the  interpolation *)

(* let () = match Sys.argv with *)
(*   | [|_; in_file; out_file; len; dh|] -> *)
(*       let dh = float_of_string dh in *)
(*       let len = int_of_string len in *)
(*       let new_len = int_of_float (float_of_int (len - 1) /. dh) in *)
(*       let ar = read_array2d in_file len in *)
(*       let out_ar = Array2.create float64 c_layout new_len new_len in *)
(*       let ipfun = make_interp_func_lin ar in *)
(*       for i = 0 to (new_len - 1) do *)
(*         for j = 0 to (new_len - 1) do *)
(*           Array2.set out_ar i j (ipfun (float_of_int i *. dh) (float_of_int j *. dh)) *)
(*         done *)
(*       done; *)
(*       printf "new size: %d x %d\n" new_len new_len; *)
(*       save_array2d out_file out_ar *)
(*   | _ -> invalid_arg "infile, outfile, len, dh" *)









(* let ar = read_array2d "abs_sk_499" 499;; *)
(* let at_k k = integrate_fun (fun phi -> (make_interp_func ar) (k *. cos phi) (k *. sin phi)) 0. (2. *. pi);; *)
(* save_points "sk_499" [? List : let x = (float_of_int n) /. 2.0 in (x, at_k x) | n <- 1 -- 800 ?];; *)

(* [? List : let x = (float_of_int n) /. 2.0 in (x, at_k x) | n <- 1 -- 800 ?];; *)
(* let ar = read_array2d "testmat_3x3_double" 3;; *)
(* let at_k k = integrate_fun (fun phi -> (make_interp_func ar) (k *. cos phi) (k *. sin phi)) 0. (2. *. pi);; *)
(* List.map at_k [? List : (float_of_int n) | n <- 1 -- 10 ?];; *)
(* [? List : let x = (float_of_int n) /. 10.0 in (x, at_k x) | n <- 1 -- 20 ?];; *)
(* save_points "sk_499" [? List : let x = (float_of_int n) /. 2.0 in (x, at_k x) | n <- 1 -- 800 ?];; *)


(* ocamlfind ocamlopt -package gsl -linkpkg rand.ml -inline 1000 -o rand *)
