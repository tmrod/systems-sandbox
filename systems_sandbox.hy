(import sympy)
(import sympy.physics.control :as control)
(import sympy.physics.control [control-plots])
(import sympy.physics.control.lti [TransferFunction Feedback Series Parallel])
(import sympy.abc *)
(import networkx :as nx)
(import functools [reduce])
(import operator [mul add])

;; some special symbols
(setv j (sympy.sqrt -1))
(setv i (sympy.sqrt -1))

;; declaration of systems
(defmacro declare [S #* H]
  `(setx ~S (~@H)))

(defn rational [[gain 1] [zeros []] [poles []]] `(system (rat ~gain ~zeros ~poles)))
(defn rational-polynomial [[numerator "1"] [denominator "1"]]
  `(system (rat-poly ~(eval numerator) ~(eval denominator))))

(defn gain [[gain 1]] `(system (rat ~gain [] [])))
(defn derivative [] `(system (rat 1 [0] [])))
(defn integrate [] `(system (rat 1 [] [0])))

;; compositional forms
(defn compose [#* l] (+ ['compose] (list l)))
(defn feedback [S T] ['feedback S T])
(defn parallel [#* l] (+ ['sum] (list l)))
(defn sum [#* l] (+ ['sum] (list l)))

;; compute symbolic transfer function of system
(defn cancel [tf] (TransferFunction.from_rational_expression (.
							      (. tf (to-expr))
							      (cancel)) s))

(defn transfer-function [S [bind {}]] (cancel (-transfer-function S bind)))

(defn -transfer-function [S [bind {}]]
  (cond
   ;; basic system
   (= (get S 0) 'system)
   (cond
    ;; rational system
    (= (get (get S 1) 0) 'rat)
    (let [G (get S 1)
	    k (get G 1)
	    zs (get G 2)
	    ps (get G 3)]
      (. (TransferFunction (reduce mul (map (fn [z] (- s z)) zs) k)
			   (reduce mul (map (fn [p] (- s p)) ps) 1)
			   s)
	 (subs bind)))

    (= (get (get S 1) 0) 'rat-poly)
    (let [num (get (get S 1) 1)
	      den (get (get S 1) 2)]
      (. (TransferFunction.from-rational-expression (/ num den) s)
	 (subs bind)))

    True
    (raise ValueError))

   ;; feedback system
   (= (get S 0) 'feedback)
   (. (Feedback (transfer-function (get S 1) bind)
		(transfer-function (get S 2) bind))
      (doit))

   ;; composed systems
   (= (get S 0) 'compose)
   (. (Series (unpack-iterable (map (fn [T] (transfer-function T bind))
				    (cut S 1 None))))
      (doit))

   ;; both parallel and generic summation
   (= (get S 0) 'sum)
   (. (Parallel (unpack-iterable (map (fn [T] (transfer-function T bind))
				      (cut S 1 None))))
      (doit))))

;; system properties
(defn zeros [S [bind {}]]
  (. (transfer-function S bind) (zeros)))
(defn poles [S [bind {}]]
  (. (transfer-function S bind) (poles)))

(defn roc-lower-bound [S [bind {}]]
  (let [ps (poles S bind)]
    (cond
     (= ps [])
     (- sympy.oo)

     True
     (sympy.Max (unpack-iterable (map sympy.re ps))))))

(defn region-of-convergence [S [bind {}]]
  (let [lb (roc-lower-bound S bind)]
     (sympy.StrictGreaterThan (sympy.re s) lb)))

;; q & a
(defmacro ?is [S #* p]
  (cond
   (= p [])
   "Ill-posed query."

   (in (get p 0) `(stable convergent-at))
   `(~@p ~S)

   True
   "Ill-posed query."))

(defn convergent-at [s S [bind {}]]
  (let [lb (roc-lower-bound S bind)]
    (cond
     (is lb None)
     True

     True
     (sympy.StrictGreaterThan (sympy.re s) lb))))

;; note: stable "reimplements" roc-lower-bound
;; this is done to handle things like
;;   0>Max(-1,re(a))
;; being the same as
;;   0>re(a)
(defn stable [S [bind {}]]
  (let [ps (poles S bind)]
    (cond
     (= ps [])
     True

     True
     (sympy.StrictLessThan (sympy.Max (unpack-iterable (map
							sympy.re
							(filter
							 (fn [p] (!= (sympy.LessThan 0
										     (sympy.re p))
								     False))
							 ps))))
			   0))))

;; plotting functions
(defn pole-zero-plot [S [bind {}] #** kwargs]
  (control-plots.pole-zero-plot
   (transfer-function S bind)
   :kwargs kwargs))

(defn impulse-plot [S [bind {}] #** kwargs]
  (control-plots.impulse-response-plot
   (transfer-function S bind)
   :kwargs kwargs))

(defn step-plot [S [bind {}] #** kwargs]
  (control-plots.step-response-plot
   (transfer-function S bind)
   :kwargs kwargs))

(defn ramp-plot [S [bind {}] #** kwargs]
  (control-plots.ramp-response-plot
   (transfer-function S bind)
   :kwargs kwargs))

(defn frequency-plot [S [bind {}] [sigma0 0] #** kwargs]
  (if (convergent-at sigma0 S bind)
      (sympy.plot (abs
		   (. (transfer-function S bind)
		      (to-expr)
		      (subs {s (+ sigma0 (* 1j omega))})))
		  :kwargs kwargs)
    "System not stable for given parameters."))
