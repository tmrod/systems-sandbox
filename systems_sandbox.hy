(import sympy)
(import sympy.plotting.plot [LineOver1DRangeSeries])
(import sympy.physics.control :as control)
(import sympy.physics.control [control-plots])
(import sympy.physics.control.lti [TransferFunction Feedback Series Parallel])
(import sympy.abc *)

(import networkx :as nx)

(import functools [reduce])
(import operator [mul add])

(import numpy :as np)
(import scipy.signal [residuez])
(import scipy.special [binom])
(import matplotlib.pyplot :as plt)
(import matplotlib.pyplot [axvline axhline])
(import matplotlib [patches])

(plt.style.use "Solarize_Light2")

;; some special symbols
(setv j (sympy.sqrt -1))
(setv i (sympy.sqrt -1))

;; declaration of systems
(defmacro declare [S #* H]
  `(setx ~S (~@H)))

(defn rational [[gain 1] [zeros []] [poles []]] `(system (rat ~gain ~zeros ~poles)))
(defn rational-polynomial [[numerator "1"] [denominator "1"]]
  `(system (rat-poly ~(eval numerator) ~(eval denominator))))

;; compositional forms
(defn compose [#* l] (+ ['compose] (list l)))
(defn feedback [S T] ['feedback S T])
(defn parallel [#* l] (+ ['sum] (list l)))
(defn sum [#* l] (+ ['sum] (list l)))

;; common CT systems
(defn gain [[gain 1]] `(system (rat ~gain [] [])))
(defn derivative [] `(system (rat 1 [0] [])))
(defn integrate [] `(system (rat 1 [] [0])))

;; common DT systems
(defn delay [[m 1]] (rational-polynomial "1" f"z**{m}"))
(defn accumulate [] (rational-polynomial "z" "z-1"))

;; simple controllers
(defn pd [[Kp 1] [Kd 0] [sym s]]
  (let [dd (cond (= sym z) (delay) True (derivative))]
    (parallel (gain Kp)
	      (compose (gain Kd)
		       dd))))

(defn pid [[Kp 1] [Ki 0] [Kd 0] [sym s]]
  (let [ii (cond (= sym z) (accumulate) True (integrate))]
    (parallel (pd Kp Kd)
	      (compose (gain Ki)
		       ii))))

;; compute symbolic transfer function of system
(defn cancel [tf [sym s]] (TransferFunction.from_rational_expression (.
								      (. tf (to-expr))
								      (cancel)) sym))

(defn -transfer-function [S bind sym]
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
      (. (TransferFunction (reduce mul (map (fn [zz] (- sym zz)) zs) k)
			   (reduce mul (map (fn [pp] (- sym pp)) ps) 1)
			   sym)
	 (subs bind)))

    (= (get (get S 1) 0) 'rat-poly)
    (let [num (get (get S 1) 1)
	      den (get (get S 1) 2)]
      (. (TransferFunction.from-rational-expression (/ num den) sym)
	 (subs bind)))

    True
    (raise ValueError))

   ;; feedback system
   (= (get S 0) 'feedback)
   (. (Feedback (transfer-function (get S 1) bind sym)
		(transfer-function (get S 2) bind sym))
      (doit))

   ;; composed systems
   (= (get S 0) 'compose)
   (. (Series (unpack-iterable (map (fn [T] (transfer-function T bind sym))
				    (cut S 1 None))))
      (doit))

   ;; both parallel and generic summation
   (= (get S 0) 'sum)
   (. (Parallel (unpack-iterable (map (fn [T] (transfer-function T bind sym))
				      (cut S 1 None))))
      (doit))))

(defn transfer-function [S [bind {}] [sym s]]
  (cancel (-transfer-function S bind sym) sym))

;; system properties
(defn zeros [S [bind {}] [sym s]]
  (. (transfer-function S bind sym) (zeros)))
(defn poles [S [bind {}] [sym s]]
  (. (transfer-function S bind sym) (poles)))

;; TODO: refactor cond statements to be nicer
(defn roc-lower-bound [S [bind {}] [sym s]]
  (cond
   (= sym z)
   (let [ps (poles S bind sym)]
     (cond
      (= ps [])
      (- sympy.oo)

      True
      (sympy.Max (unpack-iterable (map sympy.Abs ps)))))

   True
   (let [ps (poles S bind sym)]
     (cond
      (= ps [])
      (- sympy.oo)

      True
      (sympy.Max (unpack-iterable (map sympy.re ps)))))))

(defn region-of-convergence [S [bind {}] [sym s]]
  (let [lb (roc-lower-bound S bind sym)]
    (cond
     (= sym z)
     (sympy.StrictGreaterThan (sympy.Abs sym) lb)

     True
     (sympy.StrictGreaterThan (sympy.re sym) lb))))

;; q & a
(defmacro ?is [S #* p]
  (cond
   (= p [])
   (print "Ill-posed query.")

   (in (get p 0) `(stable convergent-at))
   `(~@p ~S)

   True
   (print "Ill-posed query.")))

(defn convergent-at [pt S [bind {}] [sym s]]
  (let [lb (roc-lower-bound S bind sym)]
    (cond
     (is lb None)
     True

     (= sym z)
     (sympy.StrictGreaterThan (sympy.Abs pt) lb)

     True
     (sympy.StrictGreaterThan (sympy.re pt) lb))))

;; note: stable "reimplements" roc-lower-bound
;; this is done to handle things like
;;   0>Max(-1,re(a))
;; being the same as
;;   0>re(a)
(defn stable [S [bind {}] [sym s]]
  (let [ps (poles S bind sym)]
    (cond
     (= ps [])
     True

     (= sym z)
     (sympy.StrictLessThan (sympy.Max (unpack-iterable (map
							sympy.Abs
							(filter
							 (fn [p] (!= (sympy.LessThan 1
										     (sympy.Abs p))
								     False))
							 ps))))
			   1)

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
(defn pole-zero-plot [S [bind {}] [sym s]]
  (let [pzdata (control-plots.pole-zero-numerical-data (transfer-function S bind sym))
	       zdata (get pzdata 0)
	       pdata (get pzdata 1)
	       maxz (if (> (len zdata) 0) (np.max (np.abs zdata)) 1)
	       maxp (if (> (len pdata) 0) (np.max (np.abs pdata)) 1)
	       modrange (* 1.1 (np.maximum 1 (np.maximum maxz maxp)))
	       fig (plt.figure)
	       ax (. fig (add-subplot 1 1 1 :aspect "equal"))
	       unitcirc (patches.Circle #(0 0) 1
					:fill False :color "black" :ls "dashed" :alpha 0.5)]
    (axvline 0 :color "0.7")
    (axhline 0 :color "0.7")
    (. ax (set-xlim (* -1 modrange) modrange))
    (. ax (set-xlabel "Re"))
    (. ax (set-ylabel "Im"))
    (. ax (set-ylim (* -1 modrange) modrange))
    (. ax (set-title f"Poles and zeros of ${(sympy.printing.latex (transfer-function S bind sym))}$" :pad 20))
    (when (= sym z)
      (. ax (add-patch unitcirc)))
    (plt.plot (np.real pdata) (np.imag pdata)
	      "x" :markersize 9 :alpha 0.8)
    (plt.plot (np.real zdata) (np.imag zdata)
	      "o" :markersize 9 :alpha 0.8)
    (plt.show)))

(defn -plot [xs hs [title ""] [xlab ""] [ylab ""]]
  (let [fig (plt.figure)
	    ax (. fig (add-subplot 1 1 1))]
    (axvline 0 :color "0.7")
    (axhline 0 :color "0.7")
    (. ax (set-xlabel xlab))
    (. ax (set-ylabel ylab))
    (. ax (set-title title :pad 20))
    (plt.plot xs hs "-")
    (plt.show)))

(defn impulse-plot [S [bind {}] [sym s] #** kwargs]
  (cond
   (= sym z)
   (print "DT impulse response not yet implemented.")

   True
   (let [[xs hs] (control-plots.impulse-response-numerical-data
		  (transfer-function S bind sym))]
     (-plot xs hs f"Impulse response of ${(sympy.printing.latex (transfer-function S bind sym))}$" "t"))))

(defn step-plot [S [bind {}] [sym s] #** kwargs]
  (cond
   (= sym z)
   (print "DT step response not yet implemented.")

   True
   (let [[xs hs] (control-plots.step-response-numerical-data
		  (transfer-function S bind sym))]
     (-plot xs hs f"Step response of ${(sympy.printing.latex (transfer-function S bind sym))}$" "t"))))

(defn ramp-plot [S [bind {}] [sym s] #** kwargs]
  (cond
   (= sym z)
   (print "DT ramp response not yet implemented.")

   True
   (let [[xs hs] (control-plots.ramp-response-numerical-data
		  (transfer-function S bind sym))]
     (-plot xs hs f"Ramp response of ${(sympy.printing.latex (transfer-function S bind sym))}$" "t"))))

(defn frequency-plot [S [bind {}] [sigma0 0] [sym s] #** kwargs]
  (cond
   (= sym z)
   (if (convergent-at (sympy.exp sigma0) S bind sym)
       (let [[xs hs] (. (LineOver1DRangeSeries
			 (sympy.Abs (. (transfer-function S bind sym)
				       (to-expr)
				       (subs {sym (sympy.exp (+ sigma0 (* 1j omega)))})))
			 #(omega (- sympy.pi) sympy.pi)) (get-points))]
	 (-plot xs hs f"Frequency response of ${(sympy.printing.latex (transfer-function S bind sym))}$" f"$\\omega$"))
     (print "System not stable for given parameters."))

   True
   (if (convergent-at sigma0 S bind)
       (let [[xs hs] (. (LineOver1DRangeSeries
			 (sympy.Abs (. (transfer-function S bind sym)
				       (to-expr)
				       (subs {sym (+ sigma0 (* 1j omega))})))
			 #(omega -10 10)) (get-points))]
     	 (-plot xs hs f"Frequency response of ${(sympy.printing.latex (transfer-function S bind sym))}$" f"$\\omega$"))
     (print "System not stable for given parameters."))))
