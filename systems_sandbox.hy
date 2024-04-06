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
(import matplotlib.widgets [Slider])

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

;; num-funs is a list of functions
;; plot-funs is a corresponding list of plotting functions
(defn -plot-with-sliders [num-funs tf [sym s]
				   [lbnd 0] [ubnd None]
				   [xlim-init [0 10]] [ylim-init [0 1]]
				   * [plot-funs [(fn [ax x y] (ax.plot x y))]]
				   [ax-setup-fun (fn [ax] None)]
				   [annotate-stability True]]
  (assert (= (len num-funs) (len plot-funs)) "num-funs and plot-funs should be the same length")
  (let [fig (plt.figure)
	    ax (fig.add-subplot 1 1 1)
	    freesyms (tf.free-symbols.difference [sym])
	    slider-labels (list (map str freesyms))]
    ;; set up figure and slider hardware
    (ax-setup-fun ax)
    (when (> (len freesyms) 0) (fig.subplots-adjust :bottom 0.25))
    (setv slider-axs (map (fn [idx-label] (fig.add-axes [0.15
							 (* 0.05 (+ 1 (get idx-label 0)))
							 0.7 0.05]))
			  (enumerate slider-labels)))
    (setv sliders (list (map (fn [x] (Slider (unpack-iterable x) -5 5 :valinit 0))
			     (zip slider-axs slider-labels))))

    (setv lines (lfor plot-fun plot-funs (plot-fun ax
						   (np.linspace (unpack-iterable xlim-init) 3)
						   (np.linspace (unpack-iterable ylim-init) 3))))

    (setv stable-annotate (ax.annotate "" #(0.5 1) :xycoords "axes fraction"
				       :fontsize "x-large"
				       :horizontalalignment "center"
				       :verticalalignment "top"))

    (defn update-fun [[v None]]
      (let [[xmin xmax] (ax.get-xlim)
	    var-start-end #(sym (np.clip xmin lbnd ubnd) (np.clip xmax lbnd ubnd))
	    sub-dict (dfor x (zip freesyms sliders) (get x 0) (. (get x 1) val))
	    xyss (lfor num-fun num-funs (num-fun (tf.subs sub-dict) var-start-end))]

	(for [[line xys] (zip lines xyss)]
	     (. (get line 0) (set-xdata (get xys 0)))
	     (. (get line 0) (set-ydata (get xys 1))))

	(when annotate-stability
	  (if (and
	       (!= sym z)
	       (= (. tf (subs sub-dict) (is-stable)) True)
	       (= (. tf (subs sub-dict) is-strictly-proper) True))
	      (do
		  (stable-annotate.set-text "")
		  (stable-annotate.set-color "green"))
	    (do
		(stable-annotate.set-text "possibly unstable")
		(stable-annotate.set-color "red"))))

	(ax.relim)
	(ax.autoscale-view)

	(fig.canvas.draw-idle)))

    (update-fun)

    (for [sl sliders] (sl.on-changed update-fun))

    (plt.show)))

(defn response-plot [S [bind {}] [sym s]
		       [numdatafun control-plots.impulse-response-numerical-data]
		       [responsename "Impulse"]
		       #** kwargs]
  (defn ax-setup-fun [ax]
    (axvline 0 :color "0.7")
    (axhline 0 :color "0.7")
    (ax.set-xlabel "$t$")
    (ax.set-title f"{responsename} response of ${(sympy.printing.latex (transfer-function S bind sym))}$" :pad 20))

  (cond
   (= sym z)
   (print f"DT {responsename} response not yet implemented.")

   True
   (-plot-with-sliders
    [(fn [tf var-start-end] (numdatafun tf
					:lower-limit (get var-start-end 1)
					:upper-limit (get var-start-end 2)))]
    (transfer-function S :bind bind :sym sym) s
    :lbnd 0 :ubnd 50
    :ax-setup-fun ax-setup-fun)))

(defn impulse-plot [S [bind {}] [sym s] #** kwargs]
  (response-plot S bind sym
		 control-plots.impulse-response-numerical-data
		 "Impulse"
		 :kwargs kwargs))

(defn step-plot [S [bind {}] [sym s] #** kwargs]
  (response-plot S bind sym
		 control-plots.step-response-numerical-data
		 "Step"
		 :kwargs kwargs))

(defn ramp-plot [S [bind {}] [sym s] #** kwargs]
  (response-plot S bind sym
		 control-plots.ramp-response-numerical-data
		 "Ramp"
		 :kwargs kwargs))

(defn frequency-plot [S [bind {}] [sigma0 0] [sym s] #** kwargs]
  (defn ax-setup-fun [ax]
    (axvline 0 :color "0.7")
    (axhline 0 :color "0.7")
    (ax.set-xlabel "$\\omega$")
    (ax.set-ylabel "$|F|$")
    (ax.set-title f"Frequency response of ${(sympy.printing.latex (transfer-function S bind sym))}$" :pad 20))

  (let [tf (transfer-function S :bind bind :sym sym)]
    (cond
     (= sym z)
     (-plot-with-sliders
      [(fn [tf0 var-start-end] (. (LineOver1DRangeSeries
				   (sympy.Abs (. tf0 (to-expr) (subs {sym (sympy.exp (+ sigma0 (* 1j omega)))})))
				   #(omega (get var-start-end 1) (get var-start-end 2)))
				  (get-points)))]
      tf sym
      :lbnd (* -6 np.pi) :ubnd (* 6 np.pi)
      :xlim-init [(- np.pi) (np.pi)] :ylim-init [0 1]
      :ax-setup-fun ax-setup-fun)

     True
     (-plot-with-sliders
      [(fn [tf var-start-end] (. (LineOver1DRangeSeries
				  (sympy.Abs (. tf (to-expr) (subs {sym (+ sigma0 (* 1j omega))})))
				  #(omega (get var-start-end 1) (get var-start-end 2)))
				 (get-points)))]
      tf sym
      :lbnd -50 :ubnd 50
      :xlim-init [-10 10] :ylim-init [0 1]
      :ax-setup-fun ax-setup-fun))))

(defn pole-zero-plot [S [bind {}] [sym s]]
  (defn poles-num-fun [tf var-start-end]
    (let [[zdata pdata] (control-plots.pole-zero-numerical-data tf)] (return [(np.real pdata) (np.imag pdata)])))
  (defn zeros-num-fun [tf var-start-end]
    (let [[zdata pdata] (control-plots.pole-zero-numerical-data tf)]
      (return [(np.real zdata) (np.imag zdata)])))
  (defn poles-plot-fun [ax xs ys]
    (ax.plot xs ys "x" :markersize 9 :alpha 0.8))
  (defn zeros-plot-fun [ax xs ys]
    (ax.plot xs ys "o" :markersize 9 :alpha 0.8 :fillstyle "none"))

  (defn ax-setup-fun [ax]
    (axvline 0 :color "0.7")
    (axhline 0 :color "0.7")
    (ax.set-aspect "equal")
    (ax.set-xlabel "Re")
    (ax.set-ylabel "Im")
    (ax.set-title f"Poles and zeros of ${(sympy.printing.latex (transfer-function S bind sym))}$" :pad 20)
    (if (= sym z)
	(ax.add-patch (patches.Circle #(0 0) 1 :fill False :color "black" :ls "dashed" :alpha 0.5))
      (ax.add-patch (patches.Circle #(0 0) 0.2 :fill False :alpha 0.0))))

  (setv num-funs [poles-num-fun zeros-num-fun])
  (setv plot-funs [poles-plot-fun zeros-plot-fun])

  (-plot-with-sliders num-funs (transfer-function S bind sym) sym
		      :plot-funs plot-funs
		      :ax-setup-fun ax-setup-fun
		      :lbnd (- np.inf) :ubnd None
		      :xlim-init [-1.2 1.2] :ylim-init [-1.2 1.2]
		      :annotate-stability False))
