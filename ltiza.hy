(require systems-sandbox *)
(import systems-sandbox *)

(import cmd)

(defclass LTIza [cmd.Cmd]

  (setv pfx-prompt "LTIza ")
  (setv sfx-prompt " => ")
  (setv prompt (+ pfx-prompt sfx-prompt))
  (setv symbols-list (list (map str [a b c d e f g h     k l m n o p q r   t u v w x y z
				     A B C D E F G H     K L M N O P Q R S T U V W X Y Z])))

  ;; bindings
  ;; system definitions
  ;; keyword list
  (defn __init__ [self S [sysname "S"] [bind {}]]
    (. (super) (__init__))
    (setv self.sysname sysname)
    (setv self.S S)
    (setv self.bind bind)
    (setv self.prompt (+ self.pfx-prompt (repr self.bind) self.sfx-prompt)))

  ;; show bound variables in prompt
  (defn postcmd [self stop line]
    (if stop
	(return True)
      (setv self.prompt (+ self.pfx-prompt (repr self.bind) self.sfx-prompt))))

  (defn emptyline [self]
    (return False))

  ;; catch exceptions without leaving repl
  (defn onecmd [self line]
    (try
     (return (. (super) (onecmd line)))
     (except [e []] (print (repr e)))))

  ;; binding values
  (defn helper-bind [self wl]
    (cond
     (< (len wl) 2)
     (return)

     (in (get wl 0) self.symbols-list)
     (do
	 (setv (. self.bind [(hy.eval (hy.read (get wl 0)))])
	       (hy.eval (hy.read (get wl 1))))
	 (self.helper-bind (cut wl 2 None)))

     True
     (self.helper-bind (cut wl 1 None))))

  (defn do-bind [self line]
    "Bind symbols to (possibly symbolic) values.

    Examples
    --------
    LTIza {} => bind a to 1 and b equal to 2 and c = 3
    LTIza {a:1, b:2, c:3} => bind c 4
    LTIza {a:1, b:2, c:4} => ..."
    (let [wl (.split line)
	  fwl (list (filter (fn [s] (not-in s ["=" "to" "equal"])) wl))]
      (self.helper-bind fwl)
      (return False)))

  ;; unbinding values
  (defn helper-unbind [self wl]
    (let [fwl (list (filter (fn [w] (in w self.symbols-list)) wl))]
      (cond
       (= (len fwl) 0)
       (do
	   (setv self.bind {})
	   (return))

       True
       (do
	   (list (map (fn [w] (. self.bind (pop (hy.eval (hy.read w)) None))) fwl))
	   (return)))))

  (defn do-unbind [self line]
    "Unbind provided symbols. If no symbols are given, unbind *all* symbols.

    Examples
    --------
    LTIza {a:1, b:2, c:3} => unbind c
    LTIza {a:1, b:2} => unbind everything
    LTIza {} => bind a 2 b 1
    LTIza {a:2, b:1} => unbind
    LTIza {} => ..."
    (let [wl (.split line)]
      (self.helper-unbind wl)
      (return False)))

  ;; show transfer function
  (defn do-transfer [self line]
    "Show transfer function of system.

    Examples
    --------
    LTIza {} => transfer function
    LTIza {} => what is the transfer function
    LTIza {a:1, b:2} => transfer"
    (print f"(transfer-function {self.sysname} {(hy.repr self.bind)})")
    (print (. (transfer-function self.S self.bind) (to-expr))))

  ;; show poles
  (defn do-poles [self line]
    "Show poles of system.

    Examples
    --------
    LTIza {a:1, b:2, c:3} => poles
    LTIza {a:1, b:2} => pole"
    (print f"(poles {self.sysname} {(hy.repr self.bind)})")
    (print (poles self.S self.bind)))

  ;; show zeros
  (defn do-zeros [self line]
    "Show zeros of system.

    Examples
    --------
    LTIza {a:1, b:2, c:3} => zeros
    LTIza {a:1, b:2} => zero"
    (print f"(zeros {self.sysname} {(hy.repr self.bind)})")
    (print (zeros self.S self.bind)))

  ;; plotting
  (defn do-plot [self line]
    "Show one of many types of plots. Possible choices are:
    - Pole-zero plot
    - Impulse response
    - Step response
    - Ramp response

    Examples
    --------
    LTIza {} => plot the impulse response
    LTIza {a:1, b:2} => plot the poles and zeros
    LTIza {} => plot step
    LTIza {} => ramp response"
    (cond
     (or (in "pole" line) (in "zero" line))
     (do
	 (print f"(pole-zero-plot {self.sysname} {(hy.repr self.bind)})")
	 (pole-zero-plot self.S self.bind))

     (in "step" line)
     (do
	 (print f"(step-plot {self.sysname} {(hy.repr self.bind)})")
	 (step-plot self.S self.bind))

     (in "ramp" line)
     (do
	 (print f"(ramp-plot {self.sysname} {(hy.repr self.bind)})")
	 (ramp-plot self.S self.bind))

     (or (in "impulse" line) (in "response" line))
     (do
	 (print f"(impulse-plot {self.sysname} {(hy.repr self.bind)})")
	 (impulse-plot self.S self.bind))

     True
     (print "Sorry, I do not understand what you want to plot. Try `help plot`")))

  (defn default [self line]
    (when (in line ["q"
		    "quit"
		    "bye"
		    "EOF"])
      (return True))

    (let [wl (.split line)]
      (cond
       (= [] wl)
       (return (. (super) (default line)))

       ;; ways that one might want to plot
       ;; without using the word "plot"
       (in (get wl 0) ["impulse" "step" "ramp" "response" "draw"])
       (self.do-plot line)

       ;; transfer function, without the word "transfer"
       (in (get wl 0) ["function"])
       (self.do-transfer line)

       ;; ways that one might want to bind/unbind variables
       ;; without using the word "bind" or "unbind"
       (in (get wl 0) ["set" "fix" "bind"])
       (self.do-bind line)

       (in (get wl 0) ["unset" "release"])
       (self.do-unbind line)

       ;; ways that one might ask for the RoC, stability
       (in (get wl 0) ["region" "convergence"])
       (do
	   (print f"(region-of-convergence {self.sysname} {(hy.repr self.bind)})")
	   (print (region-of-convergence self.S self.bind)))

       (in (get wl 0) ["stable" "stability"])
       (do
	   (print f"(stable {self.sysname} {(hy.repr self.bind)})")
	   (print (stable self.S self.bind)))

       ;; poles, zeros
       (in (get wl 0) ["pole"])
       (self.do-poles line)

       (in (get wl 0) ["zero"])
       (self.do-zeros line)

       ;; if no match
       ;; cut it, and repeat with rest of sentence
       (= (len wl) 1)
       (return (. (super) (default line)))

       True
       (do
	   (print (. " " (join (cut wl 1 None))))
	   (self.onecmd (. " " (join (cut wl 1 None)))))))))

;; TODO: figure out how to get the sysname automatically
(defmacro ask-ltiza [S [sysname "S"] [bind {}]]
  `(. (LTIza ~S ~sysname ~bind) (cmdloop)))
