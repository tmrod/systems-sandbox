# systems-sandbox

A domain specific language (DSL) for learning linear systems analysis.

*Note:* The contents of this repo are under active development, and thus should not yet be treated as stable.
Function behavior, syntax, and other core features may change over time.
Use at one's own risk!

I also welcome comments, pull requests, etc.
In particular, I would find suggestions for new features and ideas about language syntax helpful.
This DSL is implemented in a LISP, so defining new syntax is not too difficult!

## Information

A DSL implemented using [Hylang](https://hylang.org/) (which is itself embedded in Python) for learning about the analysis of simple linear time-invariant systems.
The goal of this project is to provide a wrapper that allows students to compute properties of linear time-invariant systems without worrying about learning the messy details of a particular programming language.

There are three components to this project.
The first is the DSL, which is implemented in [systems_sandbox.hy](systems_sandbox.hy).
An example program written in the DSL is below.
```
(declare S (rational 1 [] [-3 -2 -1]))
(print (poles S))
(print (zeros S))
(pole-zero-plot S)

(declare S1 feedback S (compose (gain 2) (integrate)))
(pole-zero-plot S1)
```
The DSL supports symbolic computation (i.e., specifying systems with variable parameters), plotting, and basic functionality for asking questions about systems in terms of their Laplace transform.
Systems can be named, such as in the example above, as well as combined in series, parallel, and feedback configurations.

The next part is the REPL.
The REPL allows users to use the DSL interactively, much like other languages such as Python or Bash.

Finally, there is a "REPL inside of the REPL" called `LTIza`, which is a natural language chatbot designed to answer simple queries about systems specified using the DSL in a highly interactive way.
The goal of this chatbot is to allow the programming language to get out of your way, while also teaching you the formal commands used to duplicate the produced results.
An example REPL session using the chatbot is shown below.
```
=> (declare S (rational 1 [a] [-3 -2 b]))
=> (ask-ltiza S)
LTIza {} => what are the poles
[-3 -2 b]
LTIza {} => set a equal to 1 and b equal to -1
LTIza {a:1, b:-1} => what are the poles
[-3 -2 -1]
LTIza {a:1, b:-1} => plot the impulse response
LTIza {a:1, b:-1} => unbind the variable b
LTIza {a:1} => is the system stable
re(b) < 0
LTIza {a:1} => quit
=>
```
## Installation

See the appendix of the [manual](manual/manual.pdf).

## Usage

Read the [manual](manual/manual.pdf).
