## Compute the complexity of an Elementary Cellular Automaton (ECA) by running a simulation or by exact calculation using the eca_lib.py python library

The python library [eca_lib.py](https://raw.githubusercontent.com/information-fluctuation-complexity/ECA-simulation/main/eca_lib.py) is a resource to accompany [Measuring complexity using information fluctuation: a tutorial](https://www.researchgate.net/publication/340284677) ([github copy](https://github.com/information-fluctuation-complexity/documents)), so that interested readers may verify the simulation results given in the tutorial for ECA rule 110, and simulate other ECA rules as desired. The library also includes a function to calculate the state probabilities of an ECA exactly, permitting a precise computation of complexity. It requires the installation of the python *sympy* module.

The eca_lib.py library contains the following functions:
* `bin_to_int(bits)`: returns the integer corresponding to *bits* list: e.g., [1, 0] → 2
* `int_to_bin(x, len)`: returns the bit list of length *len* corresponding to integer *x* padded with leading zeros. e.g., for *len == 4*: 2 → [0, 0, 1, 0]
* `create_matrix(rule_number, length)`: returns a matrix of possible transitions between the 2<sup>*length*</sup> ECA states, each state being identified by the integer equivalent of the string of bits contained in the automaton cells, given the ECA *rule_number* and *length*
* `find_attractors(matrix)`: given a transition *matrix*, returns a list of attractors (each of which is a list of states identified by integers) by finding isolated subsets of the states that are interconnected only among themselves, in the steady-state, after eliminating transient states
* `estimate_probabilities(matrix, attractor, iterations)`: returns a list of estimated state probabilities synchronized with *attractor* after running a simulation with random stimulation, given an ECA transition *matrix*, an *attractor* list and the number of *iterations*
* `compute_probabilities(matrix, attractor, ratios = False)`: returns a list of state probabilities synchronized with *attractor* given an ECA transition *matrix* and an *attractor* list, as computed exactly using `sympy.linsolve()`; returned probabilities are floats unless the optional argument *ratios* is *False*, in which case they are fractions as strings (for display)
* `compute_complexity(probabilities)`: returns (entropy, complexity) tuple given state *probabilities* list, where *complexity* is the primary complexity measure: *state information fluctuation complexity*
* `compute_gamma_complexity(matrix, attractor, probabilities)`: returns (entropy, gammaComplexity) tuple given *matrix, attractor, probabilities* lists, the latter two lists being synchronized, where *gammaComplexity* is the secondary complexity measure: *net-information-gain fluctuation complexity*

These functions can be used in the order given below:

### 1. Create a transition matrix

Give an ECA *rule_number* (such that 0 <= *rule_number* <= 255) and *length* (number of cells):

`create_matrix(rule_number, length)` returns a matrix of possible transitions between the 2<sup>*length*</sup> ECA states, each state being identified by the integer equivalent of the string of bits contained in the automaton cells.

### 2. Find attractors in state space

`find_attractors(matrix)` finds all attractors in the state space, meaning that it finds all isolated subsets of the states that are interconnected only among themselves, in the steady-state, after eliminating transient states, returned as a list of integer lists, states being identified by integers.

Each attractor, being an isolated subsystem, has its own particular value of complexity that must be determined separately. Unfortunately, much computing power is needed to find the attractors when the length of the ECA exceeds ~10, because the search requires exponentially more time and might also overflow the system stack because it uses a recursive search algorithm. Speed and stack capacity are system dependent; eca_lib.py increases stack capacity by a factor of 100 and on some systems it can be increased further. (Alternatively, an iterative search algorithm could be used.) However, most ECA rules have only one attractor and so the second step of finding attractors can be skipped. This table shows the number of ECA rules (out of 256) which have more than one attractor as a function of automaton length:

| ECA length | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|number rules with multiple attractors | 1 | 16 | 46 | 58 | 66 | 68 | 70 | 70 | 70 | 70 | 70 | 70 | 70 |

Apparently, if a rule has only one attractor when `LENGTH == 7` it will always have one attractor. ECA rule 110 has been confirmed to have only one attractor at least up to an automaton length of 19 cells.

For ECA rules with multiple attractors, instead of computing the complexity for each attractor separately, another approach is to compute the average complexity of a large number of simulations, each one starting from a random state. In this way, a single complexity value can represent a particular ECA rule, regardless of how many attractors it has.

### 3. Determine state probabilities

Given a transition matrix and a list of attractors, determine the state probabilities for each attractor using one of two methods:

#### 3.1. Determine approximate state probabilites via simulation with random stimulation

`estimate_probabilities(matrix, attractor, iterations)`: returns a list of estimated state probabilities synchronized with *attractor* after running a simulation with random stimulation, given an ECA transition *matrix*, an *attractor* list and the number of *iterations*

#### 3.2. Solve for state probabilites exactly using `sympy.linsolve()`

`compute_probabilities(matrix, attractor, ratios = False)`: returns a list of state probabilities synchronized with *attractor* given an ECA transition *matrix* and an *attractor* list, as computed exactly using `sympy.linsolve()`; returned probabilities are floats unless the optional argument *ratios* is *False*, in which case they are fractions as strings (for display)

### 4. Compute entropy and complexity

The complexity can be computed using two different formulas:

#### 4.1. Compute complexity using the primary complexity measure *state information fluctuation complexity* or simply *information fluctuation complexity* as it is known in the litererature

`compute_complexity(probabilities)`: returns (entropy, complexity) tuple given state *probabilities* list, where *complexity* is the primary complexity measure: *state information fluctuation complexity*.

#### 4.2.  Compute complexity using the secondary complexity measure *net-information-gain fluctuation complexity* or briefly, *gamma complexity*

`compute_gamma_complexity(matrix, attractor, probabilities)`: returns (entropy, gammaComplexity) tuple given *matrix, attractor, probabilities* lists, the latter two lists being synchronized, where *gammaComplexity* is the secondary complexity measure: *net-information-gain fluctuation complexity*.

### 5. Sample program

The sample program [eca_complexity.py](https://raw.githubusercontent.com/information-fluctuation-complexity/ECA-simulation/main/eca_complexity.py) computes information fluctuation complexity for an ECA rule using eca_lib.py.

Near the beginning of eca_complexity.py are five constants set to their default values:
```python
RULE_NUMBER = 110 # possible values: 0 to 255
LENGTH = 3 # number of cells in the automaton, > 0
ITERATION_POWER = 6 # number of simulation iterations (power of 10; 6 -> 10**6)
                    # if 0, do exact computation instead of simulation
SKIP_SEARCH = False # if True, skip attractor search; only one attractor exists
NET_INFORMATION_GAIN = False # if True, use the secondary complexity measure
                             # else use the primary complexity measure
```
These values can be modified manually or the program can be invoked from a shell using the following syntax:

`python eca_complexity.py -r 110 -l 3 -i 6 -s -n`

The options `-r, -l, -i, -s, -n` can override the default values of `RULE_NUMBER, LENGTH, ITERATION_POWER, SKIP_SEARCH, NET_INFORMATION_GAIN` respectively; the first three must each be followed by an integer.

Setting `SKIP_SEARCH = True` or including the command line option `-s` for rules with only one attractor results in a much faster determination of complexity for `LENGTH > 10` with little concern for stack overflow. It has been observed that for ECA rules with only one attractor,  *linsolve()* can correctly calculate the state probabilities even if it is given a list of all states, thus obviating the need to find the attractor beforehand.

