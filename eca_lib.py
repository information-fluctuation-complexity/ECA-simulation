import sys, math, random, sympy

# bin_to_int(): returns the integer corresponding to bits list: [1, 0] -> 2
def bin_to_int(bits):
  return sum(x << i for i, x in enumerate(reversed(bits)))

# int_to_bin(): returns the bit list of length len corresponding to integer x
#   padded with leading zeros. e.g., for len == 8: 2 -> [0, 0, 0, 0, 0, 0, 1, 0]
def int_to_bin(x, len):
  return [x >> i & 1 for i in range(len - 1, -1, -1)]

# create_matrix(): returns an ECA transition matrix given rule number, length;
#   each element of the matrix is an integer representing one of the 2**length
#   ECA states
def create_matrix(rule_number, length):
  rule = int_to_bin(rule_number, 8)
  rule.reverse() # LSB of RULE_NUMBER is first in rule list
  n = 1 << length # number of states
  matrix = [[] for i in range(n)]
  for i in range(n):
    eca = int_to_bin(i, length)
    new_eca = [0] * length # temporary list for new cell values
    for left in [0, 1]:
      for right in [0, 1]:
        for j in range(length): # compute new cell values
          l = left if j == 0 else eca[j - 1]
          c = eca[j]
          r = right if j == length - 1 else eca[j + 1]
          new_eca[j] = rule[bin_to_int([l, c, r])]
        matrix[i].append(bin_to_int(new_eca))
  return(matrix)

# find_attractors(): given a transition matrix, returns a list of attractors
#   (each of which is a list of states identified by integers) by finding
#   isolated subsets of the states that are interconnected only among
#   themselves, in the steady-state, after eliminating transient states
def find_attractors(matrix):
  sys.setrecursionlimit(10**5)
  numStates = len(matrix)

  # state flags (some also used as search return values)
  INIT = 0
  ORIGIN = 1
  UNVISITED = 2
  VISITED = 3
  TRANSIENT = 4

  attractors = []
  attractor = set() # to avoid duplicates
  flags = [UNVISITED] * numStates
  collect = True # controls behavior of search

  # recursive search for origin state
  # returns ORIGIN (found), VISITED (not found), TRANSIENT (origin is transient)
  def search(s):
    if collect: attractor.add(s) # collect attractor in complete search
    f = flags[s]
    if f == INIT: flags[s] = ORIGIN # skip origin test at beginning of search
    elif f != UNVISITED: return(f)
    else: flags[s] = VISITED
    found = VISITED
    for i in set(matrix[s]): # remove duplicates with set(); search children
      c = search(i)
      if c == TRANSIENT: return(c)
      if c == ORIGIN:
        found = ORIGIN # keep track of having found origin anywhere
        if not collect: break # terminate search
    return(found)

  def reset(a):
    for i in a:
      if flags[i] == VISITED: flags[i] = UNVISITED

  for s in range(numStates):
    if flags[s] == TRANSIENT: continue
    attractor = set()
    flags[s] = INIT
    collect = True # do complete search, collect attractor states
    c = search(s)
    reset(attractor)
    if c != ORIGIN:
      flags[s] = TRANSIENT
      continue

    # origin was found again starting from origin; there is a return path
    # can origin can be found from any state in attractor?
    collect = False
    for i in attractor:
      c = search(i)
      reset(attractor)
      if c != ORIGIN:
        flags[s] = TRANSIENT
        break
    if flags[s] == TRANSIENT: continue # cul-te-sac
    attractors.append(list(attractor)) # yes, save it
    for i in attractor: flags[i] = TRANSIENT # treat as transient next iteration
  return(attractors)

# estimate_probabilities(): returns estimated state probabilities
#   after running a simulation with random stimulation, given transition matrix,
#   attractor list and number of iterations 
def estimate_probabilities(matrix, attractor, iterations):
  numStates = len(attractor)
  state = attractor[0] # start state
  probabilities = [0] * numStates # track state visits during simulation
  random.seed()
  for i in range(iterations): # run the automaton
    state = matrix[state][random.randint(0, 3)]
    probabilities[attractor.index(state)] += 1 # keep statistics
  total = sum(probabilities)
  for i in range(numStates): probabilities[i] /= total # normalize
  return(probabilities)

# compute_probabilities(): returns state probabilities given transition matrix
#   and attractor list; probabilites (synced with attractor) are returned as a
#   list of floats if ratios == False, else as a list of fractions as strings 
def compute_probabilities(matrix, attractor, ratios = False):
  length = len(attractor)
  width = len(matrix[0])
  X = [[0] * length for i in range(length + 1)]
  for i, a in enumerate(attractor):
    for j, b in enumerate(attractor): X[i][j] += matrix[b].count(a)
    X[i][i] -= width
  X[length] = [1] * length
  x = [0] * length; x.append(1)
  syms = sympy.symbols(' '.join(['p' + str(i) for i in range(length)]))
  syms = list(syms) if isinstance(syms, tuple) else [syms]
  probability_ratios = sympy.linsolve((sympy.Matrix(X), sympy.Matrix(x)), syms)
  probability_ratios = list(probability_ratios)[0] # convert to a tuple
  if ratios: return([str(prob) for prob in probability_ratios])
  return([prob.evalf() for prob in probability_ratios])

# compute_complexity(): returns (entropy, complexity) given probabilities list;
#   state information fluctuation complexity
def compute_complexity(probabilities):
  entropy = 0
  complexity = 0
  for i in range(len(probabilities)):
    if probabilities[i] == 0: continue
    information = -math.log2(probabilities[i])
    entropy += probabilities[i] * information
    complexity += probabilities[i] * information * information
  return (entropy, math.sqrt(complexity - entropy * entropy))

# compute_gamma_complexity(): returns (entropy, gammaComplexity) given
#   transition matrix, attractor and probabilities lists (synchronized);
#   net-information-gain fluctuation complexity
def compute_gamma_complexity(matrix, attractor, probabilities):
  width = len(matrix[0])
  entropy = 0
  gammaComplexity = 0
  for i, a in enumerate(attractor):
    if probabilities[i] == 0: continue
    infoPresent = -math.log2(probabilities[i])
    entropy += probabilities[i] * infoPresent
    nextStates = matrix[a]
    loop = [[attractor.index(b), nextStates.count(b)] for b in set(nextStates)]
    for j, num in loop:
      if i == j or probabilities[j] == 0: continue
      infoNext = -math.log2(probabilities[j])
      infoGain =  infoNext - infoPresent
      gammaComplexity += probabilities[i] * infoGain * infoGain * num / width
  return (entropy, math.sqrt(gammaComplexity))
