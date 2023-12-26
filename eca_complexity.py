# Computes information fluctuation complexity for elementary cellular automata

import sys, getopt, eca_lib as eca

RULE_NUMBER = 110 # possible values: 0 to 255
LENGTH = 3 # number of cells in the automaton, > 0
ITERATION_POWER = 7 # number of simulation iterations (power of 10; 7 -> 10**7)
                    # if 0, do exact computation instead of simulation
SKIP_SEARCH = False # if True, skip attractor search; only one attractor exists
NET_INFORMATION_GAIN = False # if True, use the secondary complexity measure
                             # else use the primary complexity measure
sys.tracebacklimit = 0
opts, args = getopt.getopt(sys.argv[1:], 'r:l:i:sn')
if len(args) != 0:  print('invalid argument'); sys.exit()
for o, a in opts:
  if o == '-r': RULE_NUMBER = int(a)
  elif o == '-l': LENGTH = int(a)
  elif o == '-i': ITERATION_POWER = int(a)
  else:
    SKIP_SEARCH = True if o == '-s' else False
    NET_INFORMATION_GAIN = True if o == '-n' else False
iterations = 0 if ITERATION_POWER <= 0 else 10**ITERATION_POWER
if not SKIP_SEARCH: sys.tracebacklimit = 0

print(f'Elementary cellular automaton complexity')
print(f'  rule {RULE_NUMBER:d}, {LENGTH:d} cells', end='')
if SKIP_SEARCH: print(f', skip sttractor search', end='')
if iterations == 0: print(f', exact solution', end='')
else: print(f', {iterations:d} iteration simulation', end='')
if NET_INFORMATION_GAIN: print(f', secondary', end='')
else: print(f', primary', end='')
print(f' complexity measure\n')

matrix = eca.create_matrix(RULE_NUMBER, LENGTH) # transition matrix
if SKIP_SEARCH: attractors = [list(range(1 << LENGTH))]
else: attractors = eca.find_attractors(matrix) # list of attractor lists
for attractor in attractors:
  attractor.sort() # optional, to start at the lowest state number
  if not SKIP_SEARCH:
    print(f'  Attractor of size {len(attractor):d}', end='')
    print(f' including state {attractor[0]:d}:')
  if iterations == 0: # exact solution
    probabilities = eca.compute_probabilities(matrix, attractor)
  else: # simulation using random simulation
    probabilities = eca.estimate_probabilities(matrix, attractor, iterations)
  if NET_INFORMATION_GAIN: # use secondary complexity measure
    (e, c) = eca.compute_gamma_complexity(matrix, attractor, probabilities)
  else: # use primary complexity measure
    (e, c) = eca.compute_complexity(probabilities)
  print(f'    Entropy = {e:.2f}, Complexity = {c:.2f}', end='')
  print(f', Complexity/Entropy = ', end='')
  print('indeterminate\n' if e == 0 else f'{c/e:.2f}\n')
