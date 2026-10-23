"""
Check one saved MAX identity using exact rational arithmetic.

The calculation base on the sorted chamber x_1 >= ... >= x_N

Passing proves the identity. The retained 'Invalid Identity' message means
this coefficient check failed;
"""

import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations



MAX_ERRORS = 10 # Just controls how many non-zero ambiguous residuals are printed after a failed check


def read_side(text):
    side = []
    
    for word in text.split():
        if len(word) != 2 or not word.isascii() or not word.isalpha():
            raise ValueError("Must have exactly 2 letters, such as 'aa' 'bc'")
        
        side.append( "".join(sorted(word)) ) # As max(x_a, x_b) = max(x_b, x_a)
    
    return tuple(sorted(side)) # Tuples of strings can be used as dictionary keys later



def read_identity(filename): 
    ''' Processes the input file and returns key infomation '''
    lines = []
    
    with open(filename, encoding="utf-8") as file:
        for line in file:
            if line:
                lines.append(line)

    if not lines:
        raise ValueError("File is empty")


    header = lines[0].split()
    if len(header) != 3 or header[2] != "=":
        raise ValueError("Header (1st) line must look like: 2 MAX_N = ...")
    

    name, n_text = header[1].split("_")
    
    if name != "MAX" or not n_text.isdigit(): raise ValueError("Not MAX_N in 1st line")
    
    
    target = Fraction(header[0])
    n = int(n_text)
    
    if n < 1 or target == 0: raise ValueError("N > 0 and the target != 0")


    k = None
    
    coefficients = {} # Dictionary, elem belike: (("aa","bc"), ("ab","cc")): Fraction(3, 2)
    
    for line in lines[1:]:
        
        # Example term: 13/42 * T[aa aa aa | aa aa aa]
        coefficient_text, pattern = line.split("*")
        coefficient = Fraction(coefficient_text.strip() )
        pattern = pattern.strip()
        
        if not pattern.startswith("T[") or not pattern.endswith("]"):
            raise ValueError("Each term must have the form: coeff * T[left | right]")
        
        left_text, right_text = pattern[2:-1].split("|")
        left = read_side(left_text)
        right = read_side(right_text)
        
        if len(left) != len(right) or not left:
            raise ValueError("Both sides must have the same number of words")
        
        if k is None: k = len(left)
        
        if len(left) != k:
            raise ValueError("Each term must have the same number of words per side.")
        if len(set("".join(left + right))) > n:
            raise ValueError("No term uses more than N distinct letters")


        # Swapping sides does not change the max by symmetry (combine repeated terms)
        if right < left:
            left, right = right, left
            
        pair = (left, right)
        if pair in coefficients: # coeff's 2/3 and -1/3 for same pattern become 1/3
            coefficients[pair] = coefficients[pair] + coefficient
        else:
            coefficients[pair] = Fraction(0) + coefficient


    terms = []
    for left, right in sorted(coefficients):
        coefficient = coefficients[left, right]
        if coefficient != 0:
            terms.append((coefficient, left, right))
    
    if not terms:
        raise ValueError("There are no nonzero terms to check")
    
    return n, k, target, terms
   
''' Notice: this function handles some symmetries: on words (m_ab = m_ba), 
    on sides (k=2: m_ab+m_cd = m_cd+m_ab), on patterns (swapping sides)
    But not equivalence under relabelling letters '''



def side_vector(side, assignment, m):
    ''' Turns a side into the coefficient vector of a linear expression '''
    vector = [0] * m
    
    for a, b in side:
        # Smaller rank means larger coordinate on the sorted chamber x_1 >= ... >= x_N
        rank = min(assignment[a], assignment[b]) # As max(x_a, x_b) = x_min(a, b)
        vector[rank] += 1
    
    return tuple(vector)



def dominates(left, right):
    ''' Related to ambiguous & unambiguous definition in paper,
    For example, x_1 + x_2 >= x_2 + x_3 always on sorted chamber '''
    # The two sides have equal total coefficient k.
    running = 0
    
    for i in range(len(left) - 1):
        running += left[i] - right[i]
        
        if running < 0:
            return False
    
    return True



def pair_order(pair):
    left, right = pair
    
    return tuple(reversed(left)), tuple(reversed(right))



def place_vector(vector, positions, n):
    # embeds an m-coordinate vector into an N-coordinate vector
    result = [0] * n
    
    for i, position in enumerate(positions):
        result[position] = vector[i]
    
    return tuple(result)



def count_pattern(left, right, n):
    ''' The central calculation for one atom A = max(\sum max's, \sum max's) '''
    letters = sorted(set("".join(left + right)))
    m = len(letters)
    
    left_pairs = []
    right_pairs = []
    
    # Convert letter names into integer indices/positions
    for a, b in left:
        left_pairs.append( (letters.index(a), letters.index(b)) )
    for a, b in right:
        right_pairs.append( (letters.index(a), letters.index(b)) )

    # First count relative orders of the letters actually used in this term
    counts = Counter()
    for assignment in permutations(range(m)):
    # Different assignments can led the same pair. Counter records how often that happens
        left_vector = side_vector(left_pairs, assignment, m)
        right_vector = side_vector(right_pairs, assignment, m)
        counts[(left_vector, right_vector)] += 1


    fixed = [0] * n
    ambiguous = Counter()
    
    for left_vector, right_vector in sorted(counts, key=pair_order):
        count = counts[left_vector, right_vector]
        
    # If one side dominates(unambiguous/fixed), its maximum becomes a linear expression 
    # Otherwise (ambiguous)
        if dominates(left_vector, right_vector):
            winner = left_vector
        elif dominates(right_vector, left_vector):
            winner = right_vector
        else:
            winner = None

        # Choose actual ranks for these letters. Together with the relative
        # orders above, this counts every injection exactly once.
        for positions in combinations(range(n), m):
            if winner is not None:
            # If unambiguous/fixed, adds its coefficients to the global linear vector
                for i, position in enumerate(positions):
                    fixed[position] += count * winner[i]
            else:
                first = place_vector(left_vector, positions, n)
                second = place_vector(right_vector, positions, n)
                
                if tuple(reversed(first)) > tuple(reversed(second)):
                    first, second = second, first
                
                ambiguous[first, second] += count
    
    return fixed, ambiguous



def print_result(n, target, fixed, ambiguous):
    entries = []
    
    for value in fixed:
        entries.append(str(value))
    print("\nFixed residual vector M * alpha - b = [" + ", ".join(entries) + "]")
    
    print(f"Nonzero ambiguous residual classes: {len(ambiguous):,}")

    # 'not any(fixed)' means every linear residual coefficient is 0
    # also we need no nonzero ambiguous coefficients remain: 'not ambiguous'
    valid = not any(fixed) and not ambiguous
    
    if valid:
        print("\nValid Identity")
        print(f"Exact verification passed: the displayed expression equals {target} MAX_{n}.")
    else:
        print("\nInvalid Identity")
        print("At least one exact rank-atom equation has a nonzero residual")
        
        if ambiguous and MAX_ERRORS > 0:
            largest = sorted(ambiguous.items(), key=lambda item: abs(item[1]), reverse=True)
            largest = largest[:MAX_ERRORS]
            
            print(f"\nFirst {len(largest)} ambiguous residuals:")
            
            for (left, right), value in largest:
                print(f"  {value} * ambiguous[{left} | {right}]")
    
    return valid



def verify(filename):
    n, k, target, terms = read_identity(filename)
    
    print(f"N: {n}", flush=True)
    print(f"Summands per side k: {k}", flush=True)
    
    fixed = [Fraction(0)] * n
    ambiguous = {}
    
    for coefficient, left, right in terms:
        linear_counts, ambiguous_counts = count_pattern(left, right, n)
        
        for i in range(n):
            fixed[i] += coefficient * linear_counts[i] # accumulated linear contribution
        
        for pair, count in ambiguous_counts.items():
            value = ambiguous.get(pair, Fraction(0)) + coefficient * count
            if value == 0:
                ambiguous.pop(pair, None)
            else:
                ambiguous[pair] = value
    fixed[0] -= target
    
    return print_result(n, target, fixed, ambiguous)



if __name__ == "__main__":
    try:
        if len(sys.argv) == 1:
            filename = input("Identity file: ").strip()
        elif len(sys.argv) == 2:
            filename = sys.argv[1]
        else:
            raise ValueError("Need input: python verify_MAX_identity.py identity_MAX_(N)_k(n).txt")
        valid = verify(filename)
    
    except (ValueError, ZeroDivisionError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
    
    raise SystemExit(0 if valid else 1)
