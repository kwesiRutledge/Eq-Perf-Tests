# Automaton Experiments

The experiments that were performed, related to Automata and Simplification Algorithms, are described by this document.

## Usage

To run Automaton Experiment Set #4 (automata_t4.m). Enter the Eq-Perf-Tests outermost directory and call automata_tests().

```
results = automata_tests(4)
```

## Experiment List

1. Testing the use of the Finite State Machine formalism from the Automatica article by De Santis et. al. (2017). Implementing the basic operators and methods described in that work related to k-backward indistinguishable sets. Used a simple automaton for testing, hereon referred to as Automaton #1.

2. Tests identical to #1. Used a different, trickier automaton which will hereafter referred to as Automaton #2.

3. Testing the use of the Finite State Automata formal structure developed by Wang et. al. in a 2007 article mentioning an 'observable automaton.' Some FSA functions implemented but not all. Defined Automaton #1 in this new context.

4. Created an algorithm to create an automaton of the form that Prof. Ozay described, without using any previous work as guidance. Introduced Automaton 3. This is Prof. Ozay's example from a meeting on May 17.