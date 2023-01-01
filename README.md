# Shor's algorithm Qiskit-based implementations

This repository provides a series of implementation of Shor's algorithms and its variants/optimisation to solve the prime factoring problem and the ECDLP according to diverse studies available in literature and other contributions from our side.

In particular, we provided the following implementations and testing scripts:

* Shor's algorithm for prime factoring with sequential QFT and in-place addition (Toffoli-based).
* Shor's algorithm adaptation to the ECDLP case: a library is provided with all the building block available. Clearly in this case testing the whole circuit on real hardware or simulators is not straightforward so we only provided the implementation and the validation of the constituent components.
* Testing scripts according metrics such as execution time, depth of the circuit, and quality of the solution.
* Testing scripts with noise models.
* Rigetti-based porting of the Shor's implementation for factoring.

We acknowledge Vito Medici for his contribution to this work.

