from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit import Aer, execute, IBMQ, BasicAer
from qiskit.aqua import QuantumInstance
import math
import math

def addition(circuit,x,y,c,N):
    for i in range(1,N):
        circuit.cx(x[i],y[i])
    circuit.cx(x[N-1],c)
    for i in range(N-2,0,-1):
        circuit.cx(x[i],x[i+1])
    for i in range(0,N-1):
        circuit.ccx(y[i],x[i],x[i+1])
    circuit.ccx(y[N-1],x[N-1],c)
    for i in range(N-1,0,-1):
        circuit.cx(x[i],y[i])
        circuit.ccx(y[i-1],x[i-1],x[i])
    for i in range(1,N-1):
        circuit.cx(x[i],x[i+1])
    for i in range(0,N):
        circuit.cx(x[i],y[i])

def subtraction(circuit,x,y,c,N):
    for i in range(N-1,-1,-1):
        circuit.cx(x[i],y[i])
    for i in range(N-2,0,-1):
        circuit.cx(x[i],x[i+1])
    for i in range(1,N):
        circuit.ccx(y[i-1],x[i-1],x[i])
        circuit.cx(x[i],y[i])
    circuit.ccx(y[N-1],x[N-1],c)
    for i in range(N-2,-1,-1):
        circuit.ccx(y[i],x[i],x[i+1])
    for i in range(1,N-1):
        circuit.cx(x[i],x[i+1])
    circuit.cx(x[N-1],c)
    for i in range(N-1,0,-1):
        circuit.cx(x[i],y[i])

def bitwise_complement(circuit,g,n):
    for i in range(0,n):
        circuit.x(g[i])

"""This function calculates just the carry bit of the addition between the value encoded into register a and the classical value c.
a in the register of n qubits to which we have to add the constant;
c is the constant;
g is the register of n-1 dirty ancilla qubits;
m is the qubit where we will save the carry."""
def carry_gate(circuit,a,g,m,c,n):
    circuit.cx(g[n-2],m)
    for i in range(n-1,1,-1):
        if(c[i]=='1'):
            circuit.cx(a[i],g[i-1])
            circuit.x(a[i])
        circuit.ccx(g[i-2],a[i],g[i-1])
    if(c[1]=='1'):
        circuit.cx(a[1],g[0])
        circuit.x(a[1])
    if(c[0]=='1'):
        circuit.ccx(a[0],a[1],g[0])
    for i in range(2,n):
        circuit.ccx(g[i-2],a[i],g[i-1])
    circuit.cx(g[n-2],m)
    for i in range(n-1,1,-1):
        circuit.ccx(g[i-2],a[i],g[i-1])
    if(c[0]=='1'):
        circuit.ccx(a[0],a[1],g[0])
    if(c[1]=='1'):
        circuit.x(a[1])
        circuit.cx(a[1],g[0])
    for i in range(2,n):
        circuit.ccx(g[i-2],a[i],g[i-1])
        if(c[i]=='1'):
            circuit.x(a[i])
            circuit.cx(a[i],g[i-1])
            
def controlled_carry_gate(circuit,ctrl1,a,g,m,c,n):
    circuit.ccx(ctrl1,g[n-2],m)
    for i in range(n-1,1,-1):
        if(c[i]=='1'):
            circuit.cx(a[i],g[i-1])
            circuit.x(a[i])
        circuit.ccx(g[i-2],a[i],g[i-1])
    if(c[1]=='1'):
        circuit.cx(a[1],g[0])
        circuit.x(a[1])
    if(c[0]=='1'):
        circuit.ccx(a[0],a[1],g[0])
    for i in range(2,n):
        circuit.ccx(g[i-2],a[i],g[i-1])
    circuit.ccx(ctrl1,g[n-2],m)
    for i in range(n-1,1,-1):
        circuit.ccx(g[i-2],a[i],g[i-1])
    if(c[0]=='1'):
        circuit.ccx(a[0],a[1],g[0])
    if(c[1]=='1'):
        circuit.x(a[1])
        circuit.cx(a[1],g[0])
    for i in range(2,n):
        circuit.ccx(g[i-2],a[i],g[i-1])
        if(c[i]=='1'):
            circuit.x(a[i])
            circuit.cx(a[i],g[i-1])
        
"""This circuit performs the increment of the value encoded into register x. We used the control qubit as the first alement of the register to increment, placing an X-gate on it at the end of the circuit, since:
-if it was 1, also the first qubit of the register to increment will toggle and so on;
-if it was 0, it will be the only one to be toggled.
The mathematical passages of the computation are: |x>|g> -> |x-g>|g> -> |x-g>|g'-1> -> |x-g-g'+1>|g'-1> -> |x+1>|g>, where:
-x is the register to increment with the control qubit as first element, so n qubits;
-g is a register with n-1 dirty qubits;
-g'-1 is the bitwise complement of g."""
def incremental_gate_controlled(circuit,x,g,n):
    if(n==2):
        circuit.cx(x[0],x[1])
    elif(n==3):
        circuit.ccx(x[0],x[1],x[2])
        circuit.cx(x[0],x[1])
    else:
        subtraction(circuit,g,x[:n-1],x[n-1],n-1)
        bitwise_complement(circuit,g,n-1)
        subtraction(circuit,g,x[:n-1],x[n-1],n-1)
        bitwise_complement(circuit,g,n-1)
        circuit.x(x[0])
        circuit.x(x[n-1])

        
"""The next function performs the addition of a constant value (c) to the value encoded in the quantum register (x) of n qubits, controlled on qubit ctrl. 
   The constant addition will need a dirty ancilla qubit g for the computation, that will be returned in the same state.
   The function is performed recursively.
   We do not need to control the incremental gate to obtain a controlled version of the constant addition present in 'ECDLP\ECC_library.py'."""
def controlled_constant_addition(circuit,ctrl,x,g,c,n):
    if(n==2):
        if(c[0]=='1'):
            circuit.ccx(ctrl,x[0],x[1])
        controlled_constant_addition(circuit,ctrl,x[:1],g,c[:1],1)
        controlled_constant_addition(circuit,ctrl,x[1:],g,c[1:],1)
    elif(n==1):
        if(c[0]=='1'):
            circuit.cx(ctrl,x[0])
    else:
        if(n%2==0):
            l=n//2
        else:
            l=n//2+1
        incremental_gate_controlled(circuit,g[0:1]+x[l:n],x[:l],n-l+1)
        for i in range(l,n):
            circuit.cx(g,x[i])
        controlled_carry_gate(circuit,ctrl,x[:l],x[l:],g,c[:l],l)
        incremental_gate_controlled(circuit,g[0:1]+x[l:n],x[:l],n-l+1)
        controlled_carry_gate(circuit,ctrl,x[:l],x[l:],g,c[:l],l)
        for i in range(l,n):
            circuit.cx(g,x[i])
        controlled_constant_addition(circuit,ctrl,x[:l],g,c[:l],l)
        controlled_constant_addition(circuit,ctrl,x[l:],g,c[l:],n-l)
    
def controlled_constant_subtraction(circuit,ctrl,x,g,c,n):
    sub_ctrl = QuantumRegister(1)
    sub_x = QuantumRegister(n)
    sub_g = QuantumRegister(1)
    sub_circ = QuantumCircuit(sub_ctrl,sub_x,sub_g, name='ctrl_adder')
    controlled_constant_addition(sub_circ,sub_ctrl,sub_x,sub_g,c,n)
    sub_inst = sub_circ.inverse().to_gate()
    circuit.append(sub_inst,ctrl[0:1]+x[0:n]+g[0:1])
    
def doubly_controlled_carry_gate(circuit,ctrl1,ctrl2,a,g,m,c,n):
    circuit.mcx(ctrl1[0:1]+ctrl2[0:1]+g[n-2:n-1],m)
    for i in range(n-1,1,-1):
        if(c[i]=='1'):
            circuit.cx(a[i],g[i-1])
            circuit.x(a[i])
        circuit.ccx(g[i-2],a[i],g[i-1])
    if(c[1]=='1'):
        circuit.cx(a[1],g[0])
        circuit.x(a[1])
    if(c[0]=='1'):
        circuit.ccx(a[0],a[1],g[0])
    for i in range(2,n):
        circuit.ccx(g[i-2],a[i],g[i-1])
    circuit.mcx(ctrl1[0:1]+ctrl2[0:1]+g[n-2:n-1],m)
    for i in range(n-1,1,-1):
        circuit.ccx(g[i-2],a[i],g[i-1])
    if(c[0]=='1'):
        circuit.ccx(a[0],a[1],g[0])
    if(c[1]=='1'):
        circuit.x(a[1])
        circuit.cx(a[1],g[0])
    for i in range(2,n):
        circuit.ccx(g[i-2],a[i],g[i-1])
        if(c[i]=='1'):
            circuit.x(a[i])
            circuit.cx(a[i],g[i-1])

"""This is the basic block of the modular addition needed, that is called by 'Shor_library.py'."""
def ccADDmodN(circuit, ctrl1, ctrl2, q, aux, g, a, N, n):
    s_a=bin(int(a))[2:].zfill(n)
    s_N=bin(int(N))[2:].zfill(n)
    N_a=bin(int(N)-int(a))[2:].zfill(n)
    """CMP(N-a)"""
    bitwise_complement(circuit,q,n)
    doubly_controlled_carry_gate(circuit,ctrl1,ctrl2,q,aux,g,N_a,n)
    bitwise_complement(circuit,q,n)
    """ADD(a)"""
    controlled_constant_addition(circuit,g[0:1],q,aux[0:1],s_a,n)
    """XOR"""
    circuit.ccx(ctrl1[0],ctrl2[0],aux[0])
    """SUB(N-a)"""
    controlled_constant_subtraction(circuit,g[0:1],q,aux[0:1],N_a,n)
    """XOR"""
    circuit.ccx(ctrl1[0],ctrl2[0],aux[0])
    """CMP(a)"""
    bitwise_complement(circuit,q,n)
    doubly_controlled_carry_gate(circuit,ctrl1,ctrl2,q,aux,g,s_a,n)
    bitwise_complement(circuit,q,n)

def ccADDmodN_inv(circuit, ctrl1, ctrl2, q, aux, g, a, N, n):
    sub_q = QuantumRegister(n)
    sub_ctrl1 = QuantumRegister(1)
    sub_ctrl2 = QuantumRegister(1)
    sub_aux = QuantumRegister(n-1)
    sub_g = QuantumRegister(1)
    sub_circ = QuantumCircuit(sub_ctrl1,sub_ctrl2,sub_q,sub_aux,sub_g, name='ccADDmodN')
    ccADDmodN(sub_circ,sub_ctrl1,sub_ctrl2,sub_q,sub_aux,sub_g, a, N, n)
    sub_inst = sub_circ.inverse().to_gate()
    circuit.append(sub_inst,ctrl1[0:1]+ctrl2[0:1]+q[0:n]+aux[0:n-1]+g[0:1])
