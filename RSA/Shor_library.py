"""The functions in this file are adapted from https://github.com/ttlion/ShorAlgQiskit"""

from typing import Optional, Union, Tuple, List
import math
import array
import fractions
import logging
import numpy as np
import time
import random
import inplace_adder
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit import Aer, execute, IBMQ, BasicAer
from qiskit.aqua import QuantumInstance


""" Function to check if N is of type q^p"""
def check_if_power(N):
    """ Check if N is a perfect power in O(n^3) time, n=ceil(logN) """
    b=2
    while (2**b) <= N:
        a = 1
        c = N
        while (c-a) >= 2:
            m = int( (a+c)/2 )

            if (m**b) < (N+1):
                p = int( (m**b) )
            else:
                p = int(N+1)

            if int(p) == int(N):
                print('N is {0}^{1}'.format(int(m),int(b)) )
                return True

            if p<N:
                a = int(m)
            else:
                c = int(m)
        b=b+1

    return False


def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)


def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

""" Function to create QFT """
def create_QFT(circuit,up_reg,n,with_swaps):
    i=n-1
    """ Apply the H gates and Cphases"""
    """ The Cphases with |angle| < threshold are not created because they do 
    nothing. The threshold is put as being 0 so all CPhases are created,
    but the clause is there so if wanted just need to change the 0 of the
    if-clause to the desired value """
    while i>=0:
        circuit.h(up_reg[i])        
        j=i-1  
        while j>=0:
            if (np.pi)/(pow(2,(i-j))) > 0:
                circuit.cu1( (np.pi)/(pow(2,(i-j))) , up_reg[i] , up_reg[j] )
                j=j-1   
        i=i-1  
    """ If specified, apply the Swaps at the end """
    if with_swaps==1:
        i=0
        while i < ((n-1)/2):
            circuit.swap(up_reg[i], up_reg[n-1-i])
            i=i+1
            
""" Function to create inverse QFT """            
def create_inverse_QFT(circuit,up_reg,n,with_swaps):
    """ If specified, apply the Swaps at the beggining"""
    if with_swaps==1:
        i=0
        while i < ((n-1)/2):
            circuit.swap(up_reg[i], up_reg[n-1-i])
            i=i+1
    """ Apply the H gates and Cphases"""
    """ The Cphases with |angle| < threshold are not created because they do 
    nothing. The threshold is put as being 0 so all CPhases are created,
    but the clause is there so if wanted just need to change the 0 of the
    if-clause to the desired value """
    i=0
    while i<n:
        circuit.h(up_reg[i])
        if i != n-1:
            j=i+1
            y=i
            while y>=0:
                 if (np.pi)/(pow(2,(j-y))) > 0:
                    circuit.cu1( - (np.pi)/(pow(2,(j-y))) , up_reg[j] , up_reg[y] )
                    y=y-1   
        i=i+1
        
"""Function that calculates the angle of a phase shift in the sequential QFT based on the binary digits of a."""
"""a represents a possile value of the classical register."""
"""Needed only for sequential QFT version."""
def getAngle(a, N):
    """convert the number a to a binary string with length N"""
    s=bin(int(a))[2:].zfill(N) 
    angle = 0
    for i in range(0, N):
        """if the digit is 1, add the corresponding value to the angle"""
        if s[N-1-i] == '1': 
            angle += math.pow(2, -(N-i))
    angle *= np.pi
    return angle

"""Function that calculates the array of angles to be used in the addition in Fourier Space"""
def getAngles(a,N):
    s=bin(int(a))[2:].zfill(N) 
    angles=np.zeros([N])
    for i in range(0, N):
        for j in range(i,N):
            if s[j]=='1':
                angles[N-i-1]+=math.pow(2, -(j-i))
        angles[N-i-1]*=np.pi
    return angles

"""Creation of a doubly controlled phase gate"""
def ccphase(circuit,angle,ctl1,ctl2,tgt):
    circuit.cu1(angle/2,ctl1,tgt)
    circuit.cx(ctl2,ctl1)
    circuit.cu1(-angle/2,ctl1,tgt)
    circuit.cx(ctl2,ctl1)
    circuit.cu1(angle/2,ctl2,tgt)
    
"""Creation of the circuit that performs addition by a in Fourier Space"""
"""Can also be used for subtraction by setting the parameter inv to a value different from 0"""    
def phiADD(circuit,q,a,N,inv):
    angle=getAngles(a,N)
    for i in range(0,N):
        if inv==0:
            circuit.u1(angle[i],q[i])
        else:
            circuit.u1(-angle[i],q[i])
            
"""Single controlled version of the phiADD circuit"""
def cphiADD(circuit,q,ctl,a,n,inv):
    angle=getAngles(a,n)
    for i in range(0,n):
        if inv==0:
            circuit.cu1(angle[i],ctl,q[i])
        else:
            circuit.cu1(-angle[i],ctl,q[i])
         
"""Doubly controlled version of the phiADD circuit"""        
def ccphiADD(circuit,q,ctl1,ctl2,a,n,inv):
    angle=getAngles(a,n)
    for i in range(0,n):
        if inv==0:
            ccphase(circuit,angle[i],ctl1,ctl2,q[i])
        else:
            ccphase(circuit,-angle[i],ctl1,ctl2,q[i])
            
"""Circuit that implements doubly controlled modular addition by a"""            
def ccphiADDmodN(circuit, q, ctl1, ctl2, aux, a, N, n):
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 0)
    phiADD(circuit, q, N, n, 1)
    create_inverse_QFT(circuit, q, n, 0)
    circuit.cx(q[n-1],aux)
    create_QFT(circuit,q,n,0)
    cphiADD(circuit, q, aux, N, n, 0)
    
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 1)
    create_inverse_QFT(circuit, q, n, 0)
    circuit.x(q[n-1])
    circuit.cx(q[n-1], aux)
    circuit.x(q[n-1])
    create_QFT(circuit,q,n,0)
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 0)
    
"""Circuit that implements the inverse of doubly controlled modular addition by a"""    
def ccphiADDmodN_inv(circuit, q, ctl1, ctl2, aux, a, N, n):
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 1)
    create_inverse_QFT(circuit, q, n, 0)
    circuit.x(q[n-1])
    circuit.cx(q[n-1],aux)
    circuit.x(q[n-1])
    create_QFT(circuit, q, n, 0)
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 0)
    cphiADD(circuit, q, aux, N, n, 1)
    create_inverse_QFT(circuit, q, n, 0)
    circuit.cx(q[n-1], aux)
    create_QFT(circuit, q, n, 0)
    phiADD(circuit, q, N, n, 0)
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 1)
    
"""Circuit that implements the controlled modular multiplication by a."""    
def Ua_basic(circuit, ctl, q, aux, a, N, n):
    create_QFT(circuit,aux,n+1,0)
    for i in range(0, n):
        ccphiADDmodN(circuit, aux, q[i], ctl, aux[n+1], (2**i)*a % N, N, n+1)
    create_inverse_QFT(circuit, aux, n+1, 0)

    for i in range(0, n):
        circuit.cswap(ctl,q[i],aux[i])

    a_inv = modinv(a, N)
    create_QFT(circuit, aux, n+1, 0)
    i = n-1
    while i >= 0:
        ccphiADDmodN_inv(circuit, aux, q[i], ctl, aux[n+1], math.pow(2,i)*a_inv % N, N, n+1)
        i -= 1
    create_inverse_QFT(circuit, aux, n+1, 0)
    
"""Function for the in-place adder oracle"""
def Ua_inplace(circuit, ctl, q, aux, a, N, n):
    for i in range(0, n):
        #print(i)
        inplace_adder.ccADDmodN(circuit, ctl[0:1], q[i:i+1], aux[0:n], q[0:i]+q[i+1:n], aux[n:n+1], (2**i)*a % N, N, n)

    for i in range(0, n):
        circuit.cswap(ctl,q[i],aux[i])
    #print('--------')
    a_inv = modinv(a, N)
    i = n-1
    while i >= 0:
        #print(i)
        inplace_adder.ccADDmodN_inv(circuit, ctl[0:1], q[i:i+1], aux[0:n], q[0:i]+q[i+1:n], aux[n:n+1], (2**i)*a % N, N, n)
        i -= 1

"""Function to create and execute the circuit.
    Registers:
        -> aux - auxilliary quantum register used in addition and multiplication
        -> up_reg - quantum register where the QFT (sequential or not) is performed
        -> down_reg - quantum register where the multiplications are made
        -> up_classic - classical register where the measured values of the QFT are stored"""
def period_finder(N,a,number_shots,sequential,simulator,ibm,quantum_machine,method):
    n = math.ceil(math.log(N,2))
    if sequential=='False':
        if method=='Basic':
            n_qubit=4*n+2
        else:
            n_qubit=4*n+1
    else:
        if method=='Basic':
            n_qubit=2*n+3
        else:
            n_qubit=2*n+2
    #print('Total number of qubits used: {0}\n'.format(n_qubit))
    
    if method=='Basic':
        """initialize the registers and the circuit"""
        aux = QuantumRegister(n+2)
        if sequential=='False':
            up_reg = QuantumRegister(2*n)
        else:
            up_reg = QuantumRegister(1)
            c_aux = ClassicalRegister(1)
        down_reg = QuantumRegister(n)
        up_classic = ClassicalRegister(2*n)
        if sequential=='False':
            circuit = QuantumCircuit(down_reg , up_reg , aux, up_classic)
        else:
            circuit = QuantumCircuit(down_reg , up_reg , aux, up_classic,c_aux)
        """Create the circuit"""
        if sequential=='False':
            """ Initialize down register to 1 and create maximal superposition in top register """
            circuit.h(up_reg)
            circuit.x(down_reg[0])
            """ Apply the multiplication gates as showed in the report in order to create the exponentiation """
            for i in range(0, 2*n):
                Ua_basic(circuit, up_reg[i], down_reg, aux, int(pow(a, pow(2, i))), N, n)
            """ Apply inverse QFT """
            create_inverse_QFT(circuit, up_reg, 2*n ,1)
            """ Measure the top qubits, to get x value"""
            circuit.measure(up_reg,up_classic)
        else:
            """ Initialize down register to 1"""
            circuit.x(down_reg[0])
            """ Cycle to create the Sequential QFT, measuring qubits and applying the right gates according to measurements """
            for i in range(0, 2*n):
                """reset the top qubit to 0 if the previous measurement was 1"""
                circuit.x(up_reg).c_if(c_aux, 1)
                circuit.h(up_reg)
                Ua_basic(circuit, up_reg[0], down_reg, aux, a**(2**(2*n-1-i)), N, n)
                """cycle through all possible values of the classical register and apply the corresponding conditional phase shift"""
                for j in range(0, 2**i):
                    """the phase shift is applied if the value of the classical register matches j exactly"""
                    circuit.u1(getAngle(j, i), up_reg[0]).c_if(up_classic, j)
                circuit.h(up_reg)
                circuit.measure(up_reg[0], up_classic[i])
                circuit.measure(up_reg[0], c_aux[0])
    else:
        """Case with inplace_adder"""
        if sequential=='False':
            aux = QuantumRegister(n+1)                          
            up_reg = QuantumRegister(2*n)
            down_reg = QuantumRegister(n)
            up_classic = ClassicalRegister(2*n)
            circuit = QuantumCircuit(down_reg , up_reg , aux, up_classic)
            circuit.h(up_reg)
            circuit.x(down_reg[0])
            for i in range(0, 2*n):
                Ua_inplace(circuit, up_reg[i:i+1], down_reg, aux, int(pow(a, pow(2, i))), N, n)
            create_inverse_QFT(circuit, up_reg, 2*n ,1)
            circuit.measure(up_reg,up_classic)
        else:
            aux = QuantumRegister(n+1)                          
            up_reg = QuantumRegister(1)
            down_reg = QuantumRegister(n)
            up_classic = ClassicalRegister(2*n)
            c_aux = ClassicalRegister(1)
            circuit = QuantumCircuit(down_reg , up_reg , aux, up_classic,c_aux)
            circuit.x(down_reg[0])
            for i in range(0, 2*n):
                """reset the top qubit to 0 if the previous measurement was 1"""
                circuit.x(up_reg).c_if(c_aux, 1)
                circuit.h(up_reg)
                Ua_inplace(circuit, up_reg[0:1], down_reg, aux, a**(2**(2*n-1-i)), N, n)
                """cycle through all possible values of the classical register and apply the corresponding conditional phase shift"""
                for j in range(0, 2**i):
                    """the phase shift is applied if the value of the classical register matches j exactly"""
                    circuit.u1(getAngle(j, i), up_reg[0]).c_if(up_classic, j)
                circuit.h(up_reg)
                circuit.measure(up_reg[0], up_classic[i])
                circuit.measure(up_reg[0], c_aux[0])

        
    
    if simulator=='True':
        #print('Simulating the circuit {0} times for N={1} and a={2}\n'.format(number_shots,N,a))
        if ibm=='False':
            backend = BasicAer.get_backend('qasm_simulator')
        else:
            provider = IBMQ.get_provider('ibm-q')
            backend = provider.get_backend('ibmq_qasm_simulator')
        simulation = execute(circuit, backend=backend ,shots=number_shots)
        sim_result = simulation.result()
        counts_result = sim_result.get_counts(circuit)
        T = sim_result.time_taken
        #print("Time consumed in quantum working: ",T)
        return (sim_result, counts_result,T)
    else:    
        #print('Executing the circuit {0} times for N={1} and a={2} on {3}\n'.format(number_shots,N,a,quantum_machine))
        provider = IBMQ.get_provider('ibm-q')
        qcomp = provider.get_backend(quantum_machine)
        job = execute(circuit, backend=qcomp, shots=number_shots)
        result = job.result()
        counts_result = result.get_counts(circuit)
        T = result.time_taken
        #print("Time consumed in quantum working: ",T)
        return (result, counts_result,T)
    
    
def get_factors(x_value,N,a):
    t_upper = 2*(math.ceil(math.log(N,2)))
    if x_value<=0:
        return (0,0)
    """ Calculate T and x/T """
    T = pow(2,t_upper)
    x_over_T = x_value/T
    """ Cycle in which each iteration corresponds to putting one more term in the calculation of the Continued Fraction (CF) of x/T """
    """ Initialize the first values according to CF rule """
    i = 0
    b = array.array('i')
    t = array.array('f')
    b.append(math.floor(x_over_T))
    t.append(x_over_T - b[i])
    while i>=0:
        """From the 2nd iteration onwards, calculate the new terms of the CF based
        on the previous terms as the rule suggests"""
        if i>0:
            b.append( math.floor( 1 / (t[i-1]) ) )
            t.append( ( 1 / (t[i-1]) ) - b[i] )
        """print('i: {2} ; b: {0} ; t: {1}\n'.format(b[i],t[i],i))"""
        """ Calculate the CF using the known terms """
        aux = 0
        j=i
        while j>0:
            aux = 1 / ( b[j] + aux )
            j = j-1
        aux = aux + b[0]
        """Get the denominator from the value obtained"""
        frac = fractions.Fraction(aux).limit_denominator()
        r=frac.denominator
        """ Increment i for next iteration """
        i=i+1
        if i>=15:
            if (r%2) == 1:
                sqr=math.sqrt(a)
                if (sqr-int(sqr))!=0:
                    """ Return if already too much tries and numbers are huge """ 
                    return (0,0) 
                else:
                    a=int(sqr)
                    r=2*r         
        if (r%2) == 1:
            continue    
        """ If denominator even, try to get factors of N """
        """ Get the exponential a^(r/2) """
        exponential = a
        s = 1
        while s < r/2:
            exponential = (exponential * a) % N
            s=s+1
        """if r<1000:
            try:
                exponential=pow(a , (r/2))
            except:
                continue"""
        """ Check if the value is too big or not """
        if math.isinf(exponential)==1 or exponential>1000000000:
            return (0,0)
        """If the value is not too big (infinity), then get the right values and do the proper gcd()"""
        putting_plus = int(exponential + 1)
        putting_minus = int(exponential - 1)    
        one_factor = math.gcd(putting_plus,N)
        other_factor = math.gcd(putting_minus,N)    
        """ Check if the factors found are trivial factors or are the desired factors """
        if one_factor==1 or one_factor==N or other_factor==1 or other_factor==N:
            """ Check if the number has already been found, use i-1 because i was already incremented """
            if t[i-1]==0:
                return (0,0)     
            if i>=15:
                """ Return if already too much tries and numbers are huge """ 
                return (0,0)
        else:
            return (one_factor,other_factor)
        
"""Function to factorize a biprimal number
    Parameters:
        -> N - number to factorize
        -> a - random number useful for the computation
        -> number_shots - number of shots of the circuit
        -> sequential - true if the sequential QFT method must be used
        -> simulator - true if the circuit must be simulated, false if it must be executed on a real quantum device
        -> ibm - true if, in case of simulation, the ibm simulator must be used
        -> quantum_machine - name of the quantum machine that must be used in case the circuit must be executed on a real quantum device
        -> method - Basic if the method that performs the additions in the Fourier space is used, inplace if the method that uses Toffoli based additions is used
        -> runtime - True if the runtime of the circuit is wanted in the outputs
        -> prob - True if the probability of success is wanted in the outputs
    Outputs:
        -> p, q - factor of the number N (if founded) or (0, 0)
        -> T - runtime of the circuit if wanted
        -> P - probability of success if wanted"""
def Shor(N,a=0,number_shots=0,sequential='True',simulator='True',ibm='True',quantum_machine='',method='Basic',runtime='False',prob='False'):
    """ Check if N is even """
    if (N%2)==0:
        print('N is even, so does not make sense!')
        exit()    
    """ Check if N can be put in N=p^q, p>1, q>=2 """
    """ Try all numbers for p: from 2 to sqrt(N) """
    if check_if_power(N)==True:
        exit()    
    #print('Not an easy case, using the quantum circuit is necessary\n')
    if a==0:
        a = random.randint(2,N-1)
    if number_shots==0:
        number_shots=100
    if math.gcd(a,N)!=1:
        print('{0} is a factor of {1}'.format(math.gcd(a,N),N))
        exit()    
    """print('Starting quantum part')"""    
    (result,counts_result,T) = period_finder(N,a,number_shots,sequential,simulator,ibm,quantum_machine,method)
    """print('Quantum part finished')"""        
    """ Print info to user from the simulation results """
    """print('Printing the various results followed by how many times they happened (out of the {} cases):\n'.format(number_shots))"""
    """i=0
    while i < len(counts_result):
        print('Result \"{0}\" happened {1} times out of {2}'.format(list(result.get_counts().keys())[i],list(result.get_counts().values(                                                                   [i],number_shots))
        i=i+1"""
    """ An empty print just to have a good display in terminal """
    """print(' ')"""
    """ Initialize this variable """
    prob_success=0
    """ For each simulation result, print proper info to user and try to calculate the factors of N"""
    i=0
    pstar=0
    qstar=0
    while i < len(counts_result):
        """ Get the x_value from the final state qubits """
        output_desired = list(result.get_counts().keys())[i]
        if sequential=='True':
            output_desired=output_desired.split(" ")[1]
        x_value = int(output_desired, 2)
        prob_this_result = 100 * ( int( list(result.get_counts().values())[i] ) ) / (number_shots)
        """print("-----> Analysing result {0}. This result happened in {1:.4f} % of all cases\n".format(output_desired,prob_this_result))"""
        """ Print the final x_value to user """
        """print('In decimal, x_final value for this result is: {0}\n'.format(x_value))"""
        """ Get the factors using the x value obtained """   
        (p,q) = get_factors(int(x_value),int(N),int(a))
        """print('{0},{1}'.format(p,q))"""
        if p!=0 or q!=0:
            prob_success = prob_success + prob_this_result
            pstar=p
            qstar=q
        i=i+1
    P=prob_success
    #print("\nUsing a={0}, found the factors of N={1} in {2:.4f} % of the cases\n".format(a,N,P))
    if runtime!='False' and prob=='False':
        return (pstar,qstar,T)
    if runtime!='False' and prob!='False':
        return (pstar,qstar,T,P)
    if runtime=='False' and prob!='False':
        return (pstar,qstar,P)
    return (pstar,qstar)
