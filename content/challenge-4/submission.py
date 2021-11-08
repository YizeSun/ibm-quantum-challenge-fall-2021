"""
Challenge 4c solver function for the IBM Quantum Fall Challenge 2021
Author: Yize Sun <yize.sun@campus.lmu.de>
Score: 1,023,424

TODO:
1. Add author name, email and score at the top of the file
2. Write a summary of your approach in the header and highlight
the techniques that gave you the biggest improvement
3. Import all required libraries and modules
4. Print author name and score inside the `solver_function`
5. Copy code to each sub-functions (e.g.phase_return, subroutine_add_constant)
and explain the implementation in comments

Summary of the approach:
1. Integrate multiple Phase Gate(U1) into one by adding multiple phase parameters into one
2. **Incorporate multiple adding phases into one QFT pair**
3. Use different proper weights for Phase operator and Penalty dephasing

"""

from typing import List, Union
import math
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, assemble
from qiskit.compiler import transpile
from qiskit.circuit import Gate
from qiskit.circuit.library.standard_gates import *
from qiskit.circuit.library import QFT


def solver_function(L1: list, L2: list, C1: list, C2: list, C_max: int) -> QuantumCircuit:
    author = ' Yize Sun'
    score = 1023424
    print(f'{author}: {score}')

    # the number of qubits representing answers
    index_qubits = len(L1)

    # the maximum possible total cost
    max_c = sum([max(l0, l1) for l0, l1 in zip(C1, C2)])

    # the number of qubits representing data values can be defined using the maximum possible total cost as follows:
    data_qubits = math.ceil(math.log(max_c, 2)) + 1 if not max_c & (max_c - 1) == 0 else math.ceil(
        math.log(max_c, 2)) + 2

    ### Phase Operator ###

    # return part
    # so weights 
    ws = [0.2*0.78*math.pi, 0.2*1.02*math.pi, -0.2*1.02*math.pi]
    def phase_return(index_qubits: int, gamma: float, L1: list, L2: list, to_gate=True) -> Union[Gate, QuantumCircuit]:
        ##############################
        
        qr_index = QuantumRegister(index_qubits, "index")
        qc = QuantumCircuit(qr_index)
        # add return part for every line
        for i in range(index_qubits):
            qc.p(-ws[0]*gamma*(L2[i]-L1[i]), qr_index[i])
        # qc.decompose.draw()

        ##############################
        
        return qc.to_gate(label=" phase return ") if to_gate else qc
    
    # const_adder
    def const_adder(data_qubits: int, const: int, to_gate=True) -> Union[Gate, QuantumCircuit]:
        
        qr_data = QuantumRegister(data_qubits, "data")
        qc = QuantumCircuit(qr_data)

        ##############################
        ### QFT ###
        # use a effincienter adder later
        ##############################
        return qc.to_gate(label=" [ +" + str(const) + "] ") if to_gate else qc 
    
    def cost_calculation(index_qubits: int, data_qubits: int, list1: list, list2: list, to_gate = True) -> Union[Gate, QuantumCircuit]:
        qr_index = QuantumRegister(index_qubits, "index")
        qr_data = QuantumRegister(data_qubits, "data")
        qc = QuantumCircuit(qr_index, qr_data)

        # QFT for a + b
        qc.append(QFT(num_qubits=data_qubits, do_swaps=False).to_gate(), qr_data[:])

        # sum of intege               [{a1}, a2 ,.., {.}]
        # in different position       [ b1, {b2},..,  . ]
        for i, (val1, val2) in enumerate(zip(list1, list2)):
            ##############################
            ### Add val2 using const_adder controlled by i-th index register (set to 1) ###
            # Provide your code here
            const = val2
            m = 0 # sum of rotation
            pre = 1
            # check sign
            if const < 0:
                const = -const
                pre = -1

            const_b = format(const, f'0{data_qubits}b') # const in bitstr
            const_b = const_b[::-1] # in reversed formation
            # print(f'const_b: {const_b}')    
            # solution 1
            if const:
                for k in range(data_qubits): # i is 0, 1,..,n-1
                    for j in range(k+1): # j in [0],..,[0,..,data_qubits-1]
                        if int(const_b[j]): # b='00 000 000 01{0}'
                            m += 1/2**(k-j)
                    if (math.pi*m)%(2*math.pi) !=0:
                        qc.cp(pre*(math.pi*m)%(2*math.pi), qr_index[i], qr_data[k])
                        m=0

            qc.x(qr_index[i]) # flip control qubit

            ##############################
            ### Add val1 using const_adder controlled by i-th index register (set to 0) ###
            # Provide your code here
            const = val1
            if const < 0:
                const = -const
                pre = -1

            const_b = format(const, f'0{data_qubits}b')
            const_b = const_b[::-1]
            # print(f'const_b: {const_b}')    
            # solution 1
            if const:
                for k in range(data_qubits): # i is 0, 1,..,n-1
                    for j in range(k+1): # j in [0],..,[0,..,data_qubits-1]
                        if int(const_b[j]): # b='00 000 000 01{0}'
                            m += 1/2**(k-j)
                    if (math.pi*m)%(2*math.pi) !=0:
                        qc.cp(pre*(math.pi*m)%(2*math.pi), qr_index[i], qr_data[k])
                        m=0

            qc.x(qr_index[i]) # flip control qubit
        # adding QFT inverse gate to transform phase "a+b"to state of a+b
        qc.append(QFT(num_qubits=data_qubits, inverse=True, do_swaps=False).to_gate(), qr_data[:])

        return qc.to_gate(label=" Cost Calculation ") if to_gate else qc
    
    # constraint_testing
    def constraint_testing(data_qubits: int, C_max: int, to_gate = True) -> Union[Gate, QuantumCircuit]:
    
        qr_data = QuantumRegister(data_qubits, "data")
        qr_f = QuantumRegister(1, "flag")
        qc = QuantumCircuit(qr_data, qr_f)

        ##############################
        ### Set the flag register for indices with costs larger than C_max ###
        
        # QFT for adding
        qc.append(QFT(num_qubits=data_qubits, do_swaps=False).to_gate(), qr_data[:])
        # adding part works like above described
        c = (C_max).bit_length()
        gap = 2**c - C_max
        const = gap
        m = 0
        pre = 1
        if const < 0:
            const = -const
            pre = -1

        const_b = format(const, f'0{data_qubits}b')
        const_b = const_b[::-1]
        # print(f'const_b: {const_b}')    
        # solution 1
        if const:
            for k in range(data_qubits): # i is 0, 1,..,n-1
                for j in range(k+1): # j in [0],..,[0,..,data_qubits-1]
                    if int(const_b[j]): # b='00 000 000 01{0}'
                        m += 1/2**(k-j)
                if (math.pi*m)%(2*math.pi) !=0:
                    qc.p(pre*(math.pi*m)%(2*math.pi), qr_data[k])
                    m=0
        # QFT inverse=> adding ends
        # qc.append(QFT(num_qubits=data_qubits, inverse=True, do_swaps=False).to_gate(), qr_data[:])
        
        # cx gate (data_c->data_n-1, flag_reg)=(control, target)
        sub_cir = QuantumCircuit(1)
        sub_cir.x(0) # flip the gate
        qc.append(sub_cir.to_gate().control(num_ctrl_qubits=data_qubits-c, ctrl_state=0), qr_data[c:] + [qr_f[0]])

        # minus the padding part from cost
        # qc.append(QFT(num_qubits=data_qubits, do_swaps=False).to_gate(), qr_data[:])

        const = -gap
        m = 0
        pre = 1
        if const < 0:
            const = -const
            pre = -1

        const_b = format(const, f'0{data_qubits}b')
        const_b = const_b[::-1]
        # print(f'const_b: {const_b}')    
        # solution 1
        if const:
            for k in range(data_qubits): # i is 0, 1,..,n-1
                for j in range(k+1): # j in [0],..,[0,..,data_qubits-1]
                    if int(const_b[j]): # b='00 000 000 01{0}'
                        m += 1/2**(k-j)
                if (math.pi*m)%(2*math.pi) !=0:
                    qc.p(pre*(math.pi*m)%(2*math.pi), qr_data[k])
                    m=0
        qc.x(qr_f[0])
        # QFT inverse
        # qc.append(QFT(num_qubits=data_qubits, inverse=True, do_swaps=False).to_gate(), qr_data[:])
        ##############################

        return qc.to_gate(label=" Constraint Testing ") if to_gate else qc

    # penalty part
    def penalty_dephasing(data_qubits: int, alpha: float, gamma: float, to_gate = True) -> Union[Gate, QuantumCircuit]:
    
        qr_data = QuantumRegister(data_qubits, "data")
        qr_f = QuantumRegister(1, "flag")
        qc = QuantumCircuit(qr_data, qr_f)
        # c = (C_max).bit_length()

        ##############################
        ### Phase Rotation ###
        for i in range(data_qubits):
            qc.cp(ws[1]*(2**i)*alpha*gamma, qr_f[0], i)
        qc.p(ws[2]*C_max*alpha*gamma, qr_f[0])

        ##############################

        return qc.to_gate(label=" Penalty Dephasing ") if to_gate else qc
        
    def reinitialization(index_qubits: int, data_qubits: int, C1: list, C2: list, C_max: int, to_gate = True) -> Union[Gate, QuantumCircuit]:
        
        qr_index = QuantumRegister(index_qubits, "index")
        qr_data = QuantumRegister(data_qubits, "data")
        qr_f = QuantumRegister(1, "flag")
        qc = QuantumCircuit(qr_index, qr_data, qr_f)
        
        ##############################
        ### Reinitialization Circuit ###
        # use inverse method to inserve the gate
        qc.append(constraint_testing(data_qubits=data_qubits, C_max=C_max).inverse(), qr_data[:] + qr_f[:])
        qc.append(cost_calculation(index_qubits=index_qubits, data_qubits=data_qubits, list1=C1, list2=C2).inverse(), qr_index[:] + qr_data[:])
 
        ##############################
        
        return qc.to_gate(label=" Reinitialization ") if to_gate else qc

    ### Mixing Operator ###
    def mixing_operator(index_qubits: int, beta: float, to_gate = True) -> Union[Gate, QuantumCircuit]:
        
        qr_index = QuantumRegister(index_qubits, "index")
        qc = QuantumCircuit(qr_index)
        
        ##############################
        ### Mixing Operator ###
        
        for i in range(index_qubits):
            qc.rx(2*beta, qr_index[i])
        
        
        ##############################
        
        return qc.to_gate(label=" Mixing Operator ") if to_gate else qc

    qr_index = QuantumRegister(index_qubits, "index")  # index register
    qr_data = QuantumRegister(data_qubits, "data")  # data register
    qr_f = QuantumRegister(1, "flag")  # flag register
    cr_index = ClassicalRegister(index_qubits,
                                 "c_index")  # classical register storing the measurement result of index register
    qc = QuantumCircuit(qr_index, qr_data, qr_f, cr_index)

    ### initialize the index register with uniform superposition state ###
    qc.h(qr_index)

    ### DO NOT CHANGE THE CODE BELOW
    p = 5
    alpha = 1
    for i in range(p):
        ### set fixed parameters for each round ###
        beta = 1 - (i + 1) / p
        gamma = (i + 1) / p

        ### return part ###
        qc.append(phase_return(index_qubits, gamma, L1, L2), qr_index)

        ### step 1: cost calculation ###
        qc.append(cost_calculation(index_qubits, data_qubits, C1, C2), qr_index[:] + qr_data[:])

        ### step 2: Constraint testing ###
        qc.append(constraint_testing(data_qubits, C_max), qr_data[:] + qr_f[:])

        ### step 3: penalty dephasing ###
        qc.append(penalty_dephasing(data_qubits, alpha, gamma), qr_data[:] + qr_f[:])

        ### step 4: reinitialization ###
        qc.append(reinitialization(index_qubits, data_qubits, C1, C2, C_max), qr_index[:] + qr_data[:] + qr_f[:])

        ### mixing operator ###
        qc.append(mixing_operator(index_qubits, beta), qr_index)

    ### measure the index ###
    ### since the default measurement outcome is shown in big endian, it is necessary to reverse the classical bits in order to unify the endian ###
    qc.measure(qr_index, cr_index[::-1])

    return qc