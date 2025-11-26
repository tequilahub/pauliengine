import random
import time
import matplotlib.pyplot as plt
import numpy as np
import PauliEngine
from openfermion import QubitOperator
from tequila import QubitHamiltonian
from pauliarray import PauliArray
import pennylane as qml




import numpy as np
from pauliarray import PauliArray, WeightedPauliArray, Operator



def generate_pauliarray_hamiltonian(size, pauli_length):


    coeffs = np.random.choice([1.0, -1.0, 1j, -1j], size=size)
    z = np.random.randint(0, 2, size=(size, pauli_length), dtype=bool)
    x = np.random.randint(0, 2, size=(size, pauli_length), dtype=bool)
    paulis = PauliArray(z, x)
    return Operator(WeightedPauliArray(paulis, coeffs))




def multiply_pauliarrays(a, b):
    '''
    print("####A####")
    print(a.inspect())
    print("####B####")
    print(b.inspect())
    '''
    result = a * b
    '''
    print("RESULT")
    print(result.inspect())
    '''
    return result


def generate_random_paulistring(length):
    paulis = ["X", "Y", "Z", "I"]
    result = []
    for i in range(length):
        pauli = random.choice(paulis)
        if pauli != "I":
            result.append((i, pauli))
    return tuple(result)

def generate_random_qubit_hamiltonian(size, pauli_length):
    coeffs = [1.0, -1.0, 1j, -1j]
    qubit_op = QubitOperator()
    for _ in range(size):
        key = generate_random_paulistring(pauli_length)
        if key:
            qubit_op += QubitOperator(key, random.choice(coeffs))
    return QubitHamiltonian(qubit_operator=qubit_op)

def generate_paulistring(length):
    result = {}
    paulis = ["X", "Y", "Z", "I"]
    for i in range(length):
        current = random.choice(paulis)
        if current != "I":
            result[i] = current
    return result

def generate_Hamiltonian(size, pauli_length):
    coeff = [1.0, -1.0, 1j, -1j]
    result = []
    random.seed(42)
    for _ in range(size):
        result.append((random.choice(coeff), generate_paulistring(pauli_length)))
    return PauliEngine.QubitHamiltonian(result)




def generate_pennylane_hamiltonian(size, pauli_length):
    coeffs = [1.0, -1.0, 1j, -1j]
    ops = []
    
    for _ in range(size):
        pauli_dict = generate_paulistring(pauli_length)

        if len(pauli_dict) > 0:
            first_wire = next(iter(pauli_dict.keys()))
            op = qml.Identity(first_wire)
        else:
            op = qml.Identity(0)

        for index, p in pauli_dict.items():
            if p == "X": 
                op = op @ qml.PauliX(index)
            elif p == "Y": 
                op = op @ qml.PauliY(index)
            elif p == "Z": 
                op = op @ qml.PauliZ(index)

        ops.append(random.choice(coeffs) * op)

    H = ops[0]
    for o in ops[1:]:
        H = H + o
    return H




testing = []
result_pauli = []
result_ham = []

def benchmark_operation(operation_fn, *args, **kwargs):
    start = time.time()
    result = operation_fn(*args, **kwargs)
    '''
    if isinstance(result, WeightedPauliArray):
        print(result.inspect())
    
    elif isinstance(result, PauliEngine.QubitHamiltonian):
        print(result.to_string())
    '''
    return time.time() - start

def benchmark_vs_size(sizes, fixed_other, generator_fn, operation_fn, plot_name, fixed_type="pauli_length"):
    import time
    random.seed(42)
    times = []
    for size in sizes:
        start_total = time.time()
        

        start_create = time.time()
        if fixed_type == "pauli_length":
            h1 = generator_fn(size, fixed_other)
            h2 = generator_fn(size, fixed_other)
        else:
            h1 = generator_fn(fixed_other, size)
            h2 = generator_fn(fixed_other, size)
        
        end_create = time.time()
        
        # Zeit messen: Operation ausführen
        start_op = time.time()
        _ = operation_fn(h1, h2)
        end_op = time.time()
        
        end_total = time.time()
        
        create_time = end_create - start_create
        op_time = end_op - start_op
        total_time = end_total - start_total
        
        times.append(op_time)
        
        if plot_name not in testing:
            testing.append(plot_name)
        
        print(
            f"{plot_name}: {fixed_type} fix = {fixed_other}, var = {size} → "
            f"Erstellung: {create_time:.3f}s, Berechnung: {op_time:.3f}s, Gesamt: {total_time:.3f}s"
        )
    
    if fixed_type == "pauli_length":
        result_ham.append(times)
    else:
        result_pauli.append(times)
    return times



def plot_results(x_values, *curves, xlabel, title, labels=None):
    if labels is None:
        labels = [f"Variante {i+1}" for i in range(len(curves))]
    for curve, label in zip(curves, labels):
        plt.plot(x_values, curve, marker='o', label=label)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Zeit (s)")
    plt.grid(True)
    plt.legend()




hamiltonian_sizes = range(500, 5100, 500)
pauli_lengths = [100, 300, 500]
fixed_hamiltonian_size = 300
fixed_pauli_length = 300
gen_times = []

print("\nPauliEngine Benchmark")

times_engine_ham = benchmark_vs_size(
    hamiltonian_sizes,
    fixed_pauli_length,
    generate_Hamiltonian,
    lambda a, b: a * b,
    "PauliEngine",
    fixed_type="pauli_length"
)

times_engine_pauli = benchmark_vs_size(
    pauli_lengths,
    fixed_hamiltonian_size,
    generate_Hamiltonian,
    lambda a, b: a * b,
    "PauliEngine",
    fixed_type="hamiltonian_size"
)

'''
###################################################################
# RUN: PennyLane
###################################################################

print("\nPennyLane Benchmark")

times_pl_ham = benchmark_vs_size(
    hamiltonian_sizes,
    fixed_pauli_length,
    generate_pennylane_hamiltonian,
    lambda a, b: a @ b,
    "PennyLane",
    fixed_type="pauli_length"
)

times_pl_pauli = benchmark_vs_size(
    pauli_lengths,
    fixed_hamiltonian_size,
    generate_pennylane_hamiltonian,
    lambda a, b: a @ b,
    "PennyLane",
    fixed_type="hamiltonian_size"
)
'''
print("PauliArray")
times_pauliarray_ham = benchmark_vs_size(
    hamiltonian_sizes,
    fixed_pauli_length,
    generate_pauliarray_hamiltonian,
    multiply_pauliarrays,
    "PauliArray",
    fixed_type="pauli_length"
)

times_pauliarray_pauli = benchmark_vs_size(
    pauli_lengths,
    fixed_hamiltonian_size,
    generate_pauliarray_hamiltonian,
    multiply_pauliarrays,
    "PauliArray",
    fixed_type="hamiltonian_size"
)


'''
print("OpenFermion")
times_qubit_ham = benchmark_vs_size(
    hamiltonian_sizes,
    fixed_pauli_length,
    generate_random_qubit_hamiltonian,
    lambda a, b: a * b,
    "OpenFermion",
    fixed_type="pauli_length"
)

print("\n Benchmark: OpenFermion Pauli-Länge variiert")
times_qubit_pauli = benchmark_vs_size(
    pauli_lengths,
    fixed_hamiltonian_size,
    generate_random_qubit_hamiltonian,
    lambda a, b: a * b,
    "OpenFermion",
    fixed_type="hamiltonian_size"
)
'''
plt.figure(figsize=(14, 6))

plt.subplot(1, 2, 1)
plot_results(
    hamiltonian_sizes,
    *result_ham,
    xlabel="Hamiltonian-Size",
    title=f"Runtime vs Hamiltonian-Size (fix Pauli length = {fixed_pauli_length})",
    labels=testing
)
'''
plt.subplot(1, 2, 2)
plot_results(
    pauli_lengths,
    *result_pauli,
    xlabel="Pauli string length",
    title=f"Runtime vs Pauli string size (fix H size = {fixed_hamiltonian_size})",
    labels=testing
)
'''
plt.tight_layout()
plt.show()
