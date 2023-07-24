import os


def encode_fidelity(f):
    return pow(f, 5) + 5 * pow(f, 4) * (1-f)

def fidelity(x):
    f = 0.94
    return pow(f, 10-x) * pow(encode_fidelity(f), x)
def waiting_time(x):
    return 4 * x + 10

directory_path = "encode_math/"
dataName = "encode_math.ans"
filename = directory_path + dataName

f = open(filename, 'w')

print(encode_fidelity(0.94))


for x in range(0, 11):
    print(x, fidelity(x), waiting_time(x), file=f)
