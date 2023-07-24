
import sys
import os

filename = sys.argv[1]
print("start remove the last element", filename)

if not os.path.exists(filename):
    print("file doesn't exist")
    exit(1)

with open(filename, 'r', encoding='utf-8') as f:
    lines = f.readlines()
f.close()



f = open(filename, 'w') 
for line in lines:
    line = line.replace('\n','')
    data = line.split(' ')
    numOfline = len(data)
    for i in range(len(data)-2):
        print(data[i], file=f, end=" ")
    print(file=f)
print("finish remove the last element", filename)