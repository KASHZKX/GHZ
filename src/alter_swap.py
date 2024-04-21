from math import ceil
import sys
import random

filename = sys.argv[1]
min_swap_prob = float(sys.argv[2])
max_swap_prob = float(sys.argv[3])

with open(filename, 'r', encoding='utf-8') as file1:
    file1_content = file1.readline()
    num_of_node = file1_content.split()[0]
    new_content = [file1_content.split()[0]]
    for index in range(0,int(num_of_node)):
        file1_content = file1.readline()
        file1_words = file1_content.split()
        file1_words[3] = str(min(1, random.random()*(max_swap_prob - min_swap_prob) + min_swap_prob))
        new_content.append(file1_words)
    
    file1_content = file1.readline()
    num_of_link = file1_content.split()[0]
    new_content.append(file1_content.split()[0])
    for index in range(0,int(num_of_link)):
        file1_content = file1.readline()
        file1_words = file1_content.split()
        new_content.append(file1_words)
    with open(filename, 'w', encoding='utf-8') as file2:
        for index2 in range(0,len(new_content)):
            if(type(new_content[index2]) == list):
                for index3 in range(0,len(new_content[index2])):
                    print(new_content[index2][index3],end=" ",file=file2)
            else:
                print(new_content[index2],end=" ",file=file2)
            print("",file=file2)

