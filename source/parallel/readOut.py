# Using readlines() 
import re
import numpy as np
from  scipy import stats
import matplotlib.pyplot as plt
file1 = open('outputTest.txt', 'r') 
Lines = file1.readlines() 
  
mylist = [[] for i in range(26-8)]
#print(mylist)
counter = [0]*(26-8)
count = 0
# Strips the newline character 
for line in Lines: 
   # print("Line{}: {}".format(count, line.strip())) 

	x = re.search(".*_(\d+)_(\d+).*: (\d+)", line) 
	mylist[int(x.group(1))-8].append(int(x.group(3)))
	counter[int(x.group(1))-8]+=1

#for i in range(0,25-8):
#	print("Lost data "+str(i)+" "+str(100-counter[i]))
m=[]
std=[]
print("4 processor")
for l in mylist:
	nplist = np.array(l)
	m.append(nplist.mean())
	std.append(nplist.std())
	print("%.2f" %nplist.mean())


#print(mylist)
mylist = [[] for i in range(22-4)]

# Strips the newline character 
file1 = open('outputSingle.txt', 'r') 
Lines = file1.readlines() 
for line in Lines: 

	x = re.search(".*_(\d+)_(\d+).*: (\d+)", line) 

	mylist[int(x.group(1))-4].append(int(x.group(3)))


#for i in range(0,25-8):
#	print("Lost data "+str(i)+" "+str(100-counter[i]))
mSingle=[]
stdSingle=[]
print("One processor")
for l in mylist:
	nplist = np.array(l)
	mSingle.append(nplist.mean())
	stdSingle.append(nplist.std())
	print("%.2f" %nplist.mean())
plt.plot(range(4,22), mSingle,label="1 processor")
plt.plot(range(8,26), m,label="4 processor")
plt.xlabel("Problem size")
plt.ylabel("Execution time (ms)")
plt.legend()
plt.show()


plt.title("Elapsed time on number of processors")
plt.plot([1,2,4,8,16],[110505,45030,21415,14075,22665])
plt.xlabel("Number of processors")
plt.ylabel("Execution time (ms)")
plt.show()


