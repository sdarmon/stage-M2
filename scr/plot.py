import numpy as np 
import matplotlib.pyplot as plt
import sys


Arg = sys.argv[:]

if len(Arg) <3:
	print("Wrong input used.\n Use : "+Arg[0]+" output.txt format\n Where format is either 'dot' or 'his', the format of the plot.")

elif Arg[2] == "dot":
	y = []
	m = 0
	with open(Arg[1]) as f:
		for line in f:
			if (len(line)<2):
				break
			y.append(float(line[:-1]))
		x = np.arange(1,len(y)+1,1)

	xb = []
	yb = []

	xg = []
	yg = []

	for i in range(len(x)):
		if y[i] == float(0):
			yb.append(y[i])
			xb.append(x[i])
		else:
			yg.append(y[i])
			xg.append(x[i])

	fig, axs = plt.subplots(2)

	axs[0].plot(xg,yg, 'bo', label = 'Positive Frequency')
	axs[0].legend()
	axs[1].plot(xb,yb, 'ro', label = 'Null Frequency')
	axs[1].legend()
	axs[0].plot(x,y, 'green')
	axs[1].plot(x,y, 'green')
	axs[0].set(ylabel= "Frequencies (in %)")
	axs[1].set(ylabel= "Frequencies (in %)")
	plt.xlabel("Sizes of weight")
	fig.suptitle("Frequencies of the different sizes of weight (Case radius = 200)")
	plt.show()

elif Arg[2] == "his":
	y = []
	m = 0
	with open(Arg[1]) as f:
		for line in f:
			if (len(line)<2):
				break
			y.append(float(line[:-1]))
		x = np.arange(1,len(y)+1,1)

	Y=[0 for i in range(10)]
	X=[""for i in range(10)]
	for i in range(len(y)):
		Y[(10*i)//len(y)]+=y[i]
		X[(10*i)//len(y)] = i

	for i in range(9, 0,-1):
		X[i]=str(X[i-1]+1)+"->"+str(X[i])
	X[0] = "1->"+str(X[0])
	plt.bar(X[1:],Y[1:],color='green')
	plt.ylabel("Frequencies (in %)")
	plt.xlabel("Sizes of weight")
	plt.title("Histogram of the frequencies of the different sizes of weight (Case radius = 200)")
	for i, v in enumerate(Y):
		plt.text(i-1.25, v + .01, int(v*100*len(y)), color='b', fontweight='bold')
	plt.show()

elif Arg[2] == "top10":
	y=[]
	index = 0
	m = 0
	size = 1
	with open(Arg[1]) as f:
		for line in f:
			if (len(line)<2):
				break
		L = line.split('\t')
		y.append(int(L[2][:-1]))
		m = max(m,y[-1])
	x = np.arange(1,m+1,1)
	Y = [0 for i in range(m)]
	for el in y:
	    Y[el-1] +=100/len(y)
	M=0
	for el in Y:
		M+=el
		if M>=90:
			index = size
			break
		size+=1
	print(i)
	print(M)

else:
	print("Wrong format used.\n Use : "+Arg[0]+" output.txt format\n Where format is either 'dot' or 'his', the format of the plot.")