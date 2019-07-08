def f(x):
	r = (x*x+1)**(0.5)-1
	return r

def g(x):
	r = (x*x)/((x*x+1)**(0.5)+1)
	return r

print "N\tF(N)\tG(N)"
print "--------------------"
F=G=0

for i in range(1,100):
	#F = F + f(i)
	#G = F + g(i)
	F = f(8**(-1)*(i))
	G = g(8**(-1)*(i))
	print str(8**(-1)*(i))+"\t"+str(F)+"\t"+str(G)

