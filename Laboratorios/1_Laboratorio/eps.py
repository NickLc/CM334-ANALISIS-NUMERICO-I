s=1/1000000000000000
for k in range(1,101):
	s=0.5*s
	t=s+1
	print("t= ",t)
	if t <= 1:
		s=s*2
		print("k-1= ",k-1," s= ",s )
		
