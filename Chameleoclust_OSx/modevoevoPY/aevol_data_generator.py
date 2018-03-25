import math 

def Aevol_objective_function(x):
	a = (1.2*math.exp(-1*((x-0.52)*(x-0.52))/(2*0.12*0.12)) - 1.4*math.exp(-1*((x-0.5)*(x-0.5))/(2*0.07*0.07)) + 0.3*math.exp(-1*((x-0.8)*(x-0.8))/(2*0.03*0.03)))
	if a > 0:
		return(a)
	else:
		return(0)	

class aevol_data:
	def __init__(self,function = Aevol_objective_function, size = 1000):
		self.size = size
		self.aevol_function = function
		self.data = []
		self.df = []
		self.creat_Aevol_objective_function()
		
	def creat_Aevol_objective_function(self):
		self.data=[]
		self.data.append([])
		i = 0
		self.data[0].append(-1)#label
		self.data[0].append(-1)#cluster
		while i < self.size:
			self.data[0].append([i,self.aevol_function(i*1./self.size)*self.size])
			i+=1
		self.df = self.data
		return self.data
	
	def table_to_data(self):
		return 0

	