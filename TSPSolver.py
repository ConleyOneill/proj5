#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import itertools



class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''
	
	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''
	def greedy( self,time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		start_city_index = 0
		while start_city_index <= ncities-1 and time.time()-start_time < time_allowance:
			start_city = cities[start_city_index]
			route = [start_city]

			current_city = start_city
			if(start_city_index < ncities-1):
				next_city = cities[start_city_index + 1]
			keep_looking = True
			while len(route) < ncities and time.time()-start_time < time_allowance and keep_looking:
				for city in cities:
					if(city not in route and current_city.costTo(city) != np.inf):
						if(current_city.costTo(next_city) > current_city.costTo(city)):
							next_city = city
				if(current_city.costTo(next_city) == np.inf):
					keep_looking = False
				else:
					route.append(next_city)
					current_city = next_city
			if route[-1].costTo(route[0]) != np.inf and bssf == None and len(route) == ncities:
				bssf = TSPSolution(route)
				count += 1
				foundTour = True
			elif route[-1].costTo(route[0]) != np.inf and len(route) == ncities:
				newBssf = TSPSolution(route)
				if newBssf.cost < bssf.cost:
					bssf = newBssf
				count += 1
				foundTour = True
			start_city_index += 1
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		results['path'] = route
		return results
	
	
	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''
		
	def branchAndBound( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		start_time = time.time()
		bssf = self.greedy()
		bssf_soln = bssf['soln']
		bssf_cost = bssf['cost']
		bssf_path = bssf['path']
		# Creating our cost matrix
		cost_matrix = []
		for city in cities:
			temp_matrix = []
			for target_city in cities:
				temp_matrix.append(city.costTo(target_city))
			cost_matrix.append(temp_matrix)
		#print(cost_matrix)
		#print('\r\n')
		lower_bound = 0
		lower_bound, cost_matrix = self.findLowerBoundReduceMatrix(lower_bound,cost_matrix)	# gives us our initial lower bound and reduces the cost matrix
		start_city = cities[0]	# We will always start at the first city in the array
		visited_cities = [start_city]
		pq = [(lower_bound, cost_matrix, visited_cities)]
		heapq.heapify(pq)
		#print(cost_matrix)
		#print(lower_bound)
		while len(pq) != 0:
			problem = heapq.heappop(pq)
			visited_cities = problem[2]
			current_city = problem[2][-1]
			#for city in cities:
			for i in range(len(cities)):
				if(current_city.costTo(cities[i]) != np.inf and cities[i] not in visited_cities):
					visited_cities.append(cities[i])
					for j in range(len(cost_matrix)):
						cost_matrix[i][j] = np.inf
						cost_matrix[j][i] = np.inf
					new_lower_bound,new_cost_matrix = self.findLowerBoundReduceMatrix(lower_bound,cost_matrix)
					if len(visited_cities) == len(cities):
						if new_lower_bound < bssf_cost:
							bssf_cost = new_lower_bound
							#bssf_path = visited_cities
							bssf_soln = TSPSolution(visited_cities)
					if new_lower_bound < bssf_cost:
						heapq.heappush(pq, (new_lower_bound, new_cost_matrix, visited_cities))
		end_time = time.time()

		results['cost'] = bssf_cost 
		results['time'] = end_time - start_time
		results['count'] = 0#count
		results['soln'] = bssf_soln
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		#results['path'] = bssf_path
		return results


	def findLowerBoundReduceMatrix(self, lower_bound, cost_matrix):
		lowest_row_cost_list = []
		for i in range(len(cost_matrix)):
			lowest_row_cost = np.inf
			for j in range(len(cost_matrix[i])):
				if cost_matrix[i][j] < lowest_row_cost:
					lowest_row_cost = cost_matrix[i][j]	# Finding the lowest cost for each row
			for j in range(len(cost_matrix[i])):
				cost_matrix[i][j] -= lowest_row_cost
			lowest_row_cost_list.append(lowest_row_cost)
			lower_bound += lowest_row_cost	# Updating the lower bound
		
		lowest_column_cost_list = []	# Reducing the columns
		for i in range(len(cost_matrix)):
			lowest_column_cost = np.inf
			for j in range(len(cost_matrix[i])):
				if cost_matrix[j][i] < lowest_column_cost:
					lowest_column_cost = cost_matrix[j][i]	# Reducing the rows
			#print(cost_matrix)
			lowest_column_cost_list.append(lowest_column_cost)
			lower_bound += lowest_column_cost	# Updating the lower bound
			for j in range(len(cost_matrix[i])):	# Reducing the columns
				cost_matrix[j][i] -= lowest_column_cost
		return lower_bound, cost_matrix







	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''
		
	def fancy( self,time_allowance=60.0 ):
		pass
		
class PriorityQueue:
	def __init__( self):
		list = []
		heapq.heapify(list)
	
	def pop(self):
		return heapq.heappop(list)
	
	def push(self, element):
		heapq.heappush(list, element)
	
	#def compare(self, element):
	#	pass


