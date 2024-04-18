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
		bssf = self.greedy()	# Initially run the greedy approach to find the initiall bssf
		bssf_soln = bssf['soln']
		bssf_cost = bssf['cost']
		max_queue_size = 0
		pruned = 0
		count = 0
		solutions = 0
		
		cost_matrix_initial = []	# Creating our cost matrix
		for temp_city in cities:
			temp_matrix = []
			for target_city in cities:
				temp_matrix.append(temp_city.costTo(target_city))
			cost_matrix_initial.append(temp_matrix)
		pq = []
		heapq.heapify(pq)
		lower_bound = 0		# Finding the initial lower bound
		lower_bound, cost_matrix = self.findInitialLowerBoundReduceMatrix(lower_bound,cost_matrix_initial)	# gives us our initial lower bound and reduces the cost matrix
		start_city = cities[0]	# We will always start at the first city in the array
		visited_cities = [start_city]
		heapq.heappush(pq, (lower_bound, cost_matrix_initial, visited_cities))
		while len(pq) != 0:
			problem = heapq.heappop(pq)
			current_city = problem[2][-1]	# Our current city is always going to be the last element that we added to our visited cities
			visited_cities_copy = problem[2]
			current_city_index = cities.index(current_city)
			for i in range(len(cities)):	# This cycles through all of the available cities
				visited_cities = visited_cities_copy.copy()
				if(current_city.costTo(cities[i]) != np.inf and cities[i] not in visited_cities):	# If the cost to get to the next city isnt infinite and hasnt been visited, we visit
					visited_cities.append(cities[i])
					count += 1
					if len(pq) > max_queue_size:	# Updating the max queue size
						max_queue_size = len(pq)
					new_lower_bound,new_cost_matrix = self.findLowerBoundReduceMatrix(lower_bound,cost_matrix, current_city_index, i)	# Calculating the updated lower bound and cost matrix
					if len(visited_cities) == len(cities) and new_lower_bound < bssf_cost:
						bssf_cost = new_lower_bound
						bssf_soln = TSPSolution(visited_cities)	# Creating a new solution if we found one and updating our number of solutions
						solutions += 1
					elif new_lower_bound <= bssf_cost:
						heapq.heappush(pq, (new_lower_bound, new_cost_matrix, visited_cities))	# Pushing a new element into the pq
					else:
						pruned += 1	# Updates our pruned nodes
				else:	# If the cost to the city is infinite and it is not in the visited cities, we can prune it
					pruned += 1
		end_time = time.time()
		
		results['cost'] = bssf_soln.cost 
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf_soln
		results['max'] = max_queue_size
		results['total'] = solutions
		results['pruned'] = pruned
		return results


	def findInitialLowerBoundReduceMatrix(self, lower_bound, cost_matrix):
		for i in range(len(cost_matrix)):
			lowest_row_cost = np.inf
			for j in range(len(cost_matrix[i])):
				if cost_matrix[i][j] < lowest_row_cost:
					lowest_row_cost = cost_matrix[i][j]	# Finding the lowest cost for each row
			for j in range(len(cost_matrix[i])):
				if(cost_matrix[i][j] != np.inf):
					cost_matrix[i][j] -= lowest_row_cost	# Reducing the cost matrix and updating the lower bound
			lower_bound += lowest_row_cost
		for i in range(len(cost_matrix)):
			lowest_column_cost = np.inf
			for j in range(len(cost_matrix[i])):
				if cost_matrix[j][i] < lowest_column_cost:
					lowest_column_cost = cost_matrix[j][i]	# Finding the lowest column cost
			for j in range(len(cost_matrix[i])):	# Reducing the columns
				if(cost_matrix[j][i] != np.inf):
					cost_matrix[j][i] -= lowest_column_cost	# Updating the cost matrix and updating the lower bound
			lower_bound+= lowest_column_cost
		return lower_bound, cost_matrix

	def findLowerBoundReduceMatrix(self, lower_bound, cost_matrix, start_city_index, destination_city_index):
		lower_bound += cost_matrix[destination_city_index][start_city_index]
		cost_matrix[start_city_index][destination_city_index] = np.inf	# Source and destination nodes to inf in the cost matrix
		cost_matrix[destination_city_index][start_city_index] = np.inf
		for i in range(len(cost_matrix)):	# Updating the cost matrix rows and columns
			cost_matrix[start_city_index][i] = np.inf
			cost_matrix[i][destination_city_index] = np.inf
		for i in range(len(cost_matrix)):
			lowest_row_cost = np.inf
			for j in range(len(cost_matrix[i])):
				if cost_matrix[i][j] < lowest_row_cost:
					lowest_row_cost = cost_matrix[i][j]	# Finding the lowest cost for each row
			for j in range(len(cost_matrix[i])):
				if(cost_matrix[i][j] != np.inf):
					cost_matrix[i][j] -= lowest_row_cost # Updating the cost matrix
			lower_bound += lowest_row_cost
		for i in range(len(cost_matrix)):
			lowest_column_cost = np.inf
			for j in range(len(cost_matrix[i])):
				if cost_matrix[j][i] < lowest_column_cost and cost_matrix[j][i] != np.inf:
					lowest_column_cost = cost_matrix[j][i]	# Finding the lowest cost for each row
			for j in range(len(cost_matrix[i])):
				if(cost_matrix[j][i] != np.inf):
					cost_matrix[j][i] -= lowest_column_cost	# Updating the cost matrix
			if lowest_column_cost != np.inf:
				lower_bound += lowest_column_cost	# Updating the lower bound
		return lower_bound, cost_matrix
	'''
	def findLowerBoundReduceMatrix(self, lower_bound, cost_matrix, isFirstCalculation):
		lowest_row_cost_list = []
		for i in range(len(cost_matrix)):
			lowest_row_cost = np.inf
			for j in range(len(cost_matrix[i])):
				if cost_matrix[i][j] < lowest_row_cost:
					lowest_row_cost = cost_matrix[i][j]	# Finding the lowest cost for each row
			for j in range(len(cost_matrix[i])):
				if(cost_matrix[i][j] != np.inf):
					cost_matrix[i][j] -= lowest_row_cost
			lowest_row_cost_list.append(lowest_row_cost)
			if lowest_row_cost == np.inf and isFirstCalculation:
				lower_bound += lowest_row_cost	# Updating the lower bound
			elif lowest_row_cost == np.inf and isFirstCalculation == False:
				continue	# If the lowest row cost is infinity and it isnt the first time we are calculating the reduced matrix, we can just continue
			else:
				lower_bound += lowest_row_cost	# Updating the lower bound
		
		lowest_column_cost_list = []	# Reducing the columns
		for i in range(len(cost_matrix)):
			lowest_column_cost = np.inf
			for j in range(len(cost_matrix[i])):
				if cost_matrix[j][i] < lowest_column_cost:
					lowest_column_cost = cost_matrix[j][i]	# Reducing the rows
			#print(cost_matrix)
			if lowest_column_cost == np.inf and isFirstCalculation:
				lower_bound += lowest_column_cost	# Updating the lower bound
			elif lowest_column_cost == np.inf and isFirstCalculation == False:	# If the lowest column cost is infinity and 
				continue														# it isnt the first time we are calculating the reduced matrix, we can just continue
			else:
				lower_bound += lowest_column_cost
			for j in range(len(cost_matrix[i])):	# Reducing the columns
				if(cost_matrix[j][i] != np.inf):
					cost_matrix[j][i] -= lowest_column_cost
		return lower_bound, cost_matrix
	'''





	'''
	 <summary>
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


