import os,multiprocessing


##########################################

#
class Consumer(multiprocessing.Process):
	
	''' Simply a process Class'''
	
	def __init__(self,tasks_queue=multiprocessing.JoinableQueue(),results=None):
		multiprocessing.Process.__init__(self)
		self.task_queue=tasksk_queue # get things to do
		self.results_queue=results # put the results somewhere
		
	def run(self):
		process_name=self.name
		while True:
			next_task = self.task_queue.get()
			if next_task == None:
				# Poison pills: die Consumer, die!
				self.task_queue.task_done() # my last task is done, sir...
				break
			print process_name, 'working hard on', next_task.name
			answer = next_task()
			self.task_queue.task_done()
			
			if self.results_queue != None:
				# store result in results queue
				self.result_queue.put(answer)
		return


#
class Task(dict):

	''' Define a generic task.
		A target function (or callable) must be specified with
		target=callable
	'''

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self.__dict__ = self
	
	def __call__(self,*args):
		return self.target(*args)
		
		
#
class Multiprocess(object):
	
	''' Object for perform multiprocessing '''
	
	def __init__(self,n_cpus=1):
		self.tasks=multiprocessing.JoinableQueue()
		self.results=multiprocessing.Queue()
		self.n_cpus=n_cpus
		self.consumers=None
		
	def make_factory(self):
		self.consumers= [Consumer(self.tasks,self.results)
							for cpu in xrange(n_cpus)]
							
	def get_jobs(self,jobs):
		# check if you have a factory
		if self.consumers == None:
			self.make_factory()
		# get orders
		for job in jobs:
			self.task.put(job)
		# put poison pills in the brew
		for i in self.consumers:
			self.task.put(None)

	def start_production(self):
		# check if you have a factory
		if self.consumers == None:
			self.make_factory()
		# check if you have jobs to do
		if self.tasks.empty:
			print 'no jobs!'
			return 
		# do the jobs, slackers!
		for c in self.consumers:
			c.start()
		# wait the jobs to be completed...
		self.tasks.join()
		
	def get_results(self):
		# its a generator
		# if you have results
		if self.results.full():
			yield self.results.get()
		else: return
		
		




##########################################
