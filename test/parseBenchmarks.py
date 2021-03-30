import re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Select your scene
#path of the .pbrt file, name of the scene
scenes = [
	'staircase',
	'bike',
]

# Select your integrator
# options: (name, heuristic, extra)
exec_filters = [
	'bdpt',
	'obdpt_direct',	
]

p_bdpt = re.compile('Rendering: \\[\\+\\+\\]  \\(\\d+\\.\\ds\\)')
p_direct_sampling = re.compile('Drawing Samples: \\[\\+\\+\\]  \\(\\d+\\.\\ds\\)')
p_direct_solving = re.compile('Solving: \\[\\+\\+\\]  \\(\\d+\\.\\ds\\)')


# Select your min and max lengths 
min_max= [
	(2, 3),
	(2, 4),
	(2, 5),
	(2, 6),
	(2, 7),
	(2, 8),
	(2, 9),
	(2, 10),
	(2, 15),
	(2, 22),
]

# Select your number of samples (or passes to be more accurate)
numbers_of_samples = 16
iterations = 16

max_lens = [mm[1] for mm in min_max]

p_float = re.compile('\\d+\\.\\d')


def toMatlab(values):
	res = [' '.join(map(str, value)) for value in values]
	res = '; '.join(res)
	res = '[' + res + ']'
	return res

def getTime(string):
	return float(p_float.findall(string)[0])

for scene in scenes:
	bdpt_res = [[0.0]*iterations for _ in range(len(min_max))]
	direct_sampling_res = [[0.0]*iterations for _ in range(len(min_max))]
	direct_solving_res = [[0.0]*iterations for _ in range(len(min_max))]

	for d_index in range(len(min_max)):
		mm = min_max[d_index]
		folder = 'benchmarks/' + scene + ('_s%d_L%d_l%d' % (numbers_of_samples, mm[1], mm[0]))
		for iter in range(iterations):
			
			bdpt_file = folder + '/bdpt_%d.txt' % iter
			with open(bdpt_file) as f:
				content = f.read()
				m = p_bdpt.findall(content)[-1]
				t = getTime(m)
				bdpt_res[d_index][iter] = t

			direct_file = folder + '/obdpt_direct_%d.txt' % iter
			with open(direct_file) as f:
				content = f.read()
				sampling_m = p_direct_sampling.findall(content)[-1]
				sampling_t = getTime(sampling_m)
				solving_m = p_direct_solving.findall(content)[-1]
				soling_t = getTime(solving_m)
				direct_sampling_res[d_index][iter] = sampling_t
				direct_solving_res[d_index][iter] = soling_t

	#print(scene + '_bdpt = ', toMatlab(bdpt_res), ';')
	#print(scene + '_direct_sampling = ', toMatlab(direct_sampling_res), ';')
	#print(scene + '_direct_solving = ', toMatlab(direct_solving_res), ';')

	bdpt = np.array(bdpt_res)
	direct_sampling = np.array(direct_sampling_res)
	direct_solving = np.array(direct_solving_res)

	

		