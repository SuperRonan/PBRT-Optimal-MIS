import re
import os


class PBRTSceneFile:

	filename = ""

	lines = []

	integrator = ""
	sampler = ""

	tmp_filename = None

	def __init__(self, filename):
		self.filename = filename
		file = open(filename, 'r')
		tmp_lines = file.read().splitlines()
		file.close()
		self.lines = []
		for line in tmp_lines:
			if len(line) > 1:
				first_word = line.partition(' ')[0]
				if first_word == 'Integrator':
					self.integrator = line
				elif first_word == 'Sampler':
					self.sampler = line
				else:
					self.lines.append(line)
	
	def makeTmp(self):
		self.tmp_filename = self.filename[0: len(self.filename)-5] + '_TMP.pbrt'

		new_list = self.lines.copy()
		new_list.insert(0, self.integrator)
		new_list.insert(0, self.sampler)

		content = '\n'.join(new_list)

		file = open(self.tmp_filename, 'w')
		file.write(content)
		file.close()

	def finish(self):
		if self.tmp_filename is not None:
			os.remove(self.tmp_filename)
			self.tmp_filename = None
	
	



def sampler_str(n):
	return 'Sampler "random" "integer pixelsamples" %d' % (n)

def integrator_string(integrator_name, min_depth, max_depth, LT, max_opti_depth):
	return 'Integrator "%s" "integer maxdepth" [ %d ] "integer mindepth" [ %d ]' % (integrator_name, max_depth, min_depth)

def integrator_str(options, min_depth, max_depth, max_opti_depth):
	return integrator_string(options[0], min_depth, max_depth, options[1], max_opti_depth)

def set_integrator(scene, integrator_str):
	start = '##INTEGRATOR-DEF-START'
	end = '##INTEGRATOR-DEF-END'
	replacement = integrator_str
	match = re.match(r'(.+%s\s*).+?(\s*%s.+)' % (start, end), scene, re.DOTALL)
	return match.group(1) + replacement + match.group(2)

def set_number_of_samples(scene, n):
	start = '##SAMPLER-DEF-START'
	end = '##SAMPLER-DEF-END'
	replacement = sampler_str(n)
	match = re.match(r'(.+%s\s*).+?(\s*%s.+)' % (start, end), scene, re.DOTALL)
	return match.group(1) + replacement + match.group(2)

def filter_name(options):
	res = options[0]
	if(options[1] is not None):
		res = res + options[1]
	return res

