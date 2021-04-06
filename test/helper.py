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
	
	



def sampler_str(s, n):
	return 'Sampler "%s" "integer pixelsamples" %d' % (s, n)

def integrator_string(integrator_name, min_depth, max_depth, heuristic, max_opti_depth, is_conservative):
	conservative_str = 'true' if is_conservative else 'false'
	res = 'Integrator "%s" "integer maxdepth" [ %d ] "integer mindepth" [ %d ] "string heuristic" "%s" "bool conservative" "%s"' % (integrator_name, max_depth, min_depth, heuristic, conservative_str)
	if max_opti_depth is not None:
		res += '"integer maxoptidepth" [ %d ]' % max_opti_depth
	return res


def techniques_string(techniques):
	res = ''
	for tech in techniques:
		tech_str = ' "integer %s" [ %d ]' % tech
		res += tech_str
	return res

def techniques_name(techniques):
	res = ''
	for tech in techniques:
		tech_str = '_%s_%d' % tech
		res += tech_str
	return res

def integrator_str(options, min_depth, max_depth, max_opti_depth):
	is_conservative = any([option == 'conservative' for option in options])

	res = integrator_string(options[0], min_depth, max_depth, options[1], max_opti_depth, is_conservative)
	if options[0] == 'opath' and len(options) >= 3:
		techniques = options[2]
		res += techniques_string(techniques)
	return res

def filter_name(options, max_opti_depth, sampler):
	res = options[0]
	if options[0] == 'obdpt' or options[0] == 'opath':
		res += '_' + options[1]
	if options[0] == 'opath' and len(options) >= 3:
		res += techniques_name(options[2]) 
	if max_opti_depth is not None:
		res += "_mo%d" % max_opti_depth
	is_conservative = any([option == 'conservative' for option in options])
	if is_conservative:
		res += '_conservative'
	if sampler != 'random':
		res += '_' + sampler
	return res

