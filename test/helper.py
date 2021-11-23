import re
import os


class PBRTSceneFile:

	filename = ""

	lines = []

	integrator = ""
	sampler = ""
	distrib = ""

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
	
	
def extract_light_strategy(options):
	light_strategy = None
	if 'uniform' in options:
		light_strategy = 'uniform'
	if 'power' in options:
		light_strategy = 'power'
	if 'spatial' in options:
		light_strategy = 'spatial'
	return light_strategy

def extract_estimation_strategy(options):
	res = None
	if 'strict' in options:
		res = 'strict'
	if 'loose' in options:
		res = 'loose'
	return res


def sampler_str(s, n):
	return 'Sampler "%s" "integer pixelsamples" %d' % (s, n)

def integrator_string(integrator_name, min_depth, max_depth, heuristic, max_opti_depth, estimation_strategy=None, light_strategy=None):
	res = 'Integrator "%s" "integer maxdepth" [ %d ] "integer mindepth" [ %d ] "string heuristic" "%s"' % (integrator_name, max_depth, min_depth, heuristic)
	if estimation_strategy is not None:
		strict_str = 'true' if estimation_strategy == 'strict' else 'false'
		res += ' "bool strict" "%s"' % strict_str
	if light_strategy is not None:
		res += '"string lightsamplestrategy" "%s"' % light_strategy
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

def extract_extra_options(options):
	return options[2:]

def integrator_str(options, min_depth, max_depth, max_opti_depth):
	extra_options = extract_extra_options(options)
	estimation_strat = extract_estimation_strategy(extra_options)
	light_strategy = extract_light_strategy(extra_options)
	res = integrator_string(options[0], min_depth, max_depth, options[1], max_opti_depth, estimation_strat, light_strategy)
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
	extra_options = extract_extra_options(options)
	light_strategy = extract_light_strategy(extra_options)
	if light_strategy is not None:
		res += '_' + light_strategy
	estimation_strat = extract_estimation_strategy(extra_options)
	if estimation_strat is not None:
		res += '_' + estimation_strat
	if sampler != 'random':
		res += '_' + sampler
	return res

