import subprocess
from colorama import Fore
from colorama import Style
import os
from helper import *
import sys


pbrt_exe = '../build/Release/pbrt.exe'
result_folder = 'results/'
pbrt_scenes_folder = '../../pbrt-v3-scenes/'
num_threads = 16


# Select your scene
#path of the .pbrt file, name of the scene
scenes = [
	#['./scenes/simple-env.pbrt', 'simple-env'], # this scene contains an environment map
	#['./scenes/simple-area.pbrt', 'simple-area'],
	#['./scenes/box.pbrt', 'box'], # a closed cornell box (diffuse)
	#['./scenes/box-glossy.pbrt', 'box-glossy'], # a closed cornell box with glossy materials on the cubes
	#['./scenes/box-phong.pbrt', 'box-phong'],
	#['./scenes/box-mirror.pbrt', 'box-mirror'], # a closed cornell box with mirrors on the cubes
	#['./scenes/box-caustic.pbrt', 'box-caustic'], # a closed cornel box with a glass sphere
	#['./scenes/box-caustic-empty.pbrt', 'box-caustic-empty'], # an open cornel box with a glass sphere
	#['./scenes/empty-box.pbrt', 'empty-box'],# an empty cornell box with the full ceiling as a light
	#['./scenes/test-mat/test-mat.pbrt', 'test-mat'],
	#['./scenes/veach-mis.pbrt', 'veach-mis'],
	#['./scenes/veach-mis-phong.pbrt', 'veach-mis-phong'],
	#['./scenes/box-three-lights.pbrt', 'box-three-lights'],
	#['./scenes/cornell-window.pbrt', 'cornell-window'],
	#['./scenes/box-sphere.pbrt', 'box-sphere'],
	#['./scenes/box-point.pbrt', 'box-point'],
	#['./scenes/box-invert.pbrt', 'box-invert'],
	#['./scenes/cornell-triangle.pbrt', 'cornell-triangle'],
	#['./scenes/cornell-large-light.pbrt', 'cornell-large'],
	#['./scenes/cornell-small-light.pbrt', 'cornell-small'],
	#[pbrt_scenes_folder + 'cornell-box/scene.pbrt', 'cornell'],
	#[pbrt_scenes_folder + 'water-caustic/scene.pbrt', 'water-caustic'],
	#[pbrt_scenes_folder + 'veach-mis/scene.pbrt', 'veach-mis'],
	#[pbrt_scenes_folder + 'veach-bidir/bidir.pbrt', 'veach-bidir'],
	#[pbrt_scenes_folder + 'caustic-glass/glass.pbrt', 'caustic-glass'],
	#[pbrt_scenes_folder + 'barcelona-pavilion/pavilion-night.pbrt', 'pavilion-night'],
	#[pbrt_scenes_folder + 'barcelona-pavilion/pavilion-day.pbrt', 'pavilion-day'],
	#[pbrt_scenes_folder + 'breakfast/breakfast.pbrt', 'breakfast'],
	#[pbrt_scenes_folder + 'pbrt-book/book.pbrt', 'pbrt-book'],
	#[pbrt_scenes_folder + 'sssdragon/dragon_50.pbrt', 'sssdragon'],
	#[pbrt_scenes_folder + 'staircase/scene.pbrt', 'staircase'],
	#[pbrt_scenes_folder + 'villa/villa-daylight.pbrt', 'villa-daylight'],
	#[pbrt_scenes_folder + 'villa/villa-lights-on.pbrt', 'villa-lights-on'],
	#[pbrt_scenes_folder + 'white-room/whiteroom-daytime.pbrt', 'whiteroom-daytime'],
	#[pbrt_scenes_folder + 'white-room/whiteroom-night.pbrt', 'whiteroom-night'],
	#[pbrt_scenes_folder + 'bunny-fur/f3-15.pbrt', 'bunny'],
	#[pbrt_scenes_folder + 'staircase2/scene.pbrt', 'staircase2'],
	#[pbrt_scenes_folder + 'staircase/scene.pbrt', 'staircase'],
	#[pbrt_scenes_folder + 'bathroom/bathroom.pbrt', 'bathroom'],
	#[pbrt_scenes_folder + 'contemporary-bathroom/contemporary-bathroom.pbrt', 'contemporary-bathroom'],
	[pbrt_scenes_folder + 'chopper-titan/chopper-titan.pbrt', 'bike'],
	#[pbrt_scenes_folder + 'chopper-titan/chopper-titan2.pbrt', 'bike2'],
	#[pbrt_scenes_folder + 'chopper-titan/chopper-titan3.pbrt', 'bike3'],
	#[pbrt_scenes_folder + 'sanmiguel/sanmiguel.pbrt', 'sanmiguel'],
	#[pbrt_scenes_folder + '2019/staircase1/scene/staircase1.pbrt', 'staircase_1_2019'],
	#[pbrt_scenes_folder + '2019/staircase2/scene/staircase2.pbrt', 'staircase_2_2019'],
	#[pbrt_scenes_folder + '2019/dining-room/scene/dining-room.pbrt', 'dining-room_2019'],
	#[pbrt_scenes_folder + '2019/veach/scene/veach.pbrt', 'veach_2019'], # Somehow corrupted and makes PBRT crash (precision errors or something)
	#['./scenes/glossy_env/scene.pbrt', 'glossy_env'],
]

# Select your integrator
# options: (name, heuristic, extra)
exec_filters = [
	#('path', ''),
	#("light", ''),			# light tracer (to reimplement)
	('bdpt', ''),

	('obdpt', 'balance'),
	('obdpt', 'power'),	
	('obdpt', 'cutoff'),	
	('obdpt', 'maximum'),
	#('obdpt', 'naive'),	
	('obdpt', 'direct', 'conservative'),	
	('obdpt', 'direct'),	

	#('opath', 'balance', [('BSDF', 1), ("Li", 1)]),
	#('opath', 'power', [('BSDF', 1), ("Li", 1)]),
	#('opath', 'direct', [('BSDF', 1), ("Li", 1)]),
	#('opath', 'direct', [('SS', 1), ('PP', 1)]),
	#('opath', 'balance', [('SS', 1), ('PP', 1)]),
	#('opath', 'direct', [('SS', 1), ('PP', 1), ('Li', 1)]),
	#('opath', 'balance', [('SS', 1), ('PP', 1), ('Li', 1)]),
	#('opath', 'balance', [('SS', 1)]),
	#('opath', 'balance', [('PP', 1)]),
	#('opath', 'balance', [('Li', 1)]),
	#('opath', 'balance', [('SP', 2)]),
	#('opath', 'balance', [('BSDF', 1)]),
	#('opath', 'progressive', [('SS', 1), ("SP", 1), ("Li", 1)]),
	#('opath', 'direct', [('SS', 1), ("SP", 1), ("Li", 1)]),
	#('opath', 'direct', [('SP', 1), ("PP", 1), ("Li", 1), ('BSDF', 1)]),
	#('opath', 'balance', [('SP', 1), ("PP", 1), ("Li", 1), ('BSDF', 1)]),
	#('opath', 'direct', [('SP', 1), ('PP', 1), ('BSDF', 1)]),
	#('opath', 'direct', [('SP', 1), ('PP', 1)]),
	#('opath', 'direct', [('SP', 1), ('PP', 1)], 'conservative'),
	#('opath', 'progressive', [('SP', 1), ('PP', 1)]),
	#('opath', 'progressive', [('SP', 1), ('PP', 1)], 'conservative'),
	#('opath', 'balance', [('SP', 1), ('PP', 1)]),
	#('opath', 'balance', [('SS', 1), ("SP", 1), ("Li", 1)]),
]


# Select your min and max lengths 
min_max= [
	#(2, 2),
	#(2, 3),
	#(2, 4),
	#(2, 5),
	#(2, 6),
	(2, 7),
	#(2, 8),
	#(2, 9),
	#(2, 10),
	#(2, 10, 5),
	#(2, 11),
	#(2, 12),
	#(2, 13),
	#(2, 14),
	#(2, 15),
	#(2, 16),
	#(2, 17),
	#(2, 18),
	#(2, 19),
	#(2, 22, 4),
	#(3, 3),
	#(4, 4),
	#(5, 5),
	#(6, 6),
	#(7, 7),
	#(5, 10),
	#(10, 10),
	#(22, 22)
]

# Select your number of samples (or passes to be more accurate)
numbers_of_samples = [
	#1,
	#2,
	#3,
	#4,
	#8, 
	16,
	#32, 
	#64, 
	#128, 
	#256, 
	#512, 
	#1024,
	#2048,
	#4096,
	#8192,
	#16384,
	#32768,
	#65536,
	#131072,
]

samplers = [
	'random',
	#'stratified',
	#'halton',
	#'02sequence',
	#'sobol',
	#'maxmindist',
	#'z',
	#'zhash',
	#'z-art',
	#'morton',
]

def main(args, i=None):

	total = 0
	passed = 0

	results = []
	filenames = []

	for mm in min_max:
		min_len = mm[0]
		max_len = mm[1]
		max_depth = max_len-2
		min_depth = min_len-2
		max_opti_depth = mm[2]-2 if len(mm) == 3 else None

		for number_of_samples in numbers_of_samples:  
			for scene_info in scenes:
				sub_folder = scene_info[1] + "_s%d_L%d_l%d/" % (number_of_samples, max_len, min_len)

				scene_path = scene_info[0]  
				
				pbrt_scene = PBRTSceneFile(scene_path)

				for exec_filter in exec_filters:
					for sampler in samplers:
						print('#' * 200)
						print("depths: ", mm)
						print(number_of_samples, " samples per pixel")  
						print("Scene name: ", scene_info[1])
						print("Integrator: ", exec_filter)
						print('Sampler: ', sampler)

						pbrt_scene.integrator = integrator_str(exec_filter, min_depth, max_depth, max_opti_depth)
						pbrt_scene.sampler = sampler_str(sampler, number_of_samples)
						
						pbrt_scene.makeTmp()

						name = filter_name(exec_filter, max_opti_depth, sampler)
						imgname = name + '.exr'
						filenames.append(imgname)

						if not os.path.exists(result_folder + sub_folder):
							os.makedirs(result_folder + sub_folder)				
						
						command = [pbrt_exe, pbrt_scene.tmp_filename, '--outfile', result_folder + sub_folder + imgname, '--nthreads', str(num_threads)]
						if i is not None:
							folder = 'benchmarks/' + sub_folder
							if not os.path.exists(folder):
								os.mkdir(folder)
							filename = name + '_' + str(i) + '.txt'
							f = open(folder + filename, "w")
							res = subprocess.call(command, stdout=f)
							f.close()
						else:
							res = subprocess.call(command)
						
						results.append([sub_folder+imgname, res])

						total = total + 1

						if(res == 0):
							print(('\n%s' + Fore.GREEN + ' returned %i' + Style.RESET_ALL) % (str(command), res))
							passed = passed + 1
						else:
							print(('\n%s' + Fore.YELLOW + ' returned %i' + Style.RESET_ALL) % (str(command), res))
				
				pbrt_scene.finish()


	for res in results:
		if res[1] == 0:
			color = Fore.BLUE
			color2 = Fore.GREEN
		else:
			color = Fore.RED
			color2 = Fore.YELLOW    
		print('-' + color + res[0] + ' returned ' + color2 + str(res[1]) + Style.RESET_ALL)

	if passed == total:
		color = Fore.GREEN
		color2 = Fore.BLUE
	else:
		color = Fore.YELLOW
		color2 = Fore.RED
	print('passed: ' + color + str(passed) + Style.RESET_ALL + ' / ' + color2 + str(total) + Style.RESET_ALL)

if __name__ == "__main__":
	for i in range(1):
		main(sys.argv)