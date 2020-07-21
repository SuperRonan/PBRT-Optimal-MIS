import subprocess
import re
from colorama import Fore
from colorama import Style
import os

# generates the string with the selected integrator
def set_integrator(scene, integrator_str):
	start = '##INTEGRATOR-DEF-START'
	end = '##INTEGRATOR-DEF-END'
	replacement = integrator_str
	match = re.match(r'(.+%s\s*).+?(\s*%s.+)' % (start, end), scene, re.DOTALL)
	return match.group(1) + replacement + match.group(2)

def set_number_of_samples(scene, n):
	start = '##SAMPLER-DEF-START'
	end = '##SAMPLER-DEF-END'
	replacement = 'Sampler "random" "integer pixelsamples" %d' % (n)
	match = re.match(r'(.+%s\s*).+?(\s*%s.+)' % (start, end), scene, re.DOTALL)
	return match.group(1) + replacement + match.group(2)



pbrt_exe = '../build/Release/pbrt.exe'
result_folder = 'results/'
pbrt_scenes_folder = '../../pbrt-v3-scenes/'
num_threads = 16


# Select your scene
#path of the .pbrt file, name of the scene
scenes = [
	#['./scenes/simple-env.pbrt', 'simple-env'], # this scene contains an environment map
	#['./scenes/simple-area.pbrt', 'simple-area'],
	['./scenes/box.pbrt', 'box'], # a closed cornell box (diffuse)
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
	#[pbrt_scenes_folder + 'cornell-box/scene.pbrt', 'cornell'],
	#[pbrt_scenes_folder + 'water-caustic/scene.pbrt', 'water-caustic'],
	#[pbrt_scenes_folder + 'veach-mis/scene.pbrt', 'veach-mis'],
	#[pbrt_scenes_folder + 'veach-bidir/scene.pbrt', 'veach-bidir'],
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
]

# Select your integrator
# options: (name, option)
exec_filters = [
	('path', ''),                  # path tracer
	#("light", ''),             # light tracer (to reimplement)
	('bdpt', ''),	# bdpt (vanilla)
]

def filter_name(options):
    res = options[0]
    if(options[1] == 'no'):
        res = res + '-nolt'
    return res


# Select your min and max lengths 
min_max= [
	#(2, 2),
	#(2, 3),
	#(2, 4),
	#(2, 5, 5),
	#(2, 6),
	(2, 7, 7),
	#(2, 8),
	#(2, 9),
	#(2, 10),
	#(2, 11),
	#(2, 12),
	#(2, 13),
	#(2, 14),
	#(2, 15),
	#(2, 16),
	#(2, 17),
	#(2, 18),
	#(2, 19),
	#(2, 20),
	#(3, 3),
	#(4, 4),
	#(5, 5),
	#(6, 6),
	#(7, 7),
]

# Select your number of samples
numbers_of_samples = [
	#1,
	#2,
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
	#38768,
	#65536,
]

results = []
filenames = []

total = 0
passed = 0

for mm in min_max:

	print("depths: ", mm)

	min_len = mm[0]
	max_len = mm[1]
	max_depth = max_len-2
	min_depth = min_len-2
	max_opti_depth = mm[2]-2 if len(mm) == 3 else None

    # add additional integrator parameters etc here
	def integrator_string(integrator_name, min_depth, max_depth, LT, max_opti_depth):
		return 'Integrator "%s" "integer maxdepth" [ %d ] "integer mindepth" [ %d ]' % (integrator_name, max_depth, min_depth)

	def integrator_str(options, min_depth, max_depth, max_opti_depth):
		return integrator_string(options[0], min_depth, max_depth, options[1], max_opti_depth)

	for number_of_samples in numbers_of_samples:  
		print(number_of_samples, " samples")  
		for scene_info in scenes:
			print(scene_info[1])
			sub_folder = scene_info[1] + "_s%d_L%d_l%d/" % (number_of_samples, max_len, min_len)

			scene_path = scene_info[0]  
			with open(scene_path, 'r') as f:
				scene = f.read()

			for exec_filter in exec_filters:
				print(exec_filter[0])

				integrator = integrator_str(exec_filter, min_depth, max_depth, max_opti_depth)

				sc = set_integrator(scene, integrator)
				sc = set_number_of_samples(sc, number_of_samples)
				
				# We still keep this one for debugging purpose
				with open('scene.pbrt', 'w') as f:
					f.write(sc)
				# We need to do it like that for pbrt files that require to open other files in the same folder
				tmp_scene_path = scene_path[0: len(scene_path)-5] + 'Script.pbrt'

				with open(tmp_scene_path, 'w') as f:
					f.write(sc)
				name = filter_name(exec_filter)
				imgname = name + '.exr'
				filenames.append(imgname)
				command = [pbrt_exe, tmp_scene_path, '--outfile', result_folder + sub_folder + imgname, '--nthreads', str(num_threads)]
				if not os.path.exists(result_folder + sub_folder):
					os.makedirs(result_folder + sub_folder)

				res = subprocess.call(command)
				results.append([sub_folder+imgname, res])

				total = total + 1

				if(res == 0):
					print(('\n%s' + Fore.GREEN + ' returned %i' + Style.RESET_ALL) % (str(command), res))
					passed = passed + 1
				else:
					print(('\n%s' + Fore.YELLOW + ' returned %i' + Style.RESET_ALL) % (str(command), res))
				os.remove(tmp_scene_path)


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

# open all images, assumes tev is in the path (and we are on Windows / WSL)
# TODO find portable way to call tev here (cannot remove .exe because of WSL)
#      maybe catch error and call plain 'tev' instead
#viewer = ['tev.exe']
#viewer += filenames
#print(viewer)
#subprocess.call(viewer)