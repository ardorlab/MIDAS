import os
import sys
import yaml
import numpy
from matplotlib import pyplot
import matplotlib as mpl
import random
from midas.utils.solution_types import return_triangular_string
'''
Addtional import for pareto front graph
'''
import glob
from PIL import Image


class Optimization_Metric_Toolbox(object):
    """
    Class for measuring the performance of the optimization.
    """
    def __init__(self):
        self.best_solution_over_time = {}
        self.best_solution_over_time['fitness'] = []
        self.best_value_over_time = {}
        self.worst_value_over_time = {}
        self.worst_solution_over_time = {}
        self.worst_solution_over_time['fitness'] = []
        self.worst_fitness_over_time = []
        self.average_fitness_over_time = []
        self.average_value_over_time = {}
        self.best_fitness_over_time = []
        self.all_param_over_time={}

    def check_best_worst_average(self,population):
        """
        Calculates the best worst and average of population in the generation.
        """
        max_fitness = -10000000000
        min_fitness = 1000000000
        #Determining the best worst and average values for each parameter. 
        for param in population[0].parameters:
            if param in self.best_value_over_time:
                pass
            else:
                self.best_value_over_time[param] = []
            if param in self.worst_value_over_time:
                pass
            else:
                self.worst_value_over_time[param] = []
            if param in self.average_value_over_time:
                pass
            else:
                self.average_value_over_time[param] = []
            param_list = []
            for solution in population:
                param_list.append(solution.parameters[param]['value'])
            min_ = numpy.min(param_list)
            max_ = numpy.max(param_list)
            ave_ = numpy.average(param_list)
            obj_ = population[0].parameters[param]['goal']
            if obj_ == 'maximize':                               #Might need to include all additional possible objectives here
                self.best_value_over_time[param].append(max_)
                self.worst_value_over_time[param].append(min_)
            else:
                self.best_value_over_time[param].append(min_)
                self.worst_value_over_time[param].append(max_)
            self.average_value_over_time[param].append(ave_)

        list_ = []
        for solution in population:
            list_.append(solution.fitness)
            if solution.fitness > max_fitness:
                max_fitness = solution.fitness
                self.current_best = solution
            if solution.fitness < min_fitness:
                min_fitness = solution.fitness
                self.current_worst = solution

        self.average_fitness_over_time.append(numpy.average(list_))

        foo = self.current_best.fitness
        self.best_solution_over_time['fitness'].append(foo)
        for param in self.current_best.parameters:
            if param in self.best_solution_over_time:
                pass
            else:
                self.best_solution_over_time[param] = []
            foo = self.current_best.parameters[param]["value"]
            self.best_solution_over_time[param].append(foo)

        foo = self.current_worst.fitness
        self.worst_solution_over_time['fitness'].append(foo)
        for param in self.current_worst.parameters:
            if param in self.worst_solution_over_time:
                pass
            else:
                self.worst_solution_over_time[param] = []
            foo = self.current_worst.parameters[param]["value"]
            self.worst_solution_over_time[param].append(foo)
    
    def record_all_param(self,population_class,generation_class,flag):
        self.all_param_over_time[generation_class.current] = {}
        for sol in population_class.parents:
            param_list = [param for param in sol.parameters]
            for param in sol.parameters:
                self.all_param_over_time[generation_class.current][param]=[]
            break 
        # appending values
        if flag:
            for sol in population_class.parents:
                if sol.name.find('initial') < 0:
                    for param in sol.parameters:
                        tmp = sol.parameters[param]['value']
                        self.all_param_over_time[generation_class.current][param].append(tmp)
        for sol in population_class.children:
            if sol.name.find('initial') < 0:
                for param in sol.parameters:
                    tmp = sol.parameters[param]['value']
                    self.all_param_over_time[generation_class.current][param].append(tmp)
        return

    def make_gif(self, frame_folder, inpstr):
        frames = [Image.open(image) for image in sorted(glob.glob(f"{frame_folder}/*.jpg")\
            ,key=lambda x: self.get_int(x))]
        # extract first image from iterator
        frame_one = frames[0]
        frame_one.save("pareto_front_"+inpstr+".gif", format="GIF", \
                        append_images=frames, save_all=True, duration=500, loop=0)

    def get_int(self, str):
        tmp = str.split('_')
        tmp = tmp[-1].split('.')
        num = int(tmp[0])
        return num

    def gen_plot(self, target, remain_param, c):
        ''' Generate corvar plot 
        target      : target parameter for covar plot
        remain_param: remain list of param 
        c           : color list           
        '''

        for param in remain_param:
            os.makedirs('image_dir_'+target+'_'+param, exist_ok=True)
            path = os.getcwd() + "/image_dir_"+target+"_"+param
            for gen in range(len(self.all_param_over_time)): # loop generation 
                fig,ax = pyplot.subplots(2,1, figsize=(5,5),  gridspec_kw={'height_ratios': [25,1]},  constrained_layout = True)
                for t in range (gen+1):
                    sc= ax[0].scatter(self.all_param_over_time[t][param],self.all_param_over_time[t][target] , s = 50, color=c[t])
                ax[0].set_xlabel(str(param))
                ax[0].set_ylabel(str(target))
                cmap = mpl.colors.ListedColormap(c)
                cmap.set_over('0.25')
                cmap.set_under('0.75')
                bounds = range(1,len(self.all_param_over_time)+2)
                norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
                cb2 = mpl.colorbar.ColorbarBase(ax[1], cmap=cmap,
                                                norm=norm,
                                                #boundaries=[0] + bounds + [13],
                                                extend='both',
                                                ticks=[bounds[0],max(bounds)],
                                                spacing='proportional',
                                                orientation='horizontal')
                cb2.set_label('Generation')
                name = path +'/image_'+str(gen+1)+'.jpg'
                pyplot.savefig(name)
                pyplot.close()
            self.make_gif(path+'/', param+'_'+target )

    def plotter(self):
        """
        Plots the best, worst and average fitness and parameter values for
        the optimization.
        Edit log:
        - - - -- -: Brian created
        26/06/2022: Khang - added pareto front graph (gif)
        """
        for param in self.best_solution_over_time:
            pyplot.figure()
            foo1 = range(len(self.best_solution_over_time[param]))
            foo2 = range(len(self.worst_solution_over_time[param]))
            pyplot.plot(foo1,self.best_solution_over_time[param],label='Best')
            pyplot.plot(foo2,self.worst_solution_over_time[param],label='Worst')
            if param == 'fitness':
                foo3 = range(len(self.average_fitness_over_time))
                pyplot.plot(foo3,self.average_fitness_over_time,label='Average')
            else:
                foo3 = range(len(self.average_value_over_time[param]))
                pyplot.plot(foo3,self.average_value_over_time[param],label='Average')

            pyplot.xlabel('Generation')
            pyplot.ylabel(param)
            pyplot.title("{} for Best, Worst, and Average Solution over Time".format(param))
            pyplot.legend()
            pyplot.savefig("{}_over_time".format(param)) 
        list_param = [param for param in self.best_solution_over_time]
        list_param.pop(list_param.index('fitness')) # no fitness pareto front graph
        while list_param:
            get_colors = lambda n: list(map(lambda i: "#" + "%06x" % random.randint(0, 0xFFFFFF),range(n)))
            color = get_colors(len(self.all_param_over_time)) # sample return:  ['#8af5da', '#fbc08c', '#b741d0', '#e599f1', '#bbcb59', '#a2a6c0']
            target = list_param.pop(0)
            self.gen_plot(target, list_param, color)

    def write_track_file(self,population_class,generation_class):
        """
        Writes a file making it easier to track what is occuring in the optimization.
        """
        track_file = open('optimization_track_file.txt','a')
        line_ = "Current generation of the optimization"
        line_ +=" {} \n".format(generation_class.current)
        track_file.write(line_)
        line_ = "Size of Parent Population"
        line_ += " {} \n".format(len(population_class.parents))
        track_file.write(line_)
        line_ ="Size of Child Population"
        line_ +=" {} \n".format(len(population_class.children))
        track_file.write(line_)
        if population_class.solution_front:
            line_ = "Size of Solution Front "
            line_ += "{} \n".format(len(population_class.solution_front))
            track_file.write(line_)
        #Record information on the current best solution in the optimization
        line_ = "Best Solution {}\n".format(self.current_best.name)
        track_file.write(line_)
        line_ = "Best Fitness {}\n".format(self.current_best.fitness)
        track_file.write(line_)
        for param in self.current_best.parameters:
            value = self.current_best.parameters[param]['value']
            line_ = "{} for best solution: {}  \n".format(param,value)
            track_file.write(line_)
        track_file.write("Best Solution Genome\n")
        gene_count = 0
        if type(self.current_best.genome) == list:
            for gene in self.current_best.genome:
                track_file.write("{}  ".format(gene))
                gene_count += 1
                if gene_count == 10:
                    gene_count = 0
                    track_file.write("\n")
        track_file.write("\n")

        #Record information on the current worst solution in the optimization
        line_ = "Worst Solution {}\n".format(self.current_worst.name)
        track_file.write(line_)
        line_ = "Worst Fitness {}\n".format(self.current_worst.fitness)
        track_file.write(line_)
        for param in self.current_worst.parameters:
            value = self.current_worst.parameters[param]['value']
            line_ = "{} for worst solution: {}  \n".format(param,value)
            track_file.write(line_)
        track_file.write("Worst Solution Genome\n")
        gene_count = 0
        for gene in self.current_worst.genome:
            track_file.write("{}  ".format(gene))
            gene_count += 1
            if gene_count == 10:
                gene_count = 0
                track_file.write("\n")
        track_file.write('\n')
        for param in self.best_value_over_time:
            foo = self.best_value_over_time[param][-1]
            line = "{}: Best = {}\n".format(param,foo)
            track_file.write(line)
            foo = self.average_value_over_time[param][-1]
            line = "{}: Average = {}\n".format(param,foo)
            track_file.write(line)
            foo = self.worst_value_over_time[param][-1]
            line = "{}: Worst = {}\n".format(param,foo)
            track_file.write(line)
        track_file.write("\n\n\n")


        track_file.close()

    def record_optimized_solutions(self,population_class):
        """
        Records the optimized solutions of the optimization for review or use as at a later time.
        Optimized solution results are recorded as a dictionary and stored in a yaml file.
        """

        if population_class.solution_front:
            list_ = population_class.solution_front[:]
        else:
            list_ = population_class.parents[:]

        optimized_dictionary = {}
        for solution in list_:
            optimized_dictionary[solution.name] = {}
            optimized_dictionary[solution.name]['genome'] = solution.genome
            optimized_dictionary[solution.name]['fitness'] = float(solution.fitness)
            optimized_dictionary[solution.name]['parameters'] = solution.parameters

        with open("optimized_solutions.yaml",'w') as yaml_file:
            yaml.dump(optimized_dictionary,yaml_file)

class Plotter(object):
    """
    Class for plotting solutions to optimization problems
    """
    colors = ["xkcd:blue","xkcd:green","xkcd:orange","xkcd:red","xkcd:black","xkcd:purple","xkcd:yellow","xkcd:cobalt",
              "xkcd:brown","xkcd:light blue","xkcd:teal","xkcd:light green","xkcd:magenta","xkcd:sky blue","xkcd:pink",
              "xkcd:turquoise","xkcd:lavender","xkcd:bright green","xkcd:brick red","xkcd:sea green","xkcd:mint green",
              "xkcd:royal blue","xkcd:navy blue","xkcd:lilac","xkcd:hot pink","xkcd:lime","xkcd:olive drab","xkcd:tan",
              "xkcd:indigo","xkcd:light brown","xkcd:aquamarine","xkcd:navy","xkcd:grass green","xkcd:greenish yellow",
              "xkcd:blue green","xkcd:seafoam green","xkcd:kelly green","xkcd:puke green","xkcd:pea green","xkcd:puce",
              "xkcd:pastel green","xkcd:dark orange","xkcd:robin's egg blue","xkcd:white","xkcd:leaf green","xkcd:sky",
              "xkcd:dark magenta","xkcd:red orange","xkcd:yellowish green","xkcd:purplish pink","xkcd:cornflower blue",
              "xkcd:light violet","xkcd:dusty rose","xkcd:medium blue","xkcd:tangerine","xkcd:pinkish red","xkcd:clay",
              "xkcd:purplish blue","xkcd:greyish blue","xkcd:grape","xkcd:light olive","xkcd:dark salmon","xkcd:vomit",
              "xkcd:bright red","xkcd:azure","xkcd:blue purple","xkcd:dark turquoise","xkcd:electric blue","xkcd:snot",
              "xkcd:purple blue","xkcd:aqua blue","xkcd:dark mauve","xkcd:greenish brown","xkcd:red brown","xkcd:ocre",
              "xkcd:dark tan","xkcd:green blue","xkcd:bluish green","xkcd:pastel blue","xkcd:cobalt blue","xkcd:adobe",
              "xkcd:blood red","xkcd:sage green","xkcd:terracotta","xkcd:pastel purple","xkcd:pumpkin","xkcd:charcoal",
              "xkcd:greyish brown","xkcd:soft blue","xkcd:easter green","xkcd:bluey purple","xkcd:marine","xkcd:wheat",
              "xkcd:pine green","xkcd:deep red","xkcd:dark olive","xkcd:dirty yellow","xkcd:orchid","xkcd:grey purple",
              "xkcd:aqua green","xkcd:raspberry","xkcd:greyish purple","xkcd:rose pink","xkcd:neon pink","xkcd:cerise",
              "xkcd:orange red","xkcd:dark rose","xkcd:brownish red","xkcd:pink purple","xkcd:pinky purple","xkcd:mud",
              "xkcd:burnt umber","xkcd:dull blue","xkcd:pale brown","xkcd:umber","xkcd:light sky blue","xkcd:ice blue",
              "xkcd:eggshell","xkcd:jungle green","xkcd:dark peach","xkcd:green brown","xkcd:bright teal","xkcd:ocean",
              "xkcd:bordeaux","xkcd:light blue green","xkcd:yellowish","xkcd:snot green","xkcd:red violet","xkcd:pine",
              "xkcd:teal blue","xkcd:denim blue","xkcd:dark lime green","xkcd:dull yellow","xkcd:pistachio","xkcd:tea",
              "xkcd:dark mustard","xkcd:pea soup green","xkcd:bubblegum pink","xkcd:barbie pink","xkcd:military green",
              "xkcd:amber","xkcd:mid blue","xkcd:shit brown","xkcd:hospital green","xkcd:purpleish blue","xkcd:auburn",
              "xkcd:sickly green","xkcd:melon","xkcd:dusky rose","xkcd:brown orange","xkcd:pastel yellow","xkcd:maize",
              "xkcd:primary blue","xkcd:orangey red","xkcd:pale lilac","xkcd:rust red","xkcd:easter purple","xkcd:sea",
              "xkcd:ocean green","xkcd:mustard green","xkcd:poop brown","xkcd:olive brown","xkcd:pink red","xkcd:pear",
              "xkcd:taupe","xkcd:lime green","xkcd:light purple","xkcd:violet","xkcd:dark green","xkcd:mustard yellow",
              "xkcd:dark brown","xkcd:deep purple","xkcd:chartreuse","xkcd:bright pink","xkcd:dark cyan","xkcd:celery",
              "xkcd:crimson","xkcd:fuchsia","xkcd:blue grey","xkcd:slate blue","xkcd:light olive green","xkcd:emerald",
              "xkcd:reddish pink","xkcd:yellow orange","xkcd:light turquoise","xkcd:moss","xkcd:pale red","xkcd:denim",
              "xkcd:grey brown","xkcd:off white","xkcd:pale orange","xkcd:greenish blue","xkcd:baby pink","xkcd:ivory",
              "xkcd:reddish purple","xkcd:medium green","xkcd:deep pink","xkcd:very light blue","xkcd:fire engine red",
              "xkcd:apple green","xkcd:emerald green","xkcd:powder blue","xkcd:light cyan","xkcd:lemon","xkcd:blurple",
              "xkcd:lighter green","xkcd:wine","xkcd:grey","xkcd:dull green","xkcd:deep sky blue","xkcd:dark sky blue",
              "xkcd:dark olive green","xkcd:greeny yellow","xkcd:faded purple","xkcd:green yellow","xkcd:windows blue",
              "xkcd:macaroni and cheese","xkcd:darker blue","xkcd:marine blue","xkcd:strong blue","xkcd:brownish grey",
              "xkcd:light navy blue","xkcd:very dark green","xkcd:very dark blue","xkcd:yellow ochre","xkcd:aubergine",
              "xkcd:toxic green","xkcd:faded green","xkcd:purple pink","xkcd:camo green","xkcd:dusty pink","xkcd:shit",                         
              "xkcd:light grey blue","xkcd:dull pink","xkcd:cadet blue","xkcd:pumpkin orange","xkcd:very light purple",
              "xkcd:pinkish tan","xkcd:pale olive","xkcd:minty green","xkcd:dusty green","xkcd:teal green","xkcd:dusk",
              "xkcd:purple red","xkcd:slate grey","xkcd:light tan","xkcd:evergreen","xkcd:blueberry","xkcd:watermelon",
              "xkcd:bright lavender","xkcd:brownish orange","xkcd:pinkish purple","xkcd:robin's egg","xkcd:cornflower",
              "xkcd:greyish green","xkcd:bright purple","xkcd:forest green","xkcd:aqua","xkcd:dark purple","xkcd:cyan",
              "xkcd:light lavender","xkcd:electric green","xkcd:pale lavender","xkcd:bright yellow","xkcd:bright aqua",
              "xkcd:yellow green","xkcd:blue violet","xkcd:asparagus","xkcd:pale grey","xkcd:dark blue","xkcd:apricot",
              "xkcd:pale blue","xkcd:chocolate","xkcd:baby blue","xkcd:bright orange","xkcd:scarlet","xkcd:acid green",
              "xkcd:golden yellow","xkcd:purple grey","xkcd:bluegreen","xkcd:dark red","xkcd:salmon","xkcd:dusky pink",
              "xkcd:beige","xkcd:slate","xkcd:mauve","xkcd:olive","xkcd:steel","xkcd:light lime green","xkcd:marigold",
              "xkcd:maroon","xkcd:muddy green","xkcd:dull orange","xkcd:light mauve","xkcd:mud brown","xkcd:sandstone",
              "xkcd:buff","xkcd:turquoise green","xkcd:muddy brown","xkcd:tomato","xkcd:carnation pink","xkcd:tealish",
              "xkcd:british racing green","xkcd:light sea green","xkcd:orangish brown","xkcd:pale pink","xkcd:seafoam",
              "xkcd:light forest green","xkcd:orangish red","xkcd:light yellow","xkcd:royal purple","xkcd:neon purple",
              "xkcd:spring green","xkcd:tree green","xkcd:shamrock green","xkcd:yellowish brown","xkcd:robin egg blue",
              "xkcd:fluorescent green","xkcd:electric purple","xkcd:medium purple","xkcd:light orange","xkcd:burgundy",
              "xkcd:violet blue","xkcd:faded pink","xkcd:light aqua","xkcd:dark yellow","xkcd:khaki green","xkcd:mint",
              "xkcd:army green","xkcd:dark grey","xkcd:grey blue","xkcd:burnt red","xkcd:greenish","xkcd:vomit yellow",
              "xkcd:lemon yellow","xkcd:purplish brown","xkcd:reddish brown","xkcd:canary yellow","xkcd:darkish green",
              "xkcd:dark teal","xkcd:vermillion","xkcd:soft green","xkcd:cranberry","xkcd:purpleish","xkcd:ocean blue",
              "xkcd:dark periwinkle","xkcd:bluish purple","xkcd:bluish grey","xkcd:brownish green","xkcd:orange brown",
              "xkcd:light magenta","xkcd:terra cotta","xkcd:pale turquoise","xkcd:dust","xkcd:salmon pink","xkcd:ecru",
              "xkcd:light blue grey","xkcd:leaf","xkcd:orangish","xkcd:pale olive green","xkcd:orangered","xkcd:pinky",
              "xkcd:kiwi green","xkcd:boring green","xkcd:light pastel green","xkcd:light seafoam green","xkcd:squash",
              "xkcd:electric pink","xkcd:purplish grey","xkcd:banana","xkcd:neon yellow","xkcd:greyish","xkcd:pinkish",
              "xkcd:dark mint","xkcd:light urple","xkcd:midnight purple","xkcd:pinkish orange","xkcd:yellowish orange",
              "xkcd:mud green","xkcd:dull teal","xkcd:deep lavender","xkcd:vivid blue","xkcd:fresh green","xkcd:algae",
              "xkcd:muted purple","xkcd:greyish pink","xkcd:burnt yellow","xkcd:olive green","xkcd:plum","xkcd:barney",
              "xkcd:barney purple","xkcd:burnt sienna","xkcd:puke yellow","xkcd:dark lilac","xkcd:reddish","xkcd:sand",
              "xkcd:mid green","xkcd:periwinkle","xkcd:light navy","xkcd:pale green","xkcd:ultramarine","xkcd:mustard",
              "xkcd:midnight blue","xkcd:dark pink","xkcd:sea blue","xkcd:light grey","xkcd:bright blue","xkcd:forest",
              "xkcd:moss green","xkcd:light pink","xkcd:eggplant","xkcd:purplish","xkcd:burnt orange","xkcd:goldenrod",
              "xkcd:neon green","xkcd:pea soup","xkcd:khaki","xkcd:blush","xkcd:cream","xkcd:coral","xkcd:brown green",
              "xkcd:hunter green","xkcd:steel blue","xkcd:light teal","xkcd:grey green","xkcd:deep blue","xkcd:silver",
              "xkcd:bright turquoise","xkcd:dark aqua","xkcd:avocado","xkcd:light indigo","xkcd:peach","xkcd:red pink",
              "xkcd:purpley blue","xkcd:swamp green","xkcd:dirt","xkcd:puke","xkcd:rose","xkcd:gold","xkcd:red purple",
              "xkcd:lightish blue","xkcd:golden brown","xkcd:cerulean blue","xkcd:shocking pink","xkcd:lightish green",
              "xkcd:chocolate brown","xkcd:carolina blue","xkcd:dusty purple","xkcd:yellow brown","xkcd:forrest green",
              "xkcd:light lilac","xkcd:mulberry","xkcd:rust","xkcd:sage","xkcd:brick","xkcd:wine red","xkcd:pale teal",
              "xkcd:dark forest green","xkcd:brownish yellow","xkcd:dark lavender","xkcd:dull purple","xkcd:pinky red",
              "xkcd:dark violet","xkcd:dark maroon","xkcd:drab green","xkcd:faded blue","xkcd:hot purple","xkcd:ochre",
              "xkcd:dark blue green","xkcd:yellowy green","xkcd:darkish blue","xkcd:sienna","xkcd:bronze","xkcd:grass",
              "xkcd:dull red","xkcd:bluish","xkcd:dark gold","xkcd:dirty pink","xkcd:slate green","xkcd:greenish grey",
              "xkcd:faded yellow","xkcd:bile","xkcd:viridian","xkcd:very light pink","xkcd:puke brown","xkcd:diarrhea",
              "xkcd:very light green","xkcd:cerulean","xkcd:light red","xkcd:ugly green","xkcd:dark beige","xkcd:jade",
              "xkcd:pale yellow","xkcd:pale purple","xkcd:brownish","xkcd:pastel pink","xkcd:dirty green","xkcd:mocha",
              "xkcd:periwinkle blue","xkcd:darker green","xkcd:orangey brown","xkcd:dark sea green","xkcd:pale violet",
              "xkcd:dusty blue","xkcd:neon blue","xkcd:true blue","xkcd:green grey","xkcd:prussian blue","xkcd:russet",
              "xkcd:aqua marine","xkcd:bright violet","xkcd:lighter purple","xkcd:reddish orange","xkcd:darker purple",
              "xkcd:bright magenta","xkcd:light maroon","xkcd:vomit green","xkcd:orange yellow","xkcd:bright sky blue",
              "xkcd:bright light blue","xkcd:bright lime green","xkcd:very dark purple","xkcd:poop","xkcd:murky green",
              "xkcd:brownish purple","xkcd:dark navy blue","xkcd:pinkish brown","xkcd:medium brown","xkcd:muted green",
              "xkcd:bottle green","xkcd:dirty orange","xkcd:light yellow green","xkcd:pinkish grey","xkcd:greeny blue",
              "xkcd:butter yellow","xkcd:dark khaki","xkcd:cherry red","xkcd:steel grey","xkcd:deep green","xkcd:iris",
              "xkcd:coffee","xkcd:light salmon","xkcd:kermit green","xkcd:irish green","xkcd:light lime","xkcd:copper",
              "xkcd:milk chocolate","xkcd:strawberry","xkcd:dirt brown","xkcd:jade green","xkcd:soft pink","xkcd:drab",
              "xkcd:sepia","xkcd:light rose","xkcd:light plum","xkcd:pink/purple","xkcd:medium grey","xkcd:dusky blue",
              "xkcd:purplish red","xkcd:slime green","xkcd:soft purple","xkcd:deep orange","xkcd:mahogany","xkcd:ruby",
              "xkcd:purpley pink","xkcd:bluey green","xkcd:dark indigo","xkcd:dark lime","xkcd:faded red","xkcd:camel",
              "xkcd:poo brown","xkcd:brown red","xkcd:dark navy","xkcd:baby poop","xkcd:grey pink","xkcd:dusky purple",
              "xkcd:midnight","xkcd:celadon","xkcd:heather","xkcd:ugly yellow","xkcd:sandy brown","xkcd:greenish teal",
              "xkcd:deep magenta","xkcd:greenish tan","xkcd:blue/purple","xkcd:deep violet","xkcd:rouge","xkcd:cherry",
              "xkcd:french blue","xkcd:light beige","xkcd:purply pink","xkcd:light peach","xkcd:warm grey","xkcd:fawn",
              "xkcd:lightblue","xkcd:off green","xkcd:gunmetal","xkcd:berry","xkcd:blood","xkcd:golden","xkcd:caramel",
              "xkcd:dark mint green","xkcd:rusty orange","xkcd:grassy green","xkcd:rose red","xkcd:earth","xkcd:ocher",
              "xkcd:ultramarine blue","xkcd:navy green","xkcd:seaweed","xkcd:kiwi","xkcd:fluro green","xkcd:pale aqua",
              "xkcd:mossy green","xkcd:darker pink","xkcd:sick green","xkcd:stone","xkcd:blue/green","xkcd:muted blue",
              "xkcd:turquoise blue","xkcd:light burgundy","xkcd:pale rose","xkcd:chestnut","xkcd:amethyst","xkcd:pale",
              "xkcd:sandy","xkcd:ugly pink","xkcd:dark plum","xkcd:perrywinkle","xkcd:pastel orange","xkcd:frog green",
              "xkcd:pale lime","xkcd:bright light green","xkcd:pale lime green","xkcd:peacock blue","xkcd:vivid green",
              "xkcd:baby poop green","xkcd:purpleish pink","xkcd:mustard brown","xkcd:leafy green","xkcd:almost black",
              "xkcd:bubblegum","xkcd:shamrock","xkcd:mango","xkcd:lime yellow","xkcd:purply blue","xkcd:avocado green",
              "xkcd:dull brown","xkcd:cool blue","xkcd:yellow/green","xkcd:heliotrope","xkcd:green apple","xkcd:apple",
              "xkcd:very pale green","xkcd:night blue","xkcd:lightgreen","xkcd:tomato red","xkcd:merlot","xkcd:velvet",
              "xkcd:light aquamarine","xkcd:sunshine yellow","xkcd:baby shit green","xkcd:vibrant purple","xkcd:bland",
              "xkcd:vibrant green","xkcd:bright lilac","xkcd:grape purple","xkcd:faded orange","xkcd:light periwinkle",
              "xkcd:dark seafoam","xkcd:weird green","xkcd:yellowgreen","xkcd:rust orange","xkcd:pale cyan","xkcd:pea",
              "xkcd:bright lime","xkcd:medium pink","xkcd:ugly purple","xkcd:hot magenta","xkcd:brown grey","xkcd:poo",
              "xkcd:water blue","xkcd:fern green","xkcd:dirty blue","xkcd:light sage","xkcd:hot green","xkcd:key lime",
              "xkcd:deep teal","xkcd:spearmint","xkcd:baby poo","xkcd:old pink","xkcd:seaweed green","xkcd:cool green",
              "xkcd:dark aquamarine","xkcd:dark grey blue","xkcd:pale sky blue","xkcd:light mustard","xkcd:barf green",
              "xkcd:brownish pink","xkcd:eggshell blue","xkcd:banana yellow","xkcd:lemon green","xkcd:lightish purple",
              "xkcd:dark seafoam green","xkcd:light grass green","xkcd:light mint green","xkcd:blue/grey","xkcd:royal",
              "xkcd:yellowy brown","xkcd:lavender blue","xkcd:twilight blue","xkcd:bright olive","xkcd:darkish purple",
              "xkcd:bubble gum pink","xkcd:greenish cyan","xkcd:turtle green","xkcd:sandy yellow","xkcd:duck egg blue",
              "xkcd:lipstick red","xkcd:purple brown","xkcd:brown yellow","xkcd:muddy yellow","xkcd:highlighter green",
              "xkcd:light greenish blue","xkcd:lightish red","xkcd:brick orange","xkcd:dusty red","xkcd:metallic blue",
              "xkcd:greenish turquoise","xkcd:deep rose","xkcd:greeny brown","xkcd:camouflage green","xkcd:baby green",
              "xkcd:piss yellow","xkcd:dusty orange","xkcd:purple/pink","xkcd:purpley","xkcd:dried blood","xkcd:putty",
              "xkcd:dark hot pink","xkcd:warm blue","xkcd:light khaki","xkcd:indian red","xkcd:dark cream","xkcd:camo",
              "xkcd:icky green","xkcd:greenblue","xkcd:dirty purple","xkcd:rich blue","xkcd:mushroom","xkcd:flat blue",
              "xkcd:cool grey","xkcd:canary","xkcd:booger green","xkcd:muted pink","xkcd:hazel","xkcd:dark royal blue",
              "xkcd:rich purple","xkcd:pale magenta","xkcd:light yellowish green","xkcd:indigo blue","xkcd:candy pink",
              "xkcd:orangeish","xkcd:light royal blue","xkcd:cocoa","xkcd:baby purple","xkcd:raw sienna","xkcd:burple",
              "xkcd:rosy pink","xkcd:peachy pink","xkcd:pale light green","xkcd:old rose","xkcd:fern","xkcd:dusk blue",
              "xkcd:light seafoam","xkcd:light neon green","xkcd:light bright green","xkcd:greyish teal","xkcd:spruce",
              "xkcd:dark slate blue","xkcd:dark sage","xkcd:coral pink","xkcd:true green","xkcd:tan brown","xkcd:azul",
              "xkcd:bright yellow green","xkcd:baby puke green","xkcd:poison green","xkcd:light lavendar","xkcd:toupe",
              "xkcd:dark yellow green","xkcd:greeny grey","xkcd:stormy blue","xkcd:purple/blue","xkcd:battleship grey",
              "xkcd:light light blue","xkcd:very pale blue","xkcd:vibrant blue","xkcd:darkish red","xkcd:dark fuchsia",
              "xkcd:light green blue","xkcd:dark blue grey","xkcd:deep sea blue","xkcd:reddy brown","xkcd:algae green",
              "xkcd:green again","xkcd:lemon lime","xkcd:violet red","xkcd:warm brown","xkcd:pale mauve","xkcd:claret",
              "xkcd:poop green","xkcd:off yellow","xkcd:grey/blue","xkcd:carnation","xkcd:pure blue","xkcd:shit green",
              "xkcd:bright cyan","xkcd:sunflower","xkcd:rusty red","xkcd:dandelion","xkcd:raw umber","xkcd:lawn green",
              "xkcd:dark sand","xkcd:greyblue","xkcd:wisteria","xkcd:red wine","xkcd:bluegrey","xkcd:sunflower yellow",
              "xkcd:blush pink","xkcd:twilight","xkcd:lipstick","xkcd:saffron","xkcd:butter","xkcd:petrol","xkcd:dark",
              "xkcd:pale salmon","xkcd:dark coral","xkcd:darkblue","xkcd:pastel red","xkcd:dark taupe","xkcd:neon red",
              "xkcd:blood orange","xkcd:gross green","xkcd:dodger blue","xkcd:sand yellow","xkcd:leather","xkcd:topaz",
              "xkcd:light mint","xkcd:ugly brown","xkcd:bluey grey","xkcd:flat green","xkcd:clay brown","xkcd:carmine",
              "xkcd:grapefruit","xkcd:warm pink","xkcd:sap green","xkcd:golden rod","xkcd:plum purple","xkcd:sapphire",
              "xkcd:light bluish green","xkcd:very dark brown","xkcd:orangey yellow","xkcd:electric lime","xkcd:creme",
              "xkcd:vivid purple","xkcd:racing green","xkcd:seafoam blue","xkcd:wintergreen","xkcd:purply","xkcd:rosa",
              "xkcd:violet pink","xkcd:pale peach","xkcd:green/blue","xkcd:yellow tan","xkcd:grey/green","xkcd:lichen",
              "xkcd:radioactive green","xkcd:light light green","xkcd:terracota","xkcd:tiffany blue","xkcd:clear blue",
              "xkcd:washed out green","xkcd:light pea green","xkcd:browny orange","xkcd:tealish green","xkcd:cinnamon",
              "xkcd:charcoal grey","xkcd:butterscotch","xkcd:burnt siena","xkcd:parchment","xkcd:pale gold","xkcd:ice",
              "xkcd:foam green","xkcd:light gold","xkcd:sand brown","xkcd:rust brown","xkcd:dusty teal","xkcd:manilla",
              "xkcd:greenish beige","xkcd:sickly yellow","xkcd:green/yellow","xkcd:very light brown","xkcd:sun yellow",
              "xkcd:sunny yellow","xkcd:eggplant purple","xkcd:baby shit brown","xkcd:darkish pink","xkcd:powder pink",
              "xkcd:kelley green","xkcd:browny green","xkcd:reddish grey","xkcd:bruise","xkcd:straw","xkcd:deep brown",
              "xkcd:liliac","xkcd:light grey green","xkcd:olive yellow","xkcd:purpley grey","xkcd:swamp","xkcd:cement",
              "xkcd:orange pink","xkcd:ugly blue","xkcd:warm purple","xkcd:nice blue","xkcd:tan green","xkcd:off blue",
              "xkcd:lavender pink","xkcd:dusty lavender","xkcd:desert","xkcd:deep lilac","xkcd:pig pink","xkcd:booger",
              "xkcd:bright sea green","xkcd:blue with a hint of purple","xkcd:dark green blue","xkcd:dark grass green",
              "xkcd:really light blue","xkcd:blue blue","xkcd:yellowish tan","xkcd:dark pastel green","xkcd:grey teal",
              "xkcd:green teal","xkcd:light moss green","xkcd:deep turquoise","xkcd:light eggplant","xkcd:cloudy blue",
              "xkcd:deep aqua","xkcd:nasty green","xkcd:strong pink","xkcd:darkgreen","xkcd:tea green","xkcd:custard",
              "xkcd:egg shell"]
        
    @staticmethod
    def lattice(name,genome,chromosomes):
        """
        Plots the NCSU simulator lattice Solution. Assumes 1/2 or 1/8 diagonal symmetry.
        Uses the xkcd color list for plotting genomes. This is a list of the 255 most 
        recognized colors by a public survey carried out by the xkcd webcomic. Google 
        it and it should be pretty obvious.
        """
        row = 1
        column = 0
        radius = 0.54
        pin_pitch = 1.26
        chrom_list = list(chromosomes.keys())

        total_rows = 1
        for gene in genome:
            column += 1
            if total_rows == column:
                column = 0
                total_rows += 1
        fig,ax = pyplot.subplots()
        column = 0
        used_list = []
        for gene in genome:
            gene_color = chrom_list.index(gene)
            if gene in used_list:
                pass
            else:
                gene.append(used_list)
                pyplot.plot([-10.00,-10.01],color=gene_color,label="{}".format(gene))
            circle = pyplot.Circle(([column*pin_pitch+pin_pitch/2.,(total_rows-row-1)*pin_pitch+pin_pitch/2.]),radius,color=Plotter.colors[gene_color])
            ax.add_artist(circle)
            column += 1
            if row == column:
                column = 0
                row += 1
        ax = pyplot.gca()
        ax.set_xlim((0,12))
        ax.set_ylim((0,12))
        ax.set_title(name)
        pyplot.savefig("{}_plot.png".format(name))

    @staticmethod
    def core_fixed_genomes(solution,genome_key):
        """
        Plots a 15 by 15 core lattice
        """
        loading_pattern_dictionary = {}
        for i in range(solution.core_width):
            loading_pattern_dictionary[i] = {}
            for j in range(solution.core_width):
                loading_pattern_dictionary[i][j] = "water"

        center = int(solution.core_width/2)
        if solution.core_width == 15:
            mapping = [[0,0], 
                       [1,0], [1,1], [1,2], [1,3], [1,4], [1,5], [1,6], [1,7],
                       [2,0], [2,1], [2,2], [2,3], [2,4], [2,5], [2,6],       
                       [3,0], [3,1], [3,2], [3,3], [3,4], [3,5], [3,6],       
                       [4,0], [4,1], [4,2], [4,3], [4,4], [4,5],              
                       [5,0], [5,1], [5,2], [5,3], [5,4],                     
                       [6,0], [6,1], [6,2], [6,3],                            
                       [7,0], [7,1]]

        count = 0
        break_count = len(mapping)
        assembly_list = ['water','Burned']
        for i,gene in enumerate(solution.genome):
            rfromc = mapping[i][0]
            cfromc = mapping[i][1]

            if type(genome_key[gene]['assembly']) == list:
                pass
            else:
                if genome_key[gene]['assembly'] in assembly_list:
                    pass
                else:
                    assembly_list.append(genome_key[gene]['assembly'])

            if rfromc == 0:
                if type(genome_key[gene]['assembly']) == list:
                    loading_pattern_dictionary[center][center] = 'Burned' 
                else:
                    loading_pattern_dictionary[center][center] = genome_key[gene]['assembly']
            elif cfromc == 0:
                if type(genome_key[gene]['assembly']) == list:
                    loading_pattern_dictionary[center][center - rfromc] = "Burned"
                    loading_pattern_dictionary[center + rfromc][center] = "Burned"
                    loading_pattern_dictionary[center][center + rfromc] = "Burned"
                    loading_pattern_dictionary[center - rfromc][center] = "Burned"
                else:
                    loading_pattern_dictionary[center][center - rfromc] = genome_key[gene]['assembly']
                    loading_pattern_dictionary[center + rfromc][center] = genome_key[gene]['assembly']
                    loading_pattern_dictionary[center][center + rfromc] = genome_key[gene]['assembly']
                    loading_pattern_dictionary[center - rfromc][center] = genome_key[gene]['assembly']
            else:
                if type(genome_key[gene]['assembly']) == list:
                    loading_pattern_dictionary[center - rfromc][center - cfromc] = "Burned"
                    loading_pattern_dictionary[center + rfromc][center + cfromc] = "Burned"
                    loading_pattern_dictionary[center - rfromc][center + cfromc] = "Burned"
                    loading_pattern_dictionary[center + rfromc][center - cfromc] = "Burned"
                else:
                    loading_pattern_dictionary[center - rfromc][center - cfromc] = genome_key[gene]['assembly']
                    loading_pattern_dictionary[center + rfromc][center + cfromc] = genome_key[gene]['assembly']
                    loading_pattern_dictionary[center - rfromc][center + cfromc] = genome_key[gene]['assembly']
                    loading_pattern_dictionary[center + rfromc][center - cfromc] = genome_key[gene]['assembly']
            count +=1
            if count >= break_count:
                break

        figure,ax = pyplot.subplots()
        for i in range(solution.core_width):
            for j in range(solution.core_width):
                plot_color = Plotter.colors[assembly_list.index(loading_pattern_dictionary[i][j])]
                rect = pyplot.Rectangle((i,j),1,1,color=plot_color,label=loading_pattern_dictionary[i][j])
                ax.add_artist(rect)
        pyplot.savefig("{}_plot.png".format(solution.name))

class Simulated_Annealing_Metric_Toolbox(object):
    """
    Class for recording the metrics from simulated anealing. However Johnny tracked the files
    didn't make sense/didn't get correctly updated to Github, so I am doing it the way I would
    have made him do it anyways.
    """
    def __init__(self):
        self.best_solution_over_time = {}
        self.other_solution_over_time = {}
        self.best_solution_over_time['fitness'] = []
        self.other_solution_over_time['fitness'] = []
        self.temperature_over_time = []
        self.all_param_over_time={}
        self.acceptance_probability_over_time = []
        self.generation = 0

    def record_best_and_new_solution(self,best,other,cooling_schedule):
        """
        Records the parameters and objective values for the current best and
        new solution.

        NOTE: This record is performed before the annealing process takes place,
        meaning that changes in the best solution to the optimization problem should
        be seen one step after being recorded in the other solution column. 

        Cooling schedule is also recorded so that overall probability can be judged.
        """
        for param in best.parameters:
            if param in self.best_solution_over_time:
                foo = best.parameters[param]['value']
                self.best_solution_over_time[param].append(foo)
            else:
                foo = best.parameters[param]['value']
                self.best_solution_over_time[param] = [foo]
        for param in other.parameters:
            if param in self.other_solution_over_time:
                foo = other.parameters[param]['value']
                self.other_solution_over_time[param].append(foo)
            else:
                foo = other.parameters[param]['value']
                self.other_solution_over_time[param] = [foo]
        self.best_solution_over_time['fitness'].append(best.fitness)
        self.other_solution_over_time['fitness'].append(other.fitness)
        acceptance = numpy.exp(-1*(other.fitness-best.fitness)/cooling_schedule.temperature)
        if acceptance > 1.:
            acceptance = 1.
        self.acceptance_probability_over_time.append(acceptance)
        self.temperature_over_time.append(cooling_schedule.temperature)
        self.all_param_over_time[self.generation]={}
        track_file = open('optimization_track_file.txt','a')
        track_file.write(f"Current Generation {self.generation}\n")
        track_file.write(f"Temperature {cooling_schedule.temperature}\n")
        track_file.write(f"Acceptance probability {acceptance}\n")
        track_file.write(f"{best.genome}\n")
        track_file.write(f"{other.genome}\n")
        track_file.write(f"Current best solution, {best.name}, Other solution, {other.name}\n")
        track_file.write(f"Fitness                {round(best.fitness,3)},            ,{round(other.fitness,3)}\n")
        for param in best.parameters:
            track_file.write(f"{param}, {best.parameters[param]['value']},            ,{other.parameters[param]['value']}   \n")
            self.all_param_over_time[self.generation][param]=other.parameters[param]['value']
        track_file.write("\n\n\n")
        track_file.close()
        self.generation += 1    

    def record_all_param(self,population_class,generation_class,flag):
        self.all_param_over_time[generation_class.current] = {}
        for sol in population_class.parents:
            param_list = [param for param in sol.parameters]
            for param in sol.parameters:
                self.all_param_over_time[generation_class.current][param]=[]
            break 
        # appending values
        if flag:
            for sol in population_class.parents:
                if sol.name.find('initial') < 0:
                    for param in sol.parameters:
                        tmp = sol.parameters[param]['value']
                        self.all_param_over_time[generation_class.current][param].append(tmp)
        for sol in population_class.children:
            if sol.name.find('initial') < 0:
                for param in sol.parameters:
                    tmp = sol.parameters[param]['value']
                    self.all_param_over_time[generation_class.current][param].append(tmp)
        return

    def make_gif(self, frame_folder, inpstr):
        frames = [Image.open(image) for image in sorted(glob.glob(f"{frame_folder}/*.jpg")\
            ,key=lambda x: self.get_int(x))]
        # extract first image from iterator
        frame_one = frames[0]
        frame_one.save("pareto_front_"+inpstr+".gif", format="GIF", \
                        append_images=frames, save_all=True, duration=500, loop=1)

    def get_int(self, str):
        tmp = str.split('_')
        tmp = tmp[-1].split('.')
        num = int(tmp[0])
        return num

    def gen_plot(self, target, remain_param, c):
        ''' Generate corvar plot 
        target      : target parameter for covar plot
        remain_param: remain list of param 
        c           : color list           
        '''

        for param in remain_param:
            os.makedirs('image_dir_'+target+'_'+param, exist_ok=True)
            path = os.getcwd() + "/image_dir_"+target+"_"+param
            for gen in range(len(self.all_param_over_time)): # loop generation 
                fig,ax = pyplot.subplots(2,1, figsize=(5,5),  gridspec_kw={'height_ratios': [25,1]},  constrained_layout = True)
                for t in range (gen+1):
                    sc= ax[0].scatter(self.all_param_over_time[t][param],self.all_param_over_time[t][target] , s = 50, color=c[t])
                ax[0].set_xlabel(str(param))
                ax[0].set_ylabel(str(target))
                cmap = mpl.colors.ListedColormap(c)
                cmap.set_over('0.25')
                cmap.set_under('0.75')
                bounds = range(1,len(self.all_param_over_time)+2)
                norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
                cb2 = mpl.colorbar.ColorbarBase(ax[1], cmap=cmap,
                                                norm=norm,
                                                #boundaries=[0] + bounds + [13],
                                                extend='both',
                                                ticks=[bounds[0],max(bounds)],
                                                spacing='proportional',
                                                orientation='horizontal')
                cb2.set_label('Generation')
                name = path +'/image_'+str(gen+1)+'.jpg'
                pyplot.savefig(name)
                pyplot.close()
            self.make_gif(path+'/', param+'_'+target )

    def plotter(self):
        """
        Plots the best, worst and average fitness and parameter values for
        the optimization.
        """
        pyplot.figure()
        pyplot.plot(range(len(self.acceptance_probability_over_time)),self.acceptance_probability_over_time)
        pyplot.title("Acceptance Probability")
        pyplot.xlabel("Generation")
        pyplot.ylabel("Acceptance Probability")
        pyplot.savefig("acceptance_probability.png")
        pyplot.close()

        pyplot.figure()
        pyplot.plot(range(len(self.temperature_over_time)),self.temperature_over_time)
        pyplot.title("Temperature over time")
        pyplot.xlabel("Generation")
        pyplot.ylabel("Temperature")
        pyplot.savefig("temperature.png")
        pyplot.close()
        
        for param in self.best_solution_over_time:
            pyplot.figure()
            foo1 = range(len(self.best_solution_over_time[param]))
            foo2 = range(len(self.other_solution_over_time[param]))
            pyplot.plot(foo1,self.best_solution_over_time[param],label='Best')
            pyplot.plot(foo2,self.other_solution_over_time[param],label='Other')

            pyplot.xlabel('Generation')
            pyplot.ylabel(param)
            pyplot.title("{} for Best and Other Solution over Time".format(param))
            pyplot.legend()
            pyplot.savefig("{}_over_time".format(param))            
            pyplot.close()
        # pareto front graph is not supported in SA
        list_param = [param for param in self.best_solution_over_time]
        list_param.pop(list_param.index('fitness')) # no fitness pareto front graph
        while list_param:
            get_colors = lambda n: list(map(lambda i: "#" + "%06x" % random.randint(0, 0xFFFFFF),range(n)))
            color = get_colors(len(self.all_param_over_time)) # sample return:  ['#8af5da', '#fbc08c', '#b741d0', '#e599f1', '#bbcb59', '#a2a6c0']
            target = list_param.pop(0)
            self.gen_plot(target, list_param, color)

class Lava_Metric_Toolbox(object):
    """
    Toolbox for plotting the lava flow optimization. Not sure if this is really
    necessary, or can be replaced with existing tools.
    """
    def __init__(self):
        self.best_flow_over_time = {}
        self.best_flow_over_time['fitness'] = []
        self.generation_list = []

    def update_best_values(self,flow,gen,count,solution):
        """
        Updates the dictionary of the best values of the flow.
        Not really going to record the worst and average values, because they
        don't necessarily mean anything in this optimization.
        """
        for param in solution.parameters:
            if param in self.best_flow_over_time:
                foo = solution.parameters[param]['value']
                self.best_flow_over_time[param].append(foo)
            else:
                foo = solution.parameters[param]['value']
                self.best_flow_over_time[param] = [foo]
        self.best_flow_over_time['fitness'].append(solution.fitness)
        self.generation_list.append(gen)
        track_file = open(f'Flow_{flow}_track_file.txt','a')
        track_file.write(f"{count},   {gen},   {solution.name},   ")
        for param in solution.parameters:
            track_file.write(f"{solution.parameters[param]['value']},   ")
        track_file.write("\n\n")
        track_file.close()

    def plotter(self):
        """
        Plots the best_flow_dictionary for every parameter.
        """
        for param in self.best_flow_over_time:
            pyplot.figure()
            pyplot.plot(self.generation_list,self.best_flow_over_time[param])
            pyplot.xlabel("Flow Number")
            pyplot.ylabel(param)
            pyplot.title(f"{param} over time for best solution")
            pyplot.savefig(f"{param}_over_time")
            pyplot.close()




