import os
import sys
import copy
import math
import numpy
import random
import statistics
import shutil
import time
from multiprocessing import Pool
from midas.utils.solution_types import evaluate_function
from midas.utils.metrics import Simulated_Annealing_Metric_Toolbox
from midas import logo, version
import multiprocessing


"""
This file contains a few of the classes and methods involved with 
performing an optimization via simulated aneealing.
NOTE: This file does not include the means for reproduction, but is 
written as imf the class is included. 
Written by Johnny Klemes. 3/19/2020
"""


class Cooling_Schedule(object):
    """
    Class for Simulated Annealing cooling schedules.

    THe cooling schedule sets the tolerance for accepting new solutions.
    A high initial temperature indicates accepting a more new designs even if
    they have a less favorable objective function. Thus the logarithmic cooling
    schedule is favorable for this problem.

    There are two cooling schedules defined below. In both the temperature is
    determined by the current era of a lifetime, represented by a piecewise
    function. THe second cooling schedule is identical to the first but with
    two cycles.

    """

    def __init__(self, generation):
        self.generation = generation

    @staticmethod
    def piecewise(gen, lifet):
        if gen <= 3 / 12 * lifet:
            temp = gen
            y = 2 * numpy.log(temp / (4 / 12 * lifet))
        elif gen <= 4 / 12 * lifet and gen > 3 / 12 * lifet:
            y = 1000
        elif gen > 4 / 12 * lifet and gen <= 6 / 12 * lifet:
            temp = gen - 4 / 12 * lifet
            y = numpy.log(0.8 * temp / (3 / 12 * lifet))
        elif gen > 6 / 12 * lifet and gen <= 7 / 12 * lifet:
            y = 1000
        elif gen > 7 / 12 * lifet and gen <= 10 / 12 * lifet:
            temp = gen - 7 / 12 * lifet
            y = 3 * numpy.log(0.8 * temp / (4 / 12 * lifet))
        else:
            y = 1000

        return y

    @staticmethod
    def twicewise(gen, lifet):
        if gen <= 3 / 24 * lifet:
            temp = gen
            y = 2 * numpy.log(temp / (4 / 24 * lifet))
        elif gen <= 4 / 24 * lifet and gen > 3 / 24 * lifet:
            y = 1000
        elif gen > 4 / 24 * lifet and gen <= 6 / 24 * lifet:
            temp = gen - 4 / 24 * total
            y = numpy.log(0.8 * temp / (3 / 24 * lifet))
        elif gen > 6 / 24 * lifet and gen <= 7 / 24 * lifet:
            y = 1000
        elif gen > 7 / 24 * lifet and gen <= 10 / 24 * lifet:
            temp = gen - 7 / 24 * lifet
            y = 3 * numpy.log(0.8 * temp / (4 / 24 * lifet))
        elif gen > 10 / 24 * lifet and gen <= 1 / 2 * lifet:
            y = 1000
        elif gen > 1 / 2 * lifet and gen <= 15 / 24 * lifet:
            temp = gen - 1 / 2 * lifet
            y = 2 * numpy.log(temp / (4 / 24 * lifet))
        elif gen <= 16 / 24 * lifet and gen > 15 / 24 * lifet:
            y = 1000
        elif gen <= 19 / 24 * lifet and gen > 16 / 24 * lifet:
            temp = gen - 16 / 24 * lifet
            y = 3 * numpy.log(0.8 * temp / (4 / 24 * lifet))
        else:
            y = 1000

        return y


class Exponential_Decreasing_Cooling_Schedule(object):
    """
    Logarithmic cooling schedule for simulated annealing. Implementing this because Johnny Klemes cooling schedules seem broken.
    His implementation that I pulled from Github is broken at the least, and there isn't sufficient documentation available to
    understand how to fix it.
    This cooling schedule is about as simple as it can get.
    T = T0*alpha
    Where 0.9 < alpha < 1.0 and 1 < T0 < 10
    """

    def __init__(self, alpha, temperature):
        self.alpha = 1. - alpha
        self.temperature = temperature

    def update(self):
        """
        Updates the temperature of the cooling schedule.
        """
        if self.temperature <= 0.0001:
            self.temperature = 0.0001
        else:
            self.temperature *= self.alpha


def SA(x, k, active, Buffer, BufferCost, self, opt):
    """
    Created by Jake Mikouchi
    2/10/23
    """

    def UpdateActive(Buffer, BufferCost):
        # finds the sum of the probability of each position in the Buffer
        scalingfactor = 100
        SumProb = 0
        for i in range(len(BufferCost)):
            # CurrProb = numpy.exp((-1 * BufferCost[i]) / temp)
            CurrProb = numpy.exp((-1 * BufferCost[i]) / scalingfactor)

            SumProb += CurrProb

        # finds the probability that each position of the buffer will be selected
        PositionProb = []
        for i in range(len(BufferCost)):
            # CurrProb = numpy.exp((-1 * BufferCost[i]) / temp)
            CurrProb = numpy.exp((-1 * BufferCost[i]) / scalingfactor)

            TotalProb = (CurrProb / SumProb)

            PositionProb.append(TotalProb)

        # uses the probabilities to select which lists to add to the buffer
        SumProblist = []
        for i in range(len(Buffer)):
            if i >= 1:
                SumProblist.append(SumProblist[i - 1] + PositionProb[i])
            if i == 0:
                SumProblist.append(PositionProb[i])

        randomnum = random.uniform(0, 1)

        for i in range(len(Buffer)):
            if i == 0:
                if randomnum < SumProblist[i]:
                    acceptedlist = Buffer[i]
                    acceptedlistCost = BufferCost[i]
                    break
            if i >= 1:
                if randomnum < SumProblist[i] and randomnum >= SumProblist[i - 1]:
                    acceptedlist = Buffer[i]
                    acceptedlistCost = BufferCost[i]
                    break

        return acceptedlist, acceptedlistCost

    if x == 0:
        active = active
    else:
        active, activeCost = UpdateActive(Buffer, BufferCost)

    # activeindex = BufferCost.index(min(BufferCost))
    # active = Buffer[activeindex]
    # activeCost = BufferCost[activeindex]

    temp = self.cooling_schedule.temperature
    PAR = 0
    PAR2 = 0

    one = []
    one.append(active)

    one = self.fitness.calculate(one)

    solutions = []
    solutionsfitness = []
    NewSolutionSA = [active]
    NewSolutionCost = [active.fitness]
    for number in range(self.file_settings['optimization']['population_size']):
        challenge = self.solution()
        if str(self.solution) == "<class 'midas.applications.parcs_332.MCycle_Loading_Pattern_Solution'>" or str(self.solution) == "<class 'midas.applications.parcs_332.MCycle_Grouped_Loading_Pattern_Solution'>" or str(self.solution) == "<class 'midas.applications.parcs_332.MCycle_Inventory_Loading_Pattern_Solution'>" or str(self.solution) == "<class 'midas.applications.parcs_343.MCycle_Inventory_Loading_Pattern_Solution'>":
                    challenge.genome = active.reproduce()
        else:
            challenge.genome = self.mutation.reproduce(active.genome)
        challenge.name = f"child_{x}_{k}_{number}"
        challenge.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
        challenge.add_additional_information(self.file_settings)
        challenge.evaluate()
        all_values = open('all_value_tracker.txt','a')
        all_values.write(f"{challenge.name},    ")
        for param in challenge.parameters:
            all_values.write(f"{challenge.parameters[param]['value']},    ")
        all_values.write('\n')
        all_values.close()

        test = [challenge]
        test = self.fitness.calculate(test)

        # determining which solution makes the next generation
        acceptance = numpy.exp(-1 * (challenge.fitness - active.fitness) / temp)
        if challenge.fitness < active.fitness:
            PAR += 1
            PAR2 += 1
            active = copy.deepcopy(challenge)

        elif random.uniform(0, 1) < acceptance:
            PAR += 1
            PAR2 += acceptance
            active = challenge

        if x <= self.file_settings['optimization']['number_of_generations'] / 2:
            temp = temp * 0.9

        NewSolutionSA.append(challenge)
        NewSolutionCost.append(challenge.fitness)
        solutions.append(active)
        solutionsfitness.append(active.fitness)

    return (solutions, solutionsfitness, NewSolutionCost, NewSolutionSA, PAR, PAR2)

def SA_prun(k, self):
    # creates a single initial solution
    active = self.solution()
    # active.genome = self.mutation.reproduce(active.genome)
    

    active.name = f"initial_temp_calc" + str(k)
    active.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
    # active.genome = self.mutation.reproduce(active.genome)
    active.add_additional_information(self.file_settings)
    if active.fixed_genome:
        active.generate_initial_fixed(self.file_settings['genome']['chromosomes'],
                                        self.file_settings['optimization']['fixed_groups'])
    else:
        active.generate_initial(self.file_settings['genome']['chromosomes'])

    active.evaluate()
    all_values = open('all_value_tracker.txt','a')
    all_values.write(f"{active.name},    ")
    for param in active.parameters:
        all_values.write(f"{active.parameters[param]['value']},    ")
    all_values.write('\n')
    all_values.close()
    
    one = []
    one.append(active)
    one = self.fitness.calculate(one)

    print('calculation {}, fitness = {}'.format(k,active.fitness))

    return (active)

def SetInitial(self):
    """
    Created by Jake Mikouchi
    2/10/23
    """
    # ---------------------------------------------------------------------------------------------------
    Buffer = []
    for w in range(self.file_settings['optimization']['buffer_length']):
        active = self.solution()

        active.name = f"start_solution_{w}"
        active.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
        active.add_additional_information(self.file_settings)

        active.generate_initial(self.file_settings['genome']['chromosomes'])
        active.evaluate()
        Buffer.append(active)
    Buffer = self.fitness.calculate(Buffer)
    BufferCost = [Buffer[i].fitness for i in range(len(Buffer))]

    # a > 1.0 and a < 2.0
    a = 1.5
    deviation = statistics.stdev(BufferCost)
    # find standard deviation of the costs and caluclate
    # the initaial temp from the standard deviation
    initialtemp = a * deviation
    # -----------------------------------------------------------------------------------------------------
    return (Buffer, BufferCost, initialtemp)


class SimulatedAnnealing(object):
    def __init__(self, solution,
                 population,
                 generation,
                 cooling_schedule,
                 fitness,
                 mutation,
                 num_procs,
                 file_settings
                 ):

        self.solution = solution
        self.population = population
        self.generation = generation
        self.cooling_schedule = cooling_schedule
        self.fitness = fitness
        self.mutation = mutation
        self.num_procs = num_procs
        self.file_settings = file_settings
        self.number_generations_post_cleanup = 300  # Arbitrarily chosen default.
        if 'cleanup' in file_settings['optimization']:
            if file_settings['optimization']['cleanup']['perform']:
                self.perform_cleanup = True
            else:
                self.perform_cleanup = False
        else:
            self.perform_cleanup = False

        if self.num_procs > 1:
            omethod = "Parallel Simulated Annealing"
        else:
            omethod = "Simulated Annealing"
        f = open(version.__ofile__, 'a')
        f.write('\n\n---------------------------------------GENERAL INFORMATION----------------------------------------')
        f.write('\nOptimization Method: {}'.format(omethod))
        f.write('\nProcessors: {}'.format(self.num_procs))
        f.write('\nTemparture Updates: {}'.format(file_settings['optimization']['number_of_generations']))
        f.write('\nSA Iterations: {}'.format(file_settings['optimization']['population_size']))
        f.write('\nPSA Iterations: {}'.format(self.num_procs))
        f.write('\nBuffer Size: {}'.format(file_settings['optimization']['buffer_length']))
        f.write('\nObjectives:')
        for key, value in file_settings['optimization']['objectives'].items():
            if 'target' not in value:
                f.write('\n  Name: {},  Goal: {},  Weight:  {}'.format(key,value['goal'],value['weight']))
            else:
                f.write('\n  Name: {},  Goal: {},  Target:  {},  Weight:  {}'.format(key,value['goal'],value['target'],value['weight']))
        f.write('\n--------------------------------------------------------------------------------------------------')
        f.close()

    def main_in_serial(self):
        """
        Performs optimization using simulated annealing for a population size of one
        """
        opt = Simulated_Annealing_Metric_Toolbox()

        track_file = open('optimization_track_file.txt', 'w')
        track_file.write("Beginning Optimization \n")
        track_file.close()

        # defining the active solution
        active = self.solution()
        active.name = "initial_solution"
        active.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
        active.add_additional_information(self.file_settings)
        if active.fixed_genome:
            active.generate_initial_fixed(self.file_settings['genome']['chromosomes'],
                                          self.file_settings['optimization']['fixed_groups'])
        else:
            active.generate_initial(self.file_settings['genome']['chromosomes'])
        # print(active.genome)
        active.evaluate()
        # will store the fit for active solution
        one = []
        one.append(active)

        one = self.fitness.calculate(one)

        for self.generation.current in range(self.generation.total):
            for number in range(self.population.size):
                challenge = self.solution()
                if str(self.solution) == "<class 'midas.applications.parcs_332.MCycle_Loading_Pattern_Solution'>" or str(self.solution) == "<class 'midas.applications.parcs_332.MCycle_Grouped_Loading_Pattern_Solution'>" or str(self.solution) == "<class 'midas.applications.parcs_332.MCycle_Inventory_Loading_Pattern_Solution'>"  or str(self.solution) == "<class 'midas.applications.parcs_343.MCycle_Inventory_Loading_Pattern_Solution'>":
                    challenge.genome = active.reproduce()
                else:
                    challenge.genome = self.mutation.reproduce(active.genome)
                challenge.name = "solution_{}_{}".format(self.generation.current, number)
                challenge.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
                challenge.add_additional_information(self.file_settings)
                challenge.evaluate()

                test = [challenge]
                test = self.fitness.calculate(test)

                # determining which solution makes the next generation
                opt.record_best_and_new_solution(active, challenge, self.cooling_schedule)
                acceptance = numpy.exp(-1 * (challenge.fitness - active.fitness) / self.cooling_schedule.temperature)
                if challenge.fitness < active.fitness:
                    active = copy.deepcopy(challenge)
                elif random.uniform(0, 1) < acceptance:
                    active = challenge
            self.cooling_schedule.update()

        track_file = open('optimization_track_file.txt', 'a')
        track_file.write("End of Optimization \n")
        track_file.close()

        # plot the parameters over time
        opt.plotter()
        self.cleanup()

    def main_in_parallel(self):
        def UpdateBuffer(Buffer, BufferCost, TSolutions, TSolutionsfitness):

            lenTSolutions = len(TSolutions)
            NewBuffer = []
            NewBuffCost = []

            #-----Gets rid of duplicates in Tsolutions and Buffer -------------------------------------------------
            DummyTSolutionsCost = []
            DummyTSolutions = []
            for x in TSolutionsfitness:
                minTSolutionfitnessindex = TSolutionsfitness.index(x)
                if x not in DummyTSolutionsCost:
                    DummyTSolutionsCost.append(x)
                    DummyTSolutions.append(TSolutions[minTSolutionfitnessindex])
            # replaces TSolutions with a new list with same contents except no duplicates
            TSolutions = DummyTSolutions
            TSolutionsfitness = DummyTSolutionsCost

            DummyBufferCost = []
            DummyBuffer = []
            for x in BufferCost:
                minBufferfitnessindex = BufferCost.index(x)
                if x not in DummyBufferCost:
                    DummyBufferCost.append(x)
                    DummyBuffer.append(Buffer[minBufferfitnessindex])
            if len(DummyBuffer) < self.file_settings['optimization']['buffer_length']:
                while len(DummyBuffer) < self.file_settings['optimization']['buffer_length']:
                    # this solution may be crude but it saves me a lot of coding
                    DummyBuffer.append(100000)
                    DummyBufferCost.append(100000)
            # replaces Buffer with a new list with same contents except no duplicates 
            Buffer = DummyBuffer
            BufferCost = DummyBufferCost
            # --------------------------------------------------------------------------------------------------

            if len(TSolutions) >= self.file_settings['optimization']['buffer_length']:
                for w in range(self.file_settings['optimization']['buffer_length']):
                    minTSolutionfitness = min(TSolutionsfitness)
                    minTSolutionfitnessindex = TSolutionsfitness.index(min(TSolutionsfitness))
                    minBufferfitness = min(BufferCost)
                    minBufferfitnessindex = BufferCost.index(min(BufferCost))

                    if minTSolutionfitness <= minBufferfitness:
                        NewBuffCost.append(minTSolutionfitness)
                        NewBuffer.append(TSolutions[minTSolutionfitnessindex])
                        TSolutionsfitness.pop(minTSolutionfitnessindex)
                        TSolutions.pop(minTSolutionfitnessindex)

                    if minTSolutionfitness > minBufferfitness:
                        NewBuffCost.append(minBufferfitness)
                        NewBuffer.append(Buffer[minBufferfitnessindex])
                        BufferCost.pop(minBufferfitnessindex)
                        Buffer.pop(minBufferfitnessindex)

            if lenTSolutions < self.file_settings['optimization']['buffer_length']:
                for w in range(len(TSolutions)):
                    Buffer.pop(BufferCost.index(max(BufferCost)))
                    BufferCost.pop(BufferCost.index(max(BufferCost)))
                Buffer.extend(TSolutions)
                BufferCost.extend(TSolutionsfitness)
                NewBuffer = Buffer
                NewBuffCost = BufferCost

            return NewBuffer, NewBuffCost

        def UpdateActive(temp, Buffer, BufferCost):
            # finds the sum of the probability of each position in the Buffer
            scalingfactor = 100
            SumProb = 0
            for i in range(len(BufferCost)):
                CurrProb = numpy.exp((-1 * BufferCost[i]) / scalingfactor)
                SumProb += CurrProb

            # finds the probability that each position of the buffer will be selected
            PositionProb = []
            for i in range(len(BufferCost)):
                CurrProb = numpy.exp((-1 * BufferCost[i]) / scalingfactor)
                TotalProb = (CurrProb / SumProb)
                PositionProb.append(TotalProb)

            # uses the probabilities to select which lists to add to the buffer
            SumProblist = []
            for i in range(len(Buffer)):
                if i >= 1:
                    SumProblist.append(SumProblist[i - 1] + PositionProb[i])
                if i == 0:
                    SumProblist.append(PositionProb[i])

            randomnum = random.uniform(0, 1)

            for i in range(len(Buffer)):
                if i == 0:
                    if randomnum < SumProblist[i]:
                        acceptedlist = Buffer[i]
                        acceptedlistCost = BufferCost[i]
                        break
                if i >= 1:
                    if randomnum < SumProblist[i] and randomnum >= SumProblist[i - 1]:
                        acceptedlist = Buffer[i]
                        acceptedlistCost = BufferCost[i]
                        break

            return acceptedlist, acceptedlistCost

        def InitialSolution(temp, Buffer, BufferCost):
            # finds the sum of the probability of each position in the Buffer
            SumProb = 0
            for i in range(len(Buffer)):
                CurrProb = numpy.exp((-1 * BufferCost[i]) / temp)
                SumProb += CurrProb

            # finds the probability that each position of the buffer will be selected
            PositionProb = []
            for i in range(len(Buffer)):
                CurrProb = numpy.exp((-1 * BufferCost[i]) / temp)
                TotalProb = (CurrProb / SumProb)
                PositionProb.append(TotalProb)

            # uses the probabilities to select which lists to add to the buffer
            SumProblist = []
            for i in range(len(Buffer)):
                if i >= 1:
                    SumProblist.append(SumProblist[i - 1] + PositionProb[i])
                if i == 0:
                    SumProblist.append(PositionProb[i])
            randomnum = numpy.random.rand()
            for i in range(len(Buffer)):
                if i == 0:
                    if randomnum < SumProblist[i]:
                        acceptedlist = Buffer[i]
                        break
                if i >= 1:
                    if randomnum < SumProblist[i] and randomnum > SumProblist[i - 1]:
                        acceptedlist = Buffer[i]
                        break

            return acceptedlist


        def InitialTemp(self,ctx):

            # calculates the initial temperature based on the standard deviations
            # of costs when the probability is 100%
            Buffer = []
            costs = []
            # a > 1.0 and a < 2.0
            a = 1.5
            ninit = self.file_settings['optimization']['buffer_length']
            with ctx.Pool(processes=self.num_procs, ) as p:
                data = p.starmap(SA_prun, [(k, self) for k in range(ninit)])

            for active in data:
                costs.append(active.fitness)
                Buffer.append(active)

            deviation = statistics.stdev(costs)

            # find standard deviation of the costs and caluclate
            # the initaial temp from the standard deviation
            initialtemp = a * deviation

            return (initialtemp, Buffer)

        # Lam cooling schedule
        def LAM(currtemp, deviation, MoveAcceptanceRatio):
            # quality > 1 and quality < 2
            p = MoveAcceptanceRatio
            if p == 1:
                p = 0.9
            qualityfactor = 1.1
            # calculate Gp
            G1 = (4 * p * ((1 - p) ** 2))
            G2 = ((2 - p) ** 2)
            Gp = G1 / G2
            # calculate temperature
            sk = (1 / currtemp)
            sk1 = sk + (qualityfactor * (1 / deviation) * (1 / ((sk ** 2) * (deviation ** 2)))) * Gp
            temp = 1 / sk1

            return temp

        def LAMStats(temp, BufferCostlist):
            STDCostlist = BufferCostlist

            # find standard deviation of cost at current temp
            if not STDCostlist:
                deviation = 1
            else:
                if len(STDCostlist) == 1:
                    deviation = 1
                else:
                    if statistics.stdev(STDCostlist) == 0:
                        deviation = 0.01
                    else:
                        deviation = statistics.stdev(STDCostlist)
            return deviation
        
        def BufferPrint(Buffer):
            f = open(version.__ofile__, 'a')
            f.write('\n\n------------------------------------CURRENT BUFFER RESULTS -------------------------------------')
            for i in range(len(Buffer)):
                buffer_solution = Buffer[i]
                f.write('\nSolution {}:'.format(buffer_solution.name))
                for key, value in buffer_solution.parameters.items():
                    f.write('\n  {}: {}'.format(key,value['value']))
                f.write('\n  fitness: {}'.format(buffer_solution.fitness))
                for key, value in buffer_solution.core_dict['fuel'].items():
                    f.write('\n  {}: '.format(key))
                    for kkey, vvalue in value.items():
                        f.write('{}: {}, '.format(kkey, vvalue['Value']))
            f.write('\n--------------------------------------------------------------------------------------------------')
            f.close()

        def PSAPrint(print_dict):
            f = open(version.__ofile__, 'a')
            f.write('\n\n-----------------------------------PSA ITERATION {} RESULTS ------------------------------------'.format(print_dict['iteration']))
            f.write('\nTemperature: {}'.format(print_dict['temperature']))
            f.write('\nMove Acceptance Ratio: {}'.format(print_dict['move_acceptance_ratio']))
            f.write('\nStandard Deviation: {}'.format(print_dict['standard_deviation']))
            f.write('\nTemperature: {}'.format(print_dict['temperature_update']))
            f.write('\n--------------------------------------------------------------------------------------------------')
            f.close()

        SetStart = 0
        # SetStart determines if initial buffer will be set or completely random
        # if SetStart is 0, then the buffer will be random
        # if SetStart is 1, Buffer will be predetemined

        multiprocessing.set_start_method("spawn")
        ctx = multiprocessing.get_context('spawn')

        all_value_count = 0
        all_values = open('all_value_tracker.txt','w')
        all_values.close()

        f = open(version.__ofile__, 'a')
        f.write('\n\nRUNNING BUFFER INITIALIZATION...')
        f.close()

        start = time.time()
        if SetStart == 0:
            initial_temp, Buffer = InitialTemp(self,ctx)
            Buffer = self.fitness.calculate(Buffer)
            BufferCost = [Buffer[i].fitness for i in range(len(Buffer))]

        if SetStart == 1:
            Buffer, BufferCost, initial_temp = SetInitial(self)
            self.cooling_schedule.temperature = initial_temp
        end = time.time()

        opt = Simulated_Annealing_Metric_Toolbox()
        BufferPrint(Buffer)
        f = open(version.__ofile__, 'a')
        f.write('\n\nBUFFER INITIALIZATION COMPLETED IN {}s'.format(end-start))
        f.close()

        track_file = open('optimization_track_file.txt', 'w')
        track_file.write("Beginning Optimization \n")
        track_file.close()

        active = InitialSolution(self.cooling_schedule.temperature, Buffer, BufferCost)
        BestSolution = active
        BestSolutionCost = active.fitness

        
        for x in range(self.file_settings['optimization']['number_of_generations']):
            f = open(version.__ofile__, 'a')
            f.write('\n\nRUNNING PSA TEMPERATURE ITERATION {}...'.format(x+1))
            f.close()
            print_dict={}
            print_dict['iteration']=x+1
            print_dict['temperature']=self.cooling_schedule.temperature
            
            start = time.time()
            with ctx.Pool(processes=self.num_procs, ) as p:
                data = p.starmap(SA, [(x, k, active, Buffer, BufferCost, self, opt) for k in range(self.num_procs)])

            # data collection
            # create list of the costs of the values from data
            TSolutions = []
            TSolutionsfitness = []
            NewSolutionsfitness = []
            NewSolutions = []

            TotalMoves = 0
            TotalAcceptanceProbability = 0
            for i in range(len(data)):
                TSolutions.extend(data[i][0])
                TSolutionsfitness.extend(data[i][1])
                NewSolutions.extend(data[i][3])
                NewSolutionsfitness.extend(data[i][2])
                TotalMoves += data[i][4]
                TotalAcceptanceProbability += data[i][5]

            # determines move Move Acceptance Method
            # if 0 move acceptance is determined by total number of times a move is made by probability
            # if 1 move acceptance is determined by the sum of all probabilities
            MoveAcceptanceMethod = self.file_settings['optimization']['Move_Acceptance_Method']
            if MoveAcceptanceMethod == 0:
                MoveAcceptanceRatio = TotalMoves / (self.num_procs * self.file_settings['optimization']['population_size'])
            else:
                MoveAcceptanceRatio = TotalAcceptanceProbability / (self.num_procs * self.file_settings['optimization']['population_size'])
            print_dict['move_acceptance_ratio'] = MoveAcceptanceRatio
            # update Buffer
            Buffer, BufferCost = UpdateBuffer(Buffer, BufferCost, NewSolutions, NewSolutionsfitness)
            # update active solution
            active, activeCost = UpdateActive(self.cooling_schedule.temperature, Buffer, BufferCost)


            for i in range(len(BufferCost)):
                if BestSolutionCost > BufferCost[i]:
                    BestSolution = Buffer[i]
                    BestSolutionCost = BufferCost[i]

            tempdeviation = LAMStats(self.cooling_schedule.temperature, TSolutionsfitness)
            temp = LAM(self.cooling_schedule.temperature, tempdeviation, MoveAcceptanceRatio)
            self.cooling_schedule.temperature = temp
            print_dict['temperature_update']=self.cooling_schedule.temperature
            print_dict['standard_deviation']=tempdeviation

            opt.record_best_and_new_solution(BestSolution, active, self.cooling_schedule)
            end=time.time()
            PSAPrint(print_dict)
            BufferPrint(Buffer)
            f = open(version.__ofile__, 'a')
            f.write('\n\nPSA TEMPERATURE ITERATION {} COMPLETED IN {}s'.format(x+1,end-start))
            f.close()
        track_file = open('optimization_track_file.txt', 'a')
        track_file.write("End of Optimization \n")
        track_file.close()

        # plot the parameters over time
        opt.plotter()
        self.cleanup()

    def cleanup(self):
        """
        Deletes solution results that are no longer relevant to the optimization.
        Written by Brian Andersen 10/28/2020.
        """
        if self.perform_cleanup:
            if self.generation.current <= self.number_generations_post_cleanup:
                pass
            else:
                file_name = f'initial_solution'
                if os.path.isdir(file_name):
                    delete_file = True
                    for parent, child in zip(self.population.parents, self.population.children):
                        if parent.name == file_name:
                            delete_file = False
                            break
                        elif child.name == file_name:
                            delete_file = False
                            break
                    if delete_file:
                        os.rmdir(file_name)
                for clgen in range(self.generation.current):
                    for clpop in range(self.population.size):
                        file_name = f"solution_{clgen}_{clpop}"
                        if os.path.isdir(file_name):
                            delete_file = True
                            for parent, child in zip(self.population.parents, self.population.children):
                                if parent.name == file_name:
                                    delete_file = False
                                    break
                                elif child.name == file_name:
                                    delete_file = False
                                    break
                            if delete_file:
                                os.rmdir(file_name)