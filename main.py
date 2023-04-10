import random

#Eight Queen Chromosome
#By Ella Mohanram
#April 7, 2023

#solves eight queens problem using genetic algorithms
class EightQueensChrom:
    #initializes chess board size
    board_size = 8

    #initializes number of chromosomes
    num_chroms = 7

    #initializes number of simulations
    num_sims = 1000

    #initializes number of parents
    num_parents = 2

    #intializes number of mutations
    mutate_count = 2

    #constructor function that intializes class attributes
    def __init__(self):
        self.my_chroms = []

    #calculates fitness of this chromosome (number of queens each queen can attack)
    #param: chromosome: 1D array holding list of numbers
    #return: calculated fitness value
    def get_fitness(self, chromosome):
        score = 0
        for current_index in range(len(chromosome)):
            for index in range(len(chromosome)):
                if chromosome[current_index] == chromosome[index] and current_index != index:
                    score += 1
                if abs(current_index - index) == abs(chromosome[current_index] - chromosome[index]) and current_index != index:
                    score += 1

        return score

    #creates a starting population of random chromosomes
    def general_initial_population(self):
        while len(self.my_chroms) != EightQueensChrom.num_chroms:
            my_chrom = []
            while len(my_chrom) != EightQueensChrom.board_size:
                my_chrom.append(random.randint(0,7))
            self.my_chroms.append(my_chrom)

    #creates the roulette wheel for selecting parents
    #return: array containing scaled selection probabilities based on fitness
    def set_probabilities_of_population(self):
        fitness_values = []
        scaled_probability = []
        for my_chrom in self.my_chroms:
            fitness_values.append(EightQueensChrom.get_fitness(self, my_chrom))

        sum_probs = sum(fitness_values)
        for fitness_value in fitness_values:
            scaled_probability.append(fitness_value/sum_probs)

        return scaled_probability

    #picks random parents to use for creating offspring
    #return: array containing selected parents
    def roulette_wheel_selection(self):
        population_probabilities = EightQueensChrom.set_probabilities_of_population(self)
        selected_parents = []

        while len(selected_parents) != 2:
            random_probability = random.random()
            i = 0
            while random_probability > population_probabilities[i]:
                random_probability = random_probability - population_probabilities[i]
                i += 1
            if population_probabilities[i] not in selected_parents:
                selected_parents.append(self.my_chroms[i])

        return selected_parents

    #produces children from the selected parents using two-point crossover
    #param: the_chosen: the selected parents from roulette wheel selection
    #return: the array of offspring
    def reproduce_children(self, the_chosen):
        the_offspring = []
        percent_amount = 0.25
        parent_one = the_chosen[0]
        parent_two = the_chosen[1]

        while len(the_offspring) < (EightQueensChrom.num_chroms * percent_amount):
            crossover_indices = sorted(random.sample(range(0, len(the_chosen[0])), 2))
            crossover_index_one = crossover_indices[0]
            crossover_index_two = crossover_indices[1]
            #child = parent_one[:(crossover_index_one+1)] + parent_two[crossover_index_one+1:]
            child = parent_one[:(crossover_index_one + 1)] + parent_two[(crossover_index_one + 1):crossover_index_two + 1] + parent_one[(crossover_index_two+1):]
            the_offspring.append(child)

        return the_offspring

    #mutates each offspring using single-point mutation (MUTATE_COUNT times)
    #param: the_children: children created from reproduce_children function
    #return: the array of mutated children
    def mutate_children(self, the_children):
        for child in the_children:
            mutated_indices = []
            for i in range(EightQueensChrom.mutate_count):
                mutated_indices.append(random.randint(0,7))
            for mutated_index in mutated_indices:
                mutation = random.randint(0,7)
                child[mutated_index] = mutation

        return the_children

    #sorts myChroms in order based on getFitness() values
    def sort_chroms(self):
        n = len(self.my_chroms)
        swapped = False

        for i in range(n-1):
            for j in range(n-i-1):
                if EightQueensChrom.get_fitness(self, self.my_chroms[j]) > EightQueensChrom.get_fitness(self, self.my_chroms[j+1]):
                    swapped = True
                    self.my_chroms[j], self.my_chroms[j + 1] = self.my_chroms[j + 1], self.my_chroms[j]

    #adds given children to the population by removing the weakest population members
    #param: the_children: array of new chromosomes to add to myChroms
    #return: sorted chromosomes (num_chrom chromsomes)
    def merge_population_and_children(self, the_children):
        EightQueensChrom.sort_chroms(self)
        sorted_chroms = self.my_chroms[:-(len(the_children))] + the_children

        return sorted_chroms

    #runs the simulation
    #return: the fittest chromosome
    def run_ga(self):
        EightQueensChrom.general_initial_population(self)
        best_global_fitness = EightQueensChrom.get_fitness(self, self.my_chroms[0])
        best_global_chromosome = self.my_chroms[0]
        i = 0
        while i <= EightQueensChrom.num_sims:
            if i % (EightQueensChrom.num_sims * 0.1) == 0:
                for k in self.my_chroms:
                    print(str(k) + ": " + str(EightQueensChrom.get_fitness(self, k)))
                input("Press Enter to continue...")
            for chromosome in self.my_chroms:
                if best_global_fitness > EightQueensChrom.get_fitness(self, chromosome):
                    best_global_fitness = EightQueensChrom.get_fitness(self, chromosome)
                    best_global_chromosome = chromosome
                    if best_global_fitness == 0:
                        i = EightQueensChrom.num_sims + 1
            the_chosen = EightQueensChrom.roulette_wheel_selection(self)
            the_children = EightQueensChrom.reproduce_children(self, the_chosen)
            the_children = EightQueensChrom.mutate_children(self, the_children)
            self.my_chroms = EightQueensChrom.merge_population_and_children(self, the_children)
            i += 1

        best_score = ("Best Chromosome: " + str(best_global_chromosome) + ", Score: " + str(best_global_fitness))
        print(best_score)

eightqueenchrom = EightQueensChrom()
eightqueenchrom.run_ga()
