import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import gui
import numpy as np
import matplotlib as mpl
import math
import random

class AminoAcidType:
    # H represents Hydrophobic amino acid.
    H = 'h'
    # P represents Hydrophilic amino acid.
    P = 'p'

class AngleType:
    angle_90 = 1
    angle_180 = 2
    angle_270 = 3


def get_chromosome_figure(genotype):
    G = nx.Graph()

    for node in range(genotype.chromosome_size):
        color = ''
        if genotype.phenotype[node] == AminoAcidType.H:
            color = 'black'
        elif genotype.phenotype[node] == AminoAcidType.P:
            color = 'red'
        G.add_node(node, pos=(genotype.x[node], genotype.y[node]), color=color)

        if node != 0:
            G.add_edge(node, node - 1)

    pos = nx.get_node_attributes(G, 'pos')

    colors = nx.get_node_attributes(G, 'color')

    fig = plt.figure(figsize=(4, 4))
    nx.draw(G, pos, with_labels=True, node_color=list(colors.values()))

    return fig

class Chromosome:
    def __init__(self, chromosome_size, input):
        self.chromosome_size = chromosome_size
        self.x = [0] * chromosome_size
        self.y = [0] * chromosome_size
        self.fitness = 0
        self.phenotype = input
        self.init_chain()

    def init_chain(self):
        self.x[0] = 0
        self.y[0] = 0

        self.x[1] = 1
        self.y[1] = 0


    def copy(self):
        g = Chromosome(self.chromosome_size, self.phenotype)
        g.x = list(self.x)
        g.y = list(self.y)

        return g

    def get_DataFrame(self):
        df = pd.DataFrame()
        df['x'] = self.x
        df['y'] = self.y

        return df

    def rotate(self, start_index, end_index, angle_type):
        """
            Rotate a point counterclockwise by a given angle around a given origin.

            angle = 1 means 90 degree, ie. pi /2
            angle = 2 means 180 degree, ie. 2 * pi /2
            angle = 3 means 270 degree, ie. 3 * pi /2

            The angle should be given in radians.
        """

        radian_angle = angle_type * math.pi / 2

        if start_index < self.chromosome_size and end_index < self.chromosome_size:



            dx = self.x[start_index] - 0
            dy = self.y[start_index] - 0

            # shift to origin
            for i in range(start_index, end_index + 1):
                self.x[i] = self.x[i]  -1 * dx
                self.y[i] = self.y[i] -1 * dy

            # rotate
            for i in range(start_index, end_index + 1):
                px, py = self.x[i], self.y[i]

                qx = math.cos(radian_angle) * (px) - math.sin(radian_angle) * (py)
                qy = math.sin(radian_angle) * (px) + math.cos(radian_angle) * (py)

                self.x[i] = qx
                self.y[i] = qy

            # move to original position
            for i in range(start_index, end_index + 1):
                self.x[i] = self.x[i] + dx
                self.y[i] = self.y[i] + dy

def random_orientation(chromosome):
    # select direction as follows:
    #                                             3
    #                  Select Direction as:     2 X 1
    #                                             4
    num_possible_direction = 3
    a = [0] * num_possible_direction
    ax = [0] * num_possible_direction
    ay = [0] * num_possible_direction
    valid_folding = True

    previous_direction = 1
    for i in range(2, chromosome.chromosome_size):
        if previous_direction == 1:
            a[0] = 1
            ax[0] = 1
            ay[0] = 0

            a[1] = 3
            ax[1] = 0
            ay[1] = 1

            a[2] = 4
            ax[2] = 0
            ay[2] = -1

        if previous_direction == 2:
            a[0] = 2
            ax[0] = -1
            ay[0] = 0

            a[1] = 3
            ax[1] = 0
            ay[1] = 1

            a[2] = 4
            ax[2] = 0
            ay[2] = -1

        if previous_direction == 3:
            a[0] = 1
            ax[0] = 1
            ay[0] = 0

            a[1] = 2
            ax[1] = -1
            ay[1] = 0

            a[2] = 3
            ax[2] = 0
            ay[2] = 1

        if previous_direction == 4:
            a[0] = 1
            ax[0] = 1
            ay[0] = 0

            a[1] = 2
            ax[1] = -1
            ay[1] = 0

            a[2] = 4
            ax[2] = 0
            ay[2] = -1

        choice_list = [0, 1, 2]

        index_selected = random.choice(choice_list)

        present_direction = a[index_selected]
        X = chromosome.x[i - 1] + ax[index_selected]
        Y = chromosome.y[i - 1] + ay[index_selected]

        if has_conflict(population=chromosome, position=i - 1, X=X, Y=Y):
            choice_list.remove(index_selected)

            index_selected = random.choice(choice_list)

            present_direction = a[index_selected]
            X = chromosome.x[i - 1] + ax[index_selected]
            Y = chromosome.y[i - 1] + ay[index_selected]

            if has_conflict(population=chromosome, position=i - 1, X=X, Y=Y):
                choice_list.remove(index_selected)

                index_selected = random.choice(choice_list)

                present_direction = a[index_selected]
                X = chromosome.x[i - 1] + ax[index_selected]
                Y = chromosome.y[i - 1] + ay[index_selected]

                if has_conflict(population=chromosome, position=i - 1, X=X, Y=Y):
                    valid_folding = False
                    print "invalid"
                    break

        if valid_folding:
            chromosome.x[i] = X
            chromosome.y[i] = Y
            previous_direction = present_direction

    return valid_folding


def has_conflict(population, position, X, Y):
    for j in range(0, position + 1):
        if X == population.x[j] and Y == population.y[j]:
            return True

    return False


def compute_full_fitness(chromosome):
    fitness = 0
    for i in range(chromosome.chromosome_size - 2):
        for j in range(i + 2, chromosome.chromosome_size):
            fitness += get_fitness(population=chromosome, i=i, j=j)

    return fitness

def get_fitness(population, i, j):
    # If both are hydrophobic, then we will consider to check the next condition otherwise return 0.
    if population.phenotype[i] == AminoAcidType.H and population.phenotype[j] == AminoAcidType.H:
        # if j is one of four neighbours of i then return -1.
        #   j
        # j i j
        #   j
        #
        if (population.x[i] + 1 == population.x[j] and population.y[i] == population.y[j]) or (population.x[i] - 1 == population.x[j] and population.y[i] == population.y[j]) or (population.x[i] == population.x[j] and population.y[i] + 1 == population.y[j]) or (population.x[i] == population.x[j] and population.y[i] - 1 == population.y[j]):
            return -1

    return 0

def mutation(chromosome, n):
    mutation_successful = False

    choice_list = [AngleType.angle_90, AngleType.angle_180, AngleType.angle_270]

    while len(choice_list) > 0:
        angle = random.choice(choice_list)
        choice_list.remove(angle)

        new_chromosome = chromosome.copy()
        new_chromosome.rotate(start_index = n, end_index = new_chromosome.chromosome_size - 1, angle_type = angle)

        if is_overlapped(new_chromosome) == False:
            return True, new_chromosome

    return mutation_successful, new_chromosome


# def mutation(chromosome, n):
#     mutation_successful = False
#
#     new_chromosome = chromosome.copy()
#
#     ary = [1,2,3]
#
#     a = chromosome.x[n]
#     b = chromosome.y[n]
#
#     collision = False
#
#     while len(ary) != 0:
#         p = random.choice(ary)
#         ary.remove(p)
#
#         for k in range(n+1, chromosome.chromosome_size):
#             if p == 1:
#                 new_chromosome.x[k] = a + b - chromosome.y[k] #  X = (a+b)-Y
#                 new_chromosome.y[k] = chromosome.x[k] + b - a # Y = (X+b)-a
#             elif p == 2:
#                 new_chromosome.x[k] = 2 * a - chromosome.x[k]               # X = (2a - X)
#                 new_chromosome.y[k] = 2 * b - chromosome.y[k]               # Y = (2b - Y)
#             elif p == 3:
#                 new_chromosome.x[k] = chromosome.y[k] + a - b # X =  Y+a-b
#                 new_chromosome.y[k] = a + b - chromosome.x[k] # Y =  (a+b)-X
#
#             for z in range(n + 1):
#                 if new_chromosome.x[k] == chromosome.x[z] and new_chromosome.y[k] == chromosome.y[z]:
#                     collision = True
#                     break
#
#             if collision == True:
#                 break
#         if collision == False:
#             break
#
#         collision = False
#
#     if collision == False:
#         mutation_successful = True
#
#     return mutation_successful, new_chromosome

# cross_over on population_i and population_h, n = cross point
# first part of population i and 2nd part of poulation j.
def cross_over(chromosome_i, chromosome_j, n):
    new_chromosome_j = chromosome_j.copy()

    prev_direction = -1

    ax = [0] * 3
    ay = [0] * 3

    if chromosome_i.x[n] == chromosome_i.x[n-1]:
        p = chromosome_i.y[n - 1] - chromosome_i.y[n]

        if p == 1:
            prev_direction = 3
        else:
            prev_direction = 4

    else:
        p = chromosome_i.x[n - 1] - chromosome_i.x[n]

        if p == 1:
            prev_direction = 1
        else:
            prev_direction = 2

    if prev_direction == 1:
        ax[0] = -1
        ay[0] = 0

        ax[1] = 0
        ay[1] = 1

        ax[2] = 0
        ay[2] = -1
    elif prev_direction == 2:
        ax[0] = 1
        ay[0] = 0

        ax[1] = 0
        ay[1] = 1

        ax[2] = 0
        ay[2] = -1
    elif prev_direction == 3:
        ax[0] = 1
        ay[0] = 0

        ax[1] = -1
        ay[1] = 0

        ax[2] = 0
        ay[2] = -1
    elif prev_direction == 4:
        ax[0] = 1
        ay[0] = 0

        ax[1] = -1
        ay[1] = 0

        ax[2] = 0
        ay[2] = 1

    choice_list = [0, 1, 2]

    is_cross_over_i_j_successfull = False

    while (len(choice_list) > 0):
        index_selected = random.choice(choice_list)
        choice_list.remove(index_selected)

        new_chromosome_j.x[n + 1] = chromosome_i.x[n] + ax[index_selected]
        new_chromosome_j.y[n + 1] = chromosome_i.y[n] + ay[index_selected]

        dx = new_chromosome_j.x[n + 1] - chromosome_j.x[n + 1]
        dy = new_chromosome_j.y[n + 1] - chromosome_j.y[n + 1]

        is_valid_shift, new_chromosome_j = get_valid_shift(population_i = chromosome_i,
                                                           population_j = chromosome_j,
                                                           start_index = n + 1,
                                                           end_index =chromosome_i.chromosome_size - 1,
                                                           dx = dx,
                                                           dy = dy)

        if is_valid_shift:
            for k in range(n + 1):
                new_chromosome_j.x[k] = chromosome_i.x[k]
                new_chromosome_j.y[k] = chromosome_i.y[k]
            return True, new_chromosome_j

        # if is_valid_shift == False:
        angle_choice_list = [AngleType.angle_90, AngleType.angle_180, AngleType.angle_270]

        while len(angle_choice_list) > 0:
            angle_selected = random.choice(angle_choice_list)

            angle_choice_list.remove(angle_selected)

            chromosome_j_copy = chromosome_j.copy()

            chromosome_j_copy.rotate(start_index = n + 1, end_index  = new_chromosome_j.chromosome_size - 1, angle_type = angle_selected)

            is_valid_shift, new_chromosome_j = get_valid_shift(population_i=chromosome_i,
                                                               population_j=chromosome_j_copy,
                                                               start_index=n + 1,
                                                               end_index=chromosome_i.chromosome_size - 1,
                                                               dx=dx,
                                                               dy=dy)

            if is_valid_shift:
                for k in range(n + 1):
                    new_chromosome_j.x[k] = chromosome_i.x[k]
                    new_chromosome_j.y[k] = chromosome_i.y[k]
                return True, new_chromosome_j


    # if is_overlapped(new_chromosome_j):
    return False, new_chromosome_j

def get_valid_shift(population_i, population_j, start_index, end_index, dx ,dy):
    new_population_j = population_j.copy()

    is_shift_valid = True
    for k in range(start_index, end_index + 1):
        new_population_j.x[k] = population_j.x[k] + dx
        new_population_j.y[k] = population_j.y[k] + dy

        if has_conflict(population = population_i, position = start_index - 1, X = new_population_j.x[k], Y = new_population_j.y[k]):
            is_shift_valid = False
            return is_shift_valid, new_population_j


    return is_shift_valid, new_population_j

def initialization(population_size, input):
    population1 = []

    for i in range(population_size):
        chromosome = Chromosome(chromosome_size= len(input), input=input)
        valid_folding = random_orientation(chromosome)

        while is_overlapped(chromosome):
            valid_folding = random_orientation(chromosome)

        chromosome.fitness = compute_full_fitness(chromosome)
        population1.append(chromosome)

    return population1

# compute fitness for population
def compute_full_fitness_for_population(population):
    for i in range(len(population)):
        population[i].fitness = compute_full_fitness(population[i])


def run_genetic_algorithm(
        input,
        max_fitness,
        population_size,
        elite_rate,
        cross_over_rate,
        mutation_rate,
        max_gen
):
    population1 = initialization(population_size, input)

    gen = 0

    while gen < max_gen:
        population2 = []
        elite_size = int(population_size * elite_rate / 100.0)
        compute_full_fitness_for_population(population1)

        population1_sorted = sorted(population1, key=lambda chromosome: chromosome.fitness)

        if abs(population1_sorted[0].fitness) >= abs(max_fitness):
            print "maximum fitness found"
            print population1_sorted[0].fitness
            print population1_sorted[0].get_DataFrame()
            return population1_sorted[0]


        # elite
        population2 = population1_sorted[:elite_size]

        # cross over
        cross_over_size = int(population_size * cross_over_rate / 100.0)

        for i in range(cross_over_size/2):
            index1 = roulette_wheel_selection(population1_sorted)
            index2 = roulette_wheel_selection(population1_sorted)

            cross_over_index = random.randint(2, population1_sorted[i].chromosome_size - 2)
            is_cross_over_i_j_successfull, new_chromosome_j = cross_over(chromosome_i=population1_sorted[index1],
                                                                         chromosome_j=population1_sorted[index2],
                                                                         n=cross_over_index)

            is_cross_over_j_i_successfull, new_chromosome_i = cross_over(chromosome_i=population1_sorted[index2],
                                                                         chromosome_j=population1_sorted[index1],
                                                                         n=cross_over_index)

            while is_cross_over_i_j_successfull == False or is_cross_over_j_i_successfull == False:
                index1 = roulette_wheel_selection(population1_sorted)
                index2 = roulette_wheel_selection(population1_sorted)

                cross_over_index = random.randint(2, population1_sorted[i].chromosome_size - 2)
                is_cross_over_i_j_successfull, new_chromosome_j = cross_over(chromosome_i=population1_sorted[index1],
                                                                             chromosome_j=population1_sorted[index2],
                                                                             n= cross_over_index)

                is_cross_over_j_i_successfull, new_chromosome_i = cross_over(chromosome_i=population1_sorted[index2],
                                                                             chromosome_j=population1_sorted[index1],
                                                                             n=cross_over_index)

            if is_cross_over_i_j_successfull and is_cross_over_j_i_successfull:
                if is_overlapped(new_chromosome_j) or is_overlapped(new_chromosome_i):
                    print "!!! overlapped for cross"

                population2.append(new_chromosome_j)
                population2.append(new_chromosome_i)


        for i in range(len(population2), len(population1_sorted)):
            population2.append(population1_sorted[i])

        # mutation
        mutation_size = int(mutation_rate * len(population2)/ 100.0)

        for i in range(mutation_size):
            index = random.randint(elite_size + 1, len(population2) - 1)

            print "index: " + str(index)
            print "pop2: " + str(len(population2))
            mutation_index = random.randint(2, population2[index].chromosome_size - 2)
            mutation_successful, new_chromosome = mutation(chromosome=population2[index], n = mutation_index)

            while mutation_successful != True:
                index = random.randint(elite_size + 1, len(population1_sorted))
                mutation_index = random.randint(2, population2[index].chromosome_size - 2)
                mutation_successful, new_chromosome = mutation(chromosome=population2[index],
                                                               n=mutation_index)

            if is_overlapped(new_chromosome):
                print "!!! overlapped for mutation"

            if mutation_successful:
                print "mutation"

            population2[index2] = new_chromosome


        population1 = population2

        gen = gen + 1

def roulette_wheel_selection(population):
    fitness_sum = 0

    for i in range(len(population)):
        fitness_sum += abs(population[i].fitness)

    random_int_fitness_sum = random.randint(0, fitness_sum)

    cumulative_fitness_deduction = random_int_fitness_sum

    for i in range(len(population)):
        cumulative_fitness_deduction -= abs(population[i].fitness)

        if cumulative_fitness_deduction <= 0:
            return i


def is_overlapped(chromosome):
    for i in range(chromosome.chromosome_size - 1):
        for j in range(i + 1, chromosome.chromosome_size):
            if chromosome.x[i] == chromosome.x[j] and chromosome.y[i] == chromosome.y[j]:
                return True

    return False


input = 'hphpphhphpphphhpphph'
best_chromosome = run_genetic_algorithm(
    input=input,
    max_fitness=-9,
    population_size=200,
    elite_rate = 5,
    cross_over_rate = 80,
    mutation_rate = 5,
    max_gen = 1000
)
fig1 = get_chromosome_figure(best_chromosome)
plt.show()

# input = 'hphpphhphpphphhp'
# chromosome_size = len(input)
#
# initialization(3, input)


#
# g = Chromosome(chromosome_size, input)
# valid_folding = random_orientation(chromosome=g)
# print("fitness: " + str(compute_full_fitness(population=g)))
# fig1 = get_chromosome_figure(g)
# print(g.get_DataFrame())
#
# # another population
# g1 = Chromosome(chromosome_size, input)
# valid_folding = random_orientation(chromosome=g1)
# print("fitness: " + str(compute_full_fitness(population=g1)))
# fig2 = get_chromosome_figure(g1)
# print(g1.get_DataFrame())
#
# # mutation_successful, new_g = mutation(population=g, n=1)
# # print(new_g.get_DataFrame())
# #
# #
# # print "mutation: " + str(mutation_successful)
# # fig2 = get_chromosome_figure(new_g)
#
# # g.rotate(start_index=0, end_index=g.chromosome_size - 1, angle_type=AngleType.angle_180)
# # fig2 = get_chromosome_figure(g)
# # print(g.get_DataFrame())
#
#
# is_success, new_pop_1 = cross_over(g, g1, 3)
# print "Is cross over success: " + str(is_success)
# fig3 = get_chromosome_figure(new_pop_1)
#
#
#
# is_success, new_pop_1 = cross_over(g1, g, 3)
# print "Is cross over success: " + str(is_success)
# fig4 = get_chromosome_figure(new_pop_1)
# plt.show()
#
#
#
#
#
#
#
#
# # Generate some example data
# X = np.linspace(0, 2 * np.pi, 50)
# Y = np.sin(X)
#
# # Create the figure we desire to add to an existing canvas
# fig = mpl.figure.Figure(figsize=(8, 8))
# ax = fig.add_axes([0, 0, 1, 1])
# ax.plot(X, Y)
#
# # Keep this handle alive, or else figure will disappear
# gui.init_draw(fig1, fig2)