import java.io.IOException;
import java.util.Random;

import java.util.*;

/* A simple example of how the DemandPrediction class could be used as part
 * of a random search. It is not expected that you use this code as part of
 * your solution - it is just a demonstration of how the class's methods can be
 * called and how we can use two versions of the problem (here train and
 * test) to, respectively, obtain a promising set of parameters (using
 * train) and then to measure their performance (using test).
 */
public class Evo {
    public static void main(String[] args) throws IOException {
        var training_problem = new DemandPrediction("train");
        var bounds = DemandPrediction.bounds();
        var r = new Random();

        geneticAlgorithm3();
        /*
         * Generate N_TRIES random parameters and measure their MSE on the test
         * problem, saving the best parameters.
         */
        int N_TRIES = 100;
        double[] best_parameters = random_parameters(bounds, r);
        double best_training_error = training_problem.evaluate(best_parameters);
        for (int i = 0; i < N_TRIES - 1; i++) {
            var parameters = random_parameters(bounds, r);
            var training_error = training_problem.evaluate(parameters);
            if (training_error < best_training_error) {
                best_training_error = training_error;
                best_parameters = parameters;
            }
        }

    }

    public static double[] random_parameters(double[][] bounds, Random r) {
        var parameters = new double[bounds.length];
        for (int j = 0; j < bounds.length; j++) {
            parameters[j] = bounds[j][0] + r.nextDouble() * (bounds[j][1] - bounds[j][0]);
        }
        return parameters;
    }

    // get fitness of a solution by training error
    public static double getFitness(double[] parameters) throws IOException {
        var test_problem = new DemandPrediction("train");
        return test_problem.evaluate(parameters);
    }

    // get fitness of a solution by test error
    public static double getTestFitness(double[] parameters) throws IOException {
        var test_problem = new DemandPrediction("test");
        return test_problem.evaluate(parameters);
    }

    // get population of solutions
    public static double[][] getPopulation(int populationSize, double[][] bounds, Random r) {
        var population = new double[populationSize][bounds.length];
        for (int i = 0; i < populationSize; i++) {
            population[i] = random_parameters(bounds, r);
        }
        return population;
    }

    // normalize fitness
    public static double normalizeFitness(double fitness) {
        return 1 / (1 + fitness);
    }

    // tournament selection
    public static double[] tournamentSelection(double[][] population, int tournamentSize, Random r) {
        var tournament = new double[tournamentSize][population[0].length];
        for (int i = 0; i < tournamentSize; i++) {
            tournament[i] = population[r.nextInt(population.length)];
        }
        return tournament[0];
    }

    // crossover
    public static double[] crossover(double[] parent1, double[] parent2, Random r) {
        var child = new double[parent1.length];
        for (int i = 0; i < parent1.length; i++) {
            if (r.nextDouble() < 0.5) {
                child[i] = parent1[i];
            } else {
                child[i] = parent2[i];
            }
        }
        return child;
    }

    // breed via tournament selection
    public static double[][] breed(double[][] population, int tournamentSize, Random r) {
        var newPopulation = new double[population.length][population[0].length];
        for (int i = 0; i < population.length; i++) {
            var parent1 = tournamentSelection(population, tournamentSize, r);
            var parent2 = tournamentSelection(population, tournamentSize, r);
            newPopulation[i] = crossover(parent1, parent2, r);
        }
        return newPopulation;
    }

    // non tournament selection crossover and mutation for genetic algorithm
    public static double[][] breed2(double[][] population, Random r) {
        var newPopulation = new double[population.length][population[0].length];
        for (int i = 0; i < population.length; i++) {
            var parent1 = population[r.nextInt(population.length)];
            var parent2 = population[r.nextInt(population.length)];
            newPopulation[i] = crossover(parent1, parent2, r);
        }
        return newPopulation;
    }

    // every 150 generations, breed without tournament selection
    public static double[][] breed3(double[][] population, int tournamentSize, int generation, Random r) {
        var newPopulation = new double[population.length][population[0].length];
        if (generation % 150 == 0) {
            for (int i = 0; i < population.length; i++) {
                var parent1 = population[r.nextInt(population.length)];
                var parent2 = population[r.nextInt(population.length)];
                newPopulation[i] = crossover(parent1, parent2, r);
            }
        } else {
            for (int i = 0; i < population.length; i++) {
                var parent1 = tournamentSelection(population, tournamentSize, r);
                var parent2 = tournamentSelection(population, tournamentSize, r);
                newPopulation[i] = crossover(parent1, parent2, r);
            }
        }
        return newPopulation;
    }

    // mutate population for a lower fitness score 
    public static double[][] mutate(double[][] population, double[][] bounds, Random r) {
        var newPopulation = new double[population.length][population[0].length];
        for (int i = 0; i < population.length; i++) {
            for (int j = 0; j < population[i].length; j++) {
                if (r.nextDouble() < 0.01) {
                    newPopulation[i][j] = bounds[j][0] + r.nextDouble() * (bounds[j][1] - bounds[j][0]);
                } else {
                    newPopulation[i][j] = population[i][j];
                }
            }
        }
        return newPopulation;
    }

    // mutate the solutions to get for a better fitness score
    public static double[][] mutate2(double[][] population, double[][] bounds, Random r) {
        var Population = new double[population.length][bounds.length];
        
        for (int i = 0; i < population.length; i++) {
            var child = population[i];
            for (int j = 0; j < bounds.length; j++) {
                if (r.nextDouble() < 0.) {
                    child[j] = bounds[j][0] + r.nextDouble() * (bounds[j][1] - bounds[j][0]);
                }
            }
            Population[i] = child;
        }
        return Population;
    }

    // get the best solution from the population
    public static double[] getBest(double[][] population) throws IOException {
        var best = population[0];
        var bestFitness = getFitness(best);
        for (int i = 1; i < population.length; i++) {
            var fitness = getFitness(population[i]);
            if (fitness < bestFitness) {
                bestFitness = fitness;
                best = population[i];
            }
        }
        return best;
    }

    // get train error of best solution
    public static double getBestFitness(double[] best) throws IOException {
        var test_problem = new DemandPrediction("train");
        return test_problem.evaluate(best);
    }

    // run the genetic algorithm
    public static void geneticAlgorithm() throws IOException {
        var bounds = DemandPrediction.bounds();
        var r = new Random();
        var populationSize = 1000;
        var tournamentSize = 100;
        var population = getPopulation(populationSize, bounds, r);
        var best = getBest(population);
        var bestFitness = getBestFitness(best);
        for (int i = 0; i < 10000; i++) {
            population = breed(population, tournamentSize, r);
            population = mutate(population, bounds, r);
            var newBest = getBest(population);
            var newBestFitness = getBestFitness(newBest);
            if (newBestFitness < bestFitness) {
                bestFitness = newBestFitness;
                best = newBest;
                System.out.println("Best fitness: " + bestFitness);
            }
        }
        System.out.println("Best fitness: " + bestFitness);
        System.out.println("Best solution: " + Arrays.toString(best));
    }

    // run the genetic algorithm with non tournament selection
    public static void geneticAlgorithm2() throws IOException {
        var bounds = DemandPrediction.bounds();
        var r = new Random();
        var populationSize = 10000;
        var population = getPopulation(populationSize, bounds, r);
        var best = getBest(population);
        var bestFitness = getBestFitness(best);
        for (int i = 0; i < 100000; i++) {
            population = breed2(population, r);
            population = mutate(population, bounds, r);
            var newBest = getBest(population);
            var newBestFitness = getBestFitness(newBest);
            if (newBestFitness < bestFitness) {
                bestFitness = newBestFitness;
                best = newBest;
                System.out.println("Best fitness: " + bestFitness);
            }
        }
        System.out.println("Best fitness: " + bestFitness);
        System.out.println("Best solution: " + Arrays.toString(best));
    }

    // run the genetic algorithm with non tournament selection and mutation every 10
    // generations
    public static void geneticAlgorithm3() throws IOException {
        var bounds = DemandPrediction.bounds();
        var r = new Random();
        var populationSize = 50000;
        //float t= sqrt(populationSize)/2;
        var tournamentSize = 111 ;
        var population = getPopulation(populationSize, bounds, r);
        var best = getBest(population);
        var bestFitness = getBestFitness(best);
        for (int i = 0; i < 100; i++) {
            population = breed3(population, tournamentSize, i, r);
            population = mutate(population, bounds, r);
            var newBest = getBest(population);
            var newBestFitness = getBestFitness(newBest);
            System.out.println("generation: " + i);
            System.out.println("generational fitness " + newBestFitness);

            if (newBestFitness < bestFitness) {
                bestFitness = newBestFitness;
                best = newBest;
                System.out.println("Best fitness: " + bestFitness);
                System.out.println("Best solution: " + Arrays.toString(best));
                System.out.println("best test fitness: " + getTestFitness(best));
                System.out.println("generation: " + i);

            }
        }
        System.out.println("Best fitness: " + bestFitness);
        System.out.println("Best solution: " + Arrays.toString(best));
        System.out.println("best test fitness: " + getTestFitness(best));

    }
    
    
   

}
