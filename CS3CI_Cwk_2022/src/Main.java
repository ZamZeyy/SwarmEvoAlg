import java.io.IOException;
import java.lang.reflect.MalformedParameterizedTypeException;
import java.util.Random;

import java.util.*;

/* A simple example of how the DemandPrediction class could be used as part
 * of a random search. It is not expected that you use this code as part of
 * your solution - it is just a demonstration of how the class's methods can be
 * called and how we can use two versions of the problem (here train and
 * test) to, respectively, obtain a promising set of parameters (using
 * train) and then to measure their performance (using test).
 */
public class Main {
    public static void main(String[] args) throws IOException {
        var training_problem = new DemandPrediction("train");
        var bounds = DemandPrediction.bounds();
        var r = new Random();

        double[] PParameters3 = runPSO2(training_problem, bounds, 100000, 100,  0.37,  0.75,  0.75); //0.37,  0.75,  0.75
        //double[] PParameters3 = runScout(training_problem, bounds, 10000, 1000);
        var test = new DemandPrediction ("test");
        System.out.println("Test error: " + test.evaluate(PParameters3));
       //System.out.println("Test error2: " + test.evaluate(PParameters2));
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
    //as particles get closer to the target, their velocity decreases



    public static double[] random_parameters(double[][] bounds, Random r) {
        var parameters = new double[bounds.length];
        for (int j = 0; j < bounds.length; j++) {
            parameters[j] = bounds[j][0] + r.nextDouble() * (bounds[j][1] - bounds[j][0]);
        }
        return parameters;
    }

    //initialise the swarm with random values
    public static double[][] initialiseSwarm(double[][] bounds, Random r, int swarmSize) {
        var swarm = new double[swarmSize][bounds.length];
        for (int i = 0; i < swarmSize; i++) {
            for (int j = 0; j < bounds.length; j++) {
                swarm[i][j] = bounds[j][0] + r.nextDouble() * (bounds[j][1] - bounds[j][0]);
            }
        }
        return swarm;
    }
  
    //values cannot be outside the bounds
    public static double[][] checkBounds(double[][] swarm, double[][] bounds) {
        for (int i = 0; i < swarm.length; i++) {
            for (int j = 0; j < swarm[i].length; j++) {
                if (swarm[i][j] < bounds[j][0]) {
                    swarm[i][j] = bounds[j][0];
                }
                if (swarm[i][j] > bounds[j][1]) {
                    swarm[i][j] = bounds[j][1];
                }
            }
        }
        return swarm;
    }


    //update the velocity of the particles
    public static double[][] updateVelocity(double[][] swarm, double[][] velocity, double[][] pBest, double[] gBest, double w, double c1, double c2, Random r) {
        for (int i = 0; i < swarm.length; i++) {
            for (int j = 0; j < swarm[i].length; j++) {
                velocity[i][j] = w * velocity[i][j] + c1 * r.nextDouble() * (pBest[i][j] - swarm[i][j]) + c2 * r.nextDouble() * (gBest[j] - swarm[i][j]);
            }
        }
        return velocity;
    }

    //update the position of the particles and check the bounds
    public static double[][] updatePosition(double[][] swarm, double[][] velocity, double[][] bounds) {
        for (int i = 0; i < swarm.length; i++) {
            for (int j = 0; j < swarm[i].length; j++) {
                swarm[i][j] = swarm[i][j] + velocity[i][j];
            }
            
        }
        return checkBounds(swarm, bounds);
    }



    //update the position of the particles
    public static double[][] updatePosition2(double[][] swarm, double[][] velocity) {
        for (int i = 0; i < swarm.length; i++) {
            for (int j = 0; j < swarm[i].length; j++) {
                swarm[i][j] = swarm[i][j] + velocity[i][j];
            }
            
        }
        return swarm;
    }

    //update the personal best of the particles
    public static double[][] updatePBest(double[][] swarm, double[][] pBest, double[] fitness, DemandPrediction training_problem) {
        for (int i = 0; i < swarm.length; i++) {
            
            if (fitness[i] < training_problem.evaluate(pBest[i])) {
                
                pBest[i] = swarm[i];
            }
        }
        return pBest;
    }



    
 //iterate through the swarm and calculate the fitness of each particle
    public static double[] calculateFitness(double[][] swarm, DemandPrediction training_problem) {
        double[] fitness = new double[swarm.length];
        for (int i = 0; i < swarm.length; i++) {
            fitness[i] = training_problem.evaluate(swarm[i]);
        }
        return fitness;
    }


    //update the global best of the particles
    public static double[] updateGBest(double[][] swarm, double[] gBest, double[] fitness,DemandPrediction training_problem) {
       calculateFitness(swarm, training_problem);
        float p;
        float g;
        g = (float) training_problem.evaluate(gBest);
        for (int i = 0; i < swarm.length; i++) {
            //if fitness is less than the current gbest, update gbest
            p = (float) fitness[i];
            



            if (p < g) {
                
                g = p;

                 
                
                gBest = swarm[i];

                //System.out.println("New global best found: " + Arrays.toString(gBest));
                System.out.println("New global best found: " + g);
            }
        }
        return gBest;
    }

   

    //get the index of the particle with the best fitness
    public static int getBestParticle(double[] fitness) {
        int bestParticle = 0;
        for (int i = 0; i < fitness.length; i++) {
            if (fitness[i] < fitness[bestParticle]) {
                bestParticle = i;
            }
        }
        return bestParticle;
    }
    //compare gbest to test data
    public static void compareGBest(double[] gBest, DemandPrediction test_problem) {
        System.out.println("The global best found is: " + Arrays.toString(gBest));
        System.out.println("The global best found has a fitness of: " + test_problem.evaluate(gBest));
    }

    
    //run the PSO algorithm and test against test data
    public static double[] runPSO(DemandPrediction training_problem, double[][] bounds, int swarmSize, int maxIterations, double w, double c1, double c2) {
        var r = new Random();
        double[][] swarm = initialiseSwarm(bounds, r, swarmSize);
        double[][] velocity = initialiseSwarm(bounds, r, swarmSize);
        double[][] pBest = initialiseSwarm(bounds, r, swarmSize);
        double[] gBest = swarm[getBestParticle(calculateFitness(swarm, training_problem))];
        double[] fitness = calculateFitness(swarm, training_problem);
        for (int i = 0; i < maxIterations; i++) {
            velocity = updateVelocity(swarm, velocity, pBest, gBest, w, c1, c2, r);
            swarm = updatePosition(swarm, velocity, bounds);
            pBest = updatePBest(swarm, pBest, fitness, training_problem);
            gBest = updateGBest(swarm, gBest, fitness, training_problem);
            fitness = calculateFitness(swarm, training_problem);

        }
        System.out.println("The best parameters are: " + Arrays.toString(gBest));
        System.out.println("The best fitness is: " + fitness[getBestParticle(fitness)]);

        return gBest;
    }

    //run the PSO algorithm and test against test data, decreasing the inertia weight as the algorithm gets closer to the target
    public static double[] runPSO2(DemandPrediction training_problem, double[][] bounds, int swarmSize, int maxIterations, double w, double c1, double c2) {
        var r = new Random();
        double[][] swarm = initialiseSwarm(bounds, r, swarmSize);
        double[][] velocity = initialiseSwarm(bounds, r, swarmSize);
        double[][] pBest = initialiseSwarm(bounds, r, swarmSize);
        double[] gBest = swarm[getBestParticle(calculateFitness(swarm, training_problem))];
        double[] fitness = calculateFitness(swarm, training_problem);
        for (int i = 0; i < maxIterations; i++) {
            velocity = updateVelocity(swarm, velocity, pBest, gBest, w, c1, c2, r);
            swarm = updatePosition(swarm, velocity, bounds);
            fitness = calculateFitness(swarm, training_problem);
            pBest = updatePBest(swarm, pBest, fitness, training_problem);
            gBest = updateGBest(swarm, gBest, fitness, training_problem);
            
            w = (training_problem.evaluate(gBest))/1000;
            //w = w - w/maxIterations;
            //c1 = c1 - c1/maxIterations;
            //c2 = c2 - c2/maxIterations;
            
           
            
        }
                
            
           
            

            
            
        
        System.out.println("The best parameters are2: " + Arrays.toString(gBest));
        System.out.println("The best fitness is2: " + fitness[getBestParticle(fitness)]);

        return gBest;

        
    
    }
    //if pso converges to a local minima, the scout swarm will find a new global best
    
   
}