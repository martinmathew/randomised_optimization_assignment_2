
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Random;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

import dist.DiscreteDependencyTree;
import dist.DiscretePermutationDistribution;
import dist.DiscreteUniformDistribution;
import dist.Distribution;

import opt.DiscreteChangeOneNeighbor;
import opt.EvaluationFunction;
import opt.GenericHillClimbingProblem;
import opt.HillClimbingProblem;
import opt.NeighborFunction;
import opt.RandomizedHillClimbing;
import opt.SimulatedAnnealing;
import opt.SwapNeighbor;
import opt.example.*;
import opt.ga.CrossoverFunction;
import opt.ga.DiscreteChangeOneMutation;
import opt.ga.SingleCrossOver;
import opt.ga.GenericGeneticAlgorithmProblem;
import opt.ga.GeneticAlgorithmProblem;
import opt.ga.MutationFunction;
import opt.ga.NQueensFitnessFunction;
import opt.ga.StandardGeneticAlgorithm;
import opt.ga.SwapMutation;
import opt.prob.GenericProbabilisticOptimizationProblem;
import opt.prob.MIMIC;
import opt.prob.ProbabilisticOptimizationProblem;
import shared.FixedIterationTrainer;



public class RandomizedOPT {
	
	private static DecimalFormat df = new DecimalFormat("0.000");


	/** The n value **/
    private static final int NNN = 60;
    /** The t value */
    private static final int TTT = NNN / 10;
    /** The n value */
    private static final int NNNN = 10;
	
    /** The n value */
    private static final int N = 50;
    /** The n value */
    private static final int NN = 200;
    /** The t value */
    private static final int TT = NN / 5;
    public static void main(String[] args) throws IOException {
    	tsp();
    	rhc_tsp();
    	simulated_annealing_tsp();
    	ga_tsp();
    	mimic_tsp();
    
    // Four Peak Problem
    four_peaks();
    	rhc_four_peaks();
    	simulated_annealing_fourpeak();
    	ga_fourpeak();
    	mimic_four_peak();
    	
    	
    	//Continous Peaks
    	rhc_nqueens_prob();
    	simulated_annealing_nqueens();
    	ga_n_queens();
    	mimic_n_queens();
    	
    	
    		
    }
    
    
    
    public static void mimic_n_queens() throws IOException{

    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};
    	
    	int samples[] = {10,100, 200, 300, 400, 500, 600};


        	

        	
        	String listRHCeval[][] = new String[samples.length][numOfRestarts.length];
        	for(int z = 0;z<samples.length;z++) {
        		int keep = samples[z];
        	for(int i=0;i < numOfRestarts.length;i++)	
        	{
        		int iterCount = numOfRestarts[i];
        		System.out.println("Start iteration" + iterCount);
        		double evalMIMIC = 0;

        		
        		double runTimeForMIMIC = 0;
      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
          			double start, end;
        	        int[] ranges = new int[NNNN];
        	        Random random = new Random(NNNN);
        	        for (int ii = 0; ii < NNNN; ii++) {
        	        	ranges[ii] = random.nextInt();
        	        }
        	        NQueensFitnessFunction ef = new NQueensFitnessFunction();
        	        Distribution odd = new DiscretePermutationDistribution(N);
        	        NeighborFunction nf = new SwapNeighbor();
        	        MutationFunction mf = new SwapMutation();
        	        CrossoverFunction cf = new SingleCrossOver();
        	        Distribution df = new DiscreteDependencyTree(.1); 
        	        HillClimbingProblem hcp = new GenericHillClimbingProblem(ef, odd, nf);
        	        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(ef, odd, mf, cf);
        	        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(ef, odd, df);	
    		        
    		        
    		        //Mimic
    		        start = System.nanoTime();
    		        MIMIC mimic = new MIMIC(200, 20, pop);
    		        FixedIterationTrainer fit = new FixedIterationTrainer(mimic, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10, 9);
    		        runTimeForMIMIC = runTimeForMIMIC +  time;
    		        evalMIMIC = evalMIMIC + ef.value(mimic.getOptimal());
    		        		        

    	        
        		}
        		//Averaging the total optimal fitness values and total time for test runs
        		evalMIMIC /= runForAvg;
        		listRHCeval[z][i] = df.format(evalMIMIC);
        		runTimeForMIMIC /= runForAvg;
        	}
       
        }
        	 
            String fileName = "output/n_queens_result_mimic_keep_count.csv";
            writeToFile(fileName, numOfRestarts, listRHCeval);
    }
    
    
    
    public static void ga_n_queens() throws IOException {


    	int cross_overs_n_mutate[] = {20,40,60,80,100,120,140,160,180, 200};

        	
        	String listRHCeval[][] = new String[cross_overs_n_mutate.length][cross_overs_n_mutate.length];

        	int iterCount = 2000;
        	for(int i=0;i < cross_overs_n_mutate.length;i++)	
        	{
        		
        		double evalGA = 0;
        		int cross_over = cross_overs_n_mutate[i];
        		
        		double runTimeForGA = 0;
        		System.out.println("Cross Over " + cross_over);
        		for(int k=0;k< cross_overs_n_mutate.length;k++) {
        		

        		int mutate = cross_overs_n_mutate[k];
        		
        		System.out.println("Mutate " + mutate);

      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
          			double start, end;
        	        int[] ranges = new int[NNNN];
        	        Random random = new Random(NNNN);
        	        for (int ii = 0; ii < NNNN; ii++) {
        	        	ranges[ii] = random.nextInt();
        	        }
        	        NQueensFitnessFunction ef = new NQueensFitnessFunction();
        	        Distribution odd = new DiscretePermutationDistribution(N);
        	        NeighborFunction nf = new SwapNeighbor();
        	        MutationFunction mf = new SwapMutation();
        	        CrossoverFunction cf = new SingleCrossOver();
        	        Distribution df = new DiscreteDependencyTree(.1); 
        	        HillClimbingProblem hcp = new GenericHillClimbingProblem(ef, odd, nf);
        	        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(ef, odd, mf, cf);
        	        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(ef, odd, df);	
    		        
    		        
    		        //Genetic Algorithm
    		        
    		        //Genetic Algorithm
    		        start = System.nanoTime();
    		        StandardGeneticAlgorithm ga = new StandardGeneticAlgorithm(200, cross_over, mutate, gap);
    		        FixedIterationTrainer fit = new FixedIterationTrainer(ga, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10,9);
    		        runTimeForGA = runTimeForGA +  time;
    		        evalGA = evalGA + ef.value(ga.getOptimal());
    		        		        

    	        
        		}
        		
        		//Averaging the total optimal fitness values and total time for test runs
        		runTimeForGA /= runForAvg;

        		
        		evalGA /= runForAvg;
        		listRHCeval[i][k] = df.format(evalGA);
        		}
        		

        }
        	
            
            String fileName = "output/n_queens_result_ga_cross_over_mutate.csv";
            writeToFile(fileName, cross_overs_n_mutate, listRHCeval);	
    
    } 
    
    
    
    
    public static void simulated_annealing_nqueens() throws IOException{

    	double cooling_exps[] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000, 5000};

    	double start, end;

        	
        	String listRHCeval[][] = new String[cooling_exps.length][numOfRestarts.length];

        	
        	for(int i=0;i < cooling_exps.length;i++)
        		
        	{
        		double cooling_exp = cooling_exps[i];
        		for(int h=0;h < numOfRestarts.length;h++)	
            	{
        			int iterCount = numOfRestarts[h];
        		
        		System.out.println("Cooling Exponent " + cooling_exp);
        		double evalSA = 0;

        		
        		double runTimeForsa = 0;
      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
        			System.out.println("Current run" + currun);
        			
        	        int[] ranges = new int[NNNN];
        	        Random random = new Random(NNNN);
        	        for (int ii = 0; ii < NNNN; ii++) {
        	        	ranges[ii] = random.nextInt();
        	        }
        	        NQueensFitnessFunction ef = new NQueensFitnessFunction();
        	        Distribution odd = new DiscretePermutationDistribution(N);
        	        NeighborFunction nf = new SwapNeighbor();
        	        MutationFunction mf = new SwapMutation();
        	        CrossoverFunction cf = new SingleCrossOver();
        	        Distribution df = new DiscreteDependencyTree(.1); 
        	        HillClimbingProblem hcp = new GenericHillClimbingProblem(ef, odd, nf);
        	        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(ef, odd, mf, cf);
        	        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(ef, odd, df);	
    		        
    		        
    		        //Simulated Annealing
    		        start = System.nanoTime();
    		        SimulatedAnnealing sa = new SimulatedAnnealing(1E11, cooling_exp, hcp);
    		        FixedIterationTrainer fit = new FixedIterationTrainer(sa, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10,9);
    		        runTimeForsa = runTimeForsa +  time;
    		        evalSA = evalSA + ef.value(sa.getOptimal());
    		        		        
        		}
        		
        		//Averaging the total optimal fitness values and total time for test runs
        		runTimeForsa /= runForAvg;
        		evalSA /= runForAvg;

        		
        		listRHCeval[i][h] = df.format(evalSA);
              }
        	}
            
            String fileName = "output/n_queens_result_sa_coolin_fitness.csv";
            writeToFile(fileName, numOfRestarts, listRHCeval);

    } 
    
    
    
    
    
    
    
    
    
    
    
    
    public static void rhc_nqueens_prob() throws IOException{
    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};

        	String listRHCeval[] = new String[numOfRestarts.length];
        	double start, end;
        
        	for(int i=0;i < numOfRestarts.length;i++)
        		
        	{
        		int iterCount = numOfRestarts[i];
        		System.out.println("Start iteration" + iterCount);
        		double evalRhc = 0;

        		
        		double runTimeForrhc = 0;
      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
        			System.out.println("Current run" + currun);
        			
        	        int[] ranges = new int[NNNN];
        	        Random random = new Random(NNNN);
        	        for (int ii = 0; ii < NNNN; ii++) {
        	        	ranges[ii] = random.nextInt();
        	        }
        	        NQueensFitnessFunction ef = new NQueensFitnessFunction();
        	        Distribution odd = new DiscretePermutationDistribution(N);
        	        NeighborFunction nf = new SwapNeighbor();
        	        MutationFunction mf = new SwapMutation();
        	        CrossoverFunction cf = new SingleCrossOver();
        	        Distribution df = new DiscreteDependencyTree(.1); 
        	        HillClimbingProblem hcp = new GenericHillClimbingProblem(ef, odd, nf);
        	        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(ef, odd, mf, cf);
        	        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(ef, odd, df);	        
        	        
    		        //Randomized hill climbing
    		        start = System.nanoTime();
    		        RandomizedHillClimbing rhc = new RandomizedHillClimbing(hcp);      
    		        FixedIterationTrainer fit = new FixedIterationTrainer(rhc, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10,9);
    		        runTimeForrhc = runTimeForrhc +  time;
    		        evalRhc = evalRhc + ef.value(rhc.getOptimal());
    		        		  
        		}
        		
        		//Averaging the total optimal fitness values and total time for test runs
        		evalRhc /= runForAvg;

        		
        		runTimeForrhc /= runForAvg;

        		listRHCeval[i] = df.format(evalRhc);


        }
            
            String fileName = "output/nquuen_result_rhc_fitness.csv";
            writeToFile(fileName, numOfRestarts, listRHCeval);

                   
        }
        
        
    
    
    
    
    
    //Four Peaks
    public static void four_peaks() throws IOException {
        
    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};


    	
    	String rhcRunTime[] = new String[numOfRestarts.length];
    	String listSaTime[] = new String[numOfRestarts.length];
    	String listGATime[] = new String[numOfRestarts.length];
    	String listMimicTime[] = new String[numOfRestarts.length];
    	
    	String listRHCeval[] = new String[numOfRestarts.length];
    	String listSAeval[] = new String[numOfRestarts.length];
    	String listGAeval[] = new String[numOfRestarts.length];
    	String listMIMICeval[] = new String[numOfRestarts.length];
    	

    	for(int i=0;i < numOfRestarts.length;i++)
    		
    	{
    		int iterCount = numOfRestarts[i];
    		System.out.println("Start iteration" + iterCount);
    		double evalRhc = 0;
    		double evalSA = 0;
    		double evalGA = 0;
    		double evalMIMIC = 0;
    		
    		double runTimeForrhc = 0;
    		double runTimeForsa = 0;
    		double runTimeForGA = 0;
    		double runTimeForMIMIC = 0;
    		double time = 0;
        	int runForAvg = 5;
    		for (int currun= 0; currun < runForAvg; currun++)
    		{
    	    	double start, end;
    			System.out.println("Current run" + currun);
		    	int[] ranges = new int[NN];
		        Arrays.fill(ranges, 2);
		        
		        EvaluationFunction ef = new FourPeaksEvaluationFunction(TT);
		        Distribution odd = new DiscreteUniformDistribution(ranges);
		        NeighborFunction nf = new DiscreteChangeOneNeighbor(ranges);
		        MutationFunction mf = new DiscreteChangeOneMutation(ranges);
		        CrossoverFunction cf = new SingleCrossOver();
		        Distribution df = new DiscreteDependencyTree(.1, ranges); 
		        HillClimbingProblem hcp = new GenericHillClimbingProblem(ef, odd, nf);
		        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(ef, odd, mf, cf);
		        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(ef, odd, df);
		        
		        
		        //Randomized hill climbing
		        start = System.nanoTime();
		        RandomizedHillClimbing rhc = new RandomizedHillClimbing(hcp);      
		        FixedIterationTrainer fit = new FixedIterationTrainer(rhc, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10,9);
		        runTimeForrhc = runTimeForrhc +  time;
		        evalRhc = evalRhc + ef.value(rhc.getOptimal());
		        		        
		        
		        //Simulated Annealing
		        start = System.nanoTime();
		        SimulatedAnnealing sa = new SimulatedAnnealing(1E11, 0.5, hcp);
		        fit = new FixedIterationTrainer(sa, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10,9);
		        runTimeForsa = runTimeForsa +  time;
		        evalSA = evalSA + ef.value(sa.getOptimal());
		        
		        
		        //Genetic Algorithm
		        start = System.nanoTime();
		        StandardGeneticAlgorithm ga = new StandardGeneticAlgorithm(200, 40, 80, gap);
		        fit = new FixedIterationTrainer(ga, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10,9);
		        runTimeForGA = runTimeForGA +  time;
		        evalGA = evalGA + ef.value(ga.getOptimal());
		        
		        
		        //Mimic		        
		        MIMIC mimic = new MIMIC(200, 100, pop);
		        fit = new FixedIterationTrainer(mimic, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10, 9);
		        runTimeForMIMIC = runTimeForMIMIC +  time;
		        evalMIMIC = evalMIMIC + ef.value(mimic.getOptimal());
	        
    		}
    		
    		//Averaging the total optimal fitness values and total time for test runs
    		evalRhc /= runForAvg;
    		evalSA /= runForAvg;
    		evalGA /= runForAvg;
    		evalMIMIC /= runForAvg;
    		
    		runTimeForrhc /= runForAvg;
    		runTimeForsa /= runForAvg;
    		runTimeForGA /= runForAvg;
    		runTimeForMIMIC /= runForAvg;
    		
    		//Updating the values to array
    		rhcRunTime[i] = df.format(runTimeForrhc);
    		listSaTime[i] = df.format(runTimeForsa);
    		listGATime[i] = df.format(runTimeForGA);
    		listMimicTime[i] = df.format(runTimeForMIMIC);
    		
    		listRHCeval[i] = df.format(evalRhc);
    		listSAeval[i] = df.format(evalSA);
    		listGAeval[i] = df.format(evalGA);
    		listMIMICeval[i] = df.format(evalMIMIC);

    }
    	
        
        String fileName = "output/fourpeak_result_accuracy.csv";
        writeToFile(fileName, numOfRestarts, listRHCeval, listSAeval,listGAeval,listMIMICeval);

        
        
        
        String fileName1 = "output/fourpeak_result_runtime.csv";
        writeToFile(fileName1, numOfRestarts, rhcRunTime, listSaTime,listGATime,listMimicTime);

               
        

    }
    
    
    
    
    
    public static void mimic_four_peak() throws IOException{

    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};
    	
    	int samples[] = {10,100, 200, 300, 400, 500, 600};


        	

        	
        	String listRHCeval[][] = new String[samples.length][numOfRestarts.length];
        	for(int z = 0;z<samples.length;z++) {
        		int keep = samples[z];
        	for(int i=0;i < numOfRestarts.length;i++)	
        	{
        		int iterCount = numOfRestarts[i];
        		System.out.println("Start iteration" + iterCount);
        		double evalMIMIC = 0;

        		
        		double runTimeForMIMIC = 0;
      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
        	    	double start, end;
        			System.out.println("Current run" + currun);
    		    	int[] ranges = new int[NN];
    		        Arrays.fill(ranges, 2);
    		        
    		        EvaluationFunction ef = new FourPeaksEvaluationFunction(TT);
    		        Distribution odd = new DiscreteUniformDistribution(ranges);
    		        NeighborFunction nf = new DiscreteChangeOneNeighbor(ranges);
    		        MutationFunction mf = new DiscreteChangeOneMutation(ranges);
    		        CrossoverFunction cf = new SingleCrossOver();
    		        Distribution df = new DiscreteDependencyTree(.1, ranges); 
    		        HillClimbingProblem hcp = new GenericHillClimbingProblem(ef, odd, nf);
    		        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(ef, odd, mf, cf);
    		        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(ef, odd, df);
    		        
    		        //Mimic
    		        start = System.nanoTime();
    		        MIMIC mimic = new MIMIC(200, 20, pop);
    		        FixedIterationTrainer fit = new FixedIterationTrainer(mimic, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10, 9);
    		        runTimeForMIMIC = runTimeForMIMIC +  time;
    		        evalMIMIC = evalMIMIC + ef.value(mimic.getOptimal());
    		        		        

    	        
        		}
        		//Averaging the total optimal fitness values and total time for test runs
        		evalMIMIC /= runForAvg;
        		listRHCeval[z][i] = df.format(evalMIMIC);
        		runTimeForMIMIC /= runForAvg;
        	}
       
        }
        	 
            String fileName = "output/four_peak_result_mimic_keep_count.csv";
            writeToFile(fileName, numOfRestarts, listRHCeval);
    }
      
    
    
    
    
    public static void simulated_annealing_fourpeak() throws IOException{

    	double cooling_exps[] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};

    	double start, end;

        	
        	String listRHCeval[][] = new String[cooling_exps.length][numOfRestarts.length];

        	for(int i=0;i < cooling_exps.length;i++)
        		
        	{
        		double cooling_exp = cooling_exps[i];
        		System.out.println("Cooling Exponent " + cooling_exp);
        		
           		for(int h=0;h < numOfRestarts.length;h++)	
            	{
        			int iterCount = numOfRestarts[h];
        		
        		double evalSA = 0;

        		
        		double runTimeForsa = 0;
      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
        			System.out.println("Current run" + currun);
    		    	int[] ranges = new int[NN];
    		        Arrays.fill(ranges, 2);
    		        
    		        EvaluationFunction ef = new FourPeaksEvaluationFunction(TT);
    		        Distribution odd = new DiscreteUniformDistribution(ranges);
    		        NeighborFunction nf = new DiscreteChangeOneNeighbor(ranges);
    		        MutationFunction mf = new DiscreteChangeOneMutation(ranges);
    		        CrossoverFunction cf = new SingleCrossOver();
    		        Distribution df = new DiscreteDependencyTree(.1, ranges); 
    		        HillClimbingProblem hcp = new GenericHillClimbingProblem(ef, odd, nf);
    		        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(ef, odd, mf, cf);
    		        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(ef, odd, df);
    		        
    		        
    		        //Simulated Annealing
    		        start = System.nanoTime();
    		        SimulatedAnnealing sa = new SimulatedAnnealing(1E11, cooling_exp, hcp);
    		        FixedIterationTrainer fit = new FixedIterationTrainer(sa, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10,9);
    		        runTimeForsa = runTimeForsa +  time;
    		        evalSA = evalSA + ef.value(sa.getOptimal());
    		        		        

    	        
        		}
        		
        		//Averaging the total optimal fitness values and total time for test runs
        		runTimeForsa /= runForAvg;

        		
        		evalSA /= runForAvg;

        		listRHCeval[i][h] = df.format(evalSA);

        }
        	}
            
            String fileName = "output/four_peaks_result_sa_coolin_fitness.csv";
            writeToFile(fileName, numOfRestarts, listRHCeval);

    }    
    
    
    
    public static void rhc_four_peaks() throws IOException{
    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};

        	String listRHCeval[] = new String[numOfRestarts.length];
        	double start, end;
        
        	for(int i=0;i < numOfRestarts.length;i++)
        		
        	{
        		int iterCount = numOfRestarts[i];
        		System.out.println("Start iteration" + iterCount);
        		double evalRhc = 0;

        		
        		double runTimeForrhc = 0;
      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
        			System.out.println("Current run" + currun);
    		    	int[] ranges = new int[NN];
    		        Arrays.fill(ranges, 2);
    		        
    		        EvaluationFunction ef = new FourPeaksEvaluationFunction(TT);
    		        Distribution odd = new DiscreteUniformDistribution(ranges);
    		        NeighborFunction nf = new DiscreteChangeOneNeighbor(ranges);
    		        MutationFunction mf = new DiscreteChangeOneMutation(ranges);
    		        CrossoverFunction cf = new SingleCrossOver();
    		        Distribution df = new DiscreteDependencyTree(.1, ranges); 
    		        HillClimbingProblem hcp = new GenericHillClimbingProblem(ef, odd, nf);
    		        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(ef, odd, mf, cf);
    		        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(ef, odd, df);
    		        
    		        
    		        //Randomized hill climbing
    		        start = System.nanoTime();
    		        RandomizedHillClimbing rhc = new RandomizedHillClimbing(hcp);      
    		        FixedIterationTrainer fit = new FixedIterationTrainer(rhc, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10,9);
    		        runTimeForrhc = runTimeForrhc +  time;
    		        evalRhc = evalRhc + ef.value(rhc.getOptimal());
    		        		  
        		}
        		
        		//Averaging the total optimal fitness values and total time for test runs
        		evalRhc /= runForAvg;

        		
        		runTimeForrhc /= runForAvg;

        		listRHCeval[i] = df.format(evalRhc);


        }
            
            String fileName = "output/four_peaks_result_rhc_fitness.csv";
            writeToFile(fileName, numOfRestarts, listRHCeval);
            

                   
        }
        
        
    
    
 
    public static void ga_fourpeak() throws IOException {


    	int cross_overs_n_mutate[] = {20,40,60,80,100,120,140,160,180, 200};

        	
        	String listRHCeval[][] = new String[cross_overs_n_mutate.length][cross_overs_n_mutate.length];

        	int iterCount = 2000;
        	for(int i=0;i < cross_overs_n_mutate.length;i++)	
        	{
        		
        		double evalGA = 0;
        		int cross_over = cross_overs_n_mutate[i];
        		
        		double runTimeForGA = 0;
        		System.out.println("Cross Over " + cross_over);
        		for(int k=0;k< cross_overs_n_mutate.length;k++) {
        		

        		int mutate = cross_overs_n_mutate[k];
        		
        		System.out.println("Mutate " + mutate);

      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
        	    	double start, end;
        			System.out.println("Current run" + currun);
    		    	int[] ranges = new int[NN];
    		        Arrays.fill(ranges, 2);
    		        
    		        EvaluationFunction ef = new FourPeaksEvaluationFunction(TT);
    		        Distribution odd = new DiscreteUniformDistribution(ranges);
    		        NeighborFunction nf = new DiscreteChangeOneNeighbor(ranges);
    		        MutationFunction mf = new DiscreteChangeOneMutation(ranges);
    		        CrossoverFunction cf = new SingleCrossOver();
    		        Distribution df = new DiscreteDependencyTree(.1, ranges); 
    		        HillClimbingProblem hcp = new GenericHillClimbingProblem(ef, odd, nf);
    		        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(ef, odd, mf, cf);
    		        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(ef, odd, df);
    		        
    		        
    		        //Genetic Algorithm
    		        
    		        //Genetic Algorithm
    		        start = System.nanoTime();
    		        StandardGeneticAlgorithm ga = new StandardGeneticAlgorithm(200, cross_over, mutate, gap);
    		        FixedIterationTrainer fit = new FixedIterationTrainer(ga, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10,9);
    		        runTimeForGA = runTimeForGA +  time;
    		        evalGA = evalGA + ef.value(ga.getOptimal());
    		        		        

    	        
        		}
        		
        		//Averaging the total optimal fitness values and total time for test runs
        		runTimeForGA /= runForAvg;

        		
        		evalGA /= runForAvg;
        		listRHCeval[i][k] = df.format(evalGA);
        		}
        		

        }
        	
            
            String fileName = "output/four_peaks_result_ga_cross_over_mutate.csv";
            writeToFile(fileName, cross_overs_n_mutate, listRHCeval);	
    
    }  
    
    
public static void fourPeak() throws IOException {
        
    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};


    	
    	String rhcRunTime[] = new String[numOfRestarts.length];
    	String listSaTime[] = new String[numOfRestarts.length];
    	String listGATime[] = new String[numOfRestarts.length];
    	String listMimicTime[] = new String[numOfRestarts.length];
    	
    	String listRHCeval[] = new String[numOfRestarts.length];
    	String listSAeval[] = new String[numOfRestarts.length];
    	String listGAeval[] = new String[numOfRestarts.length];
    	String listMIMICeval[] = new String[numOfRestarts.length];
    	

    	for(int i=0;i < numOfRestarts.length;i++)
    		
    	{
    		int iterCount = numOfRestarts[i];
    		System.out.println("Start iteration" + iterCount);
    		double evalRhc = 0;
    		double evalSA = 0;
    		double evalGA = 0;
    		double evalMIMIC = 0;
    		
    		double runTimeForrhc = 0;
    		double runTimeForsa = 0;
    		double runTimeForGA = 0;
    		double runTimeForMIMIC = 0;
    		double time = 0;
        	int runForAvg = 5;
    		for (int currun= 0; currun < runForAvg; currun++)
    		{
    	    	double start, end;
    			System.out.println("Start run" + currun);
    			
    			Random random = new Random();
    	        // create the random points
    	        double[][] points = new double[N][2];
    	        for (int j = 0; j < points.length; j++) {
    	            points[j][0] = random.nextDouble();
    	            points[j][1] = random.nextDouble();   
    	        }
		        
    	        TravelingSalesmanEvaluationFunction tspef = new TravelingSalesmanRouteEvaluationFunction(points);
    	        Distribution odd = new DiscretePermutationDistribution(N);
    	        NeighborFunction nf = new SwapNeighbor();
    	        MutationFunction mf = new SwapMutation();
    	        CrossoverFunction cf = new TravelingSalesmanCrossOver(tspef);
    	        HillClimbingProblem hcp = new GenericHillClimbingProblem(tspef, odd, nf);
    	        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(tspef, odd, mf, cf);
		        
		        
		        //Randomized hill climbing
		        start = System.nanoTime();
		        RandomizedHillClimbing rhc = new RandomizedHillClimbing(hcp);      
		        FixedIterationTrainer fit = new FixedIterationTrainer(rhc, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10,9);
		        runTimeForrhc = runTimeForrhc +  time;
		        evalRhc = evalRhc + tspef.value(rhc.getOptimal());
		        		        
		        
		        //Simulated Annealing
		        start = System.nanoTime();
		        SimulatedAnnealing sa = new SimulatedAnnealing(1E11, 0.7, hcp);
		        fit = new FixedIterationTrainer(sa, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10,9);
		        runTimeForsa = runTimeForsa +  time;
		        evalSA = evalSA + tspef.value(sa.getOptimal());
		        
		        
		        //Genetic Algorithm
		        start = System.nanoTime();
		        StandardGeneticAlgorithm ga = new StandardGeneticAlgorithm(200, 40, 80, gap);
		        fit = new FixedIterationTrainer(ga, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10,9);
		        runTimeForGA = runTimeForGA +  time;
		        evalGA = evalGA + tspef.value(ga.getOptimal());
		        
		        
		        //Mimic
		        start = System.nanoTime();
		        tspef = new TravelingSalesmanSortEvaluationFunction(points);
		        int[] ranges = new int[N];
		        Arrays.fill(ranges, N);
		        odd = new  DiscreteUniformDistribution(ranges);
		        Distribution df = new DiscreteDependencyTree(.1, ranges); 
		        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(tspef, odd, df);
		        
		        MIMIC mimic = new MIMIC(200, 100, pop);
		        fit = new FixedIterationTrainer(mimic, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10, 9);
		        runTimeForMIMIC = runTimeForMIMIC +  time;
		        evalMIMIC = evalMIMIC + tspef.value(mimic.getOptimal());
	        
    		}
    		
    		//Averaging the total optimal fitness values and total time for test runs
    		evalRhc /= runForAvg;
    		evalSA /= runForAvg;
    		evalGA /= runForAvg;
    		evalMIMIC /= runForAvg;
    		
    		runTimeForrhc /= runForAvg;
    		runTimeForsa /= runForAvg;
    		runTimeForGA /= runForAvg;
    		runTimeForMIMIC /= runForAvg;
    		
    		//Updating the values to array
    		rhcRunTime[i] = df.format(runTimeForrhc);
    		listSaTime[i] = df.format(runTimeForsa);
    		listGATime[i] = df.format(runTimeForGA);
    		listMimicTime[i] = df.format(runTimeForMIMIC);
    		
    		listRHCeval[i] = df.format(evalRhc);
    		listSAeval[i] = df.format(evalSA);
    		listGAeval[i] = df.format(evalGA);
    		listMIMICeval[i] = df.format(evalMIMIC);

    }
    	
        
        String fileName = "output/travelling_salesman_result_accuracy.csv";
        writeToFile(fileName, numOfRestarts, listRHCeval, listSAeval,listGAeval,listMIMICeval);

        
        
        
        String fileName1 = "output/travelling_salesman_result_runtime.csv";
        writeToFile(fileName1, numOfRestarts, rhcRunTime, listSaTime,listGATime,listMimicTime);

               
        

    }
    
    
    
    //Travelling Salesman
    public static void tsp() throws IOException {
        
    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};


    	
    	String rhcRunTime[] = new String[numOfRestarts.length];
    	String listSaTime[] = new String[numOfRestarts.length];
    	String listGATime[] = new String[numOfRestarts.length];
    	String listMimicTime[] = new String[numOfRestarts.length];
    	
    	String listRHCeval[] = new String[numOfRestarts.length];
    	String listSAeval[] = new String[numOfRestarts.length];
    	String listGAeval[] = new String[numOfRestarts.length];
    	String listMIMICeval[] = new String[numOfRestarts.length];
    	

    	for(int i=0;i < numOfRestarts.length;i++)
    		
    	{
    		int iterCount = numOfRestarts[i];
    		System.out.println("Start iteration" + iterCount);
    		double evalRhc = 0;
    		double evalSA = 0;
    		double evalGA = 0;
    		double evalMIMIC = 0;
    		
    		double runTimeForrhc = 0;
    		double runTimeForsa = 0;
    		double runTimeForGA = 0;
    		double runTimeForMIMIC = 0;
    		double time = 0;
        	int runForAvg = 5;
    		for (int currun= 0; currun < runForAvg; currun++)
    		{
    	    	double start, end;
    			System.out.println("Start run" + currun);
    			
    			Random random = new Random();
    	        // create the random points
    	        double[][] points = new double[N][2];
    	        for (int j = 0; j < points.length; j++) {
    	            points[j][0] = random.nextDouble();
    	            points[j][1] = random.nextDouble();   
    	        }
		        
    	        TravelingSalesmanEvaluationFunction tspef = new TravelingSalesmanRouteEvaluationFunction(points);
    	        Distribution odd = new DiscretePermutationDistribution(N);
    	        NeighborFunction nf = new SwapNeighbor();
    	        MutationFunction mf = new SwapMutation();
    	        CrossoverFunction cf = new TravelingSalesmanCrossOver(tspef);
    	        HillClimbingProblem hcp = new GenericHillClimbingProblem(tspef, odd, nf);
    	        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(tspef, odd, mf, cf);
		        
		        
		        //Randomized hill climbing
		        start = System.nanoTime();
		        RandomizedHillClimbing rhc = new RandomizedHillClimbing(hcp);      
		        FixedIterationTrainer fit = new FixedIterationTrainer(rhc, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10,9);
		        runTimeForrhc = runTimeForrhc +  time;
		        evalRhc = evalRhc + tspef.value(rhc.getOptimal());
		        		        
		        
		        //Simulated Annealing
		        start = System.nanoTime();
		        SimulatedAnnealing sa = new SimulatedAnnealing(1E11, 0.7, hcp);
		        fit = new FixedIterationTrainer(sa, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10,9);
		        runTimeForsa = runTimeForsa +  time;
		        evalSA = evalSA + tspef.value(sa.getOptimal());
		        
		        
		        //Genetic Algorithm
		        start = System.nanoTime();
		        StandardGeneticAlgorithm ga = new StandardGeneticAlgorithm(200, 40, 80, gap);
		        fit = new FixedIterationTrainer(ga, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10,9);
		        runTimeForGA = runTimeForGA +  time;
		        evalGA = evalGA + tspef.value(ga.getOptimal());
		        
		        
		        //Mimic
		        start = System.nanoTime();
		        tspef = new TravelingSalesmanSortEvaluationFunction(points);
		        int[] ranges = new int[N];
		        Arrays.fill(ranges, N);
		        odd = new  DiscreteUniformDistribution(ranges);
		        Distribution df = new DiscreteDependencyTree(.1, ranges); 
		        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(tspef, odd, df);
		        
		        MIMIC mimic = new MIMIC(200, 100, pop);
		        fit = new FixedIterationTrainer(mimic, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10, 9);
		        runTimeForMIMIC = runTimeForMIMIC +  time;
		        evalMIMIC = evalMIMIC + tspef.value(mimic.getOptimal());
	        
    		}
    		
    		//Averaging the total optimal fitness values and total time for test runs
    		evalRhc /= runForAvg;
    		evalSA /= runForAvg;
    		evalGA /= runForAvg;
    		evalMIMIC /= runForAvg;
    		
    		runTimeForrhc /= runForAvg;
    		runTimeForsa /= runForAvg;
    		runTimeForGA /= runForAvg;
    		runTimeForMIMIC /= runForAvg;
    		
    		//Updating the values to array
    		rhcRunTime[i] = df.format(runTimeForrhc);
    		listSaTime[i] = df.format(runTimeForsa);
    		listGATime[i] = df.format(runTimeForGA);
    		listMimicTime[i] = df.format(runTimeForMIMIC);
    		
    		listRHCeval[i] = df.format(evalRhc);
    		listSAeval[i] = df.format(evalSA);
    		listGAeval[i] = df.format(evalGA);
    		listMIMICeval[i] = df.format(evalMIMIC);

    }
    	
        
        String fileName = "output/travelling_salesman_result_accuracy.csv";
        writeToFile(fileName, numOfRestarts, listRHCeval, listSAeval,listGAeval,listMIMICeval);

        
        
        
        String fileName1 = "output/travelling_salesman_result_runtime.csv";
        writeToFile(fileName1, numOfRestarts, rhcRunTime, listSaTime,listGATime,listMimicTime);

               
        

    }
    
    
    public static void ga_tsp() throws IOException {


    	int cross_overs_n_mutate[] = {20,40,60,80,100,120,140,160,180, 200};

        	
        	String listRHCeval[][] = new String[cross_overs_n_mutate.length][cross_overs_n_mutate.length];

        	int iterCount = 2000;
        	for(int i=0;i < cross_overs_n_mutate.length;i++)	
        	{
        		
        		double evalGA = 0;
        		int cross_over = cross_overs_n_mutate[i];
        		
        		double runTimeForGA = 0;
        		System.out.println("Cross Over " + cross_over);
        		for(int k=0;k< cross_overs_n_mutate.length;k++) {
        		

        		int mutate = cross_overs_n_mutate[k];
        		
        		System.out.println("Mutate " + mutate);

      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
        	    	double start, end;
        			System.out.println("Start run" + currun);
        			
        			Random random = new Random();
        	        // create the random points
        	        double[][] points = new double[N][2];
        	        for (int j = 0; j < points.length; j++) {
        	            points[j][0] = random.nextDouble();
        	            points[j][1] = random.nextDouble();   
        	        }
    		        
        	        TravelingSalesmanEvaluationFunction tspef = new TravelingSalesmanRouteEvaluationFunction(points);
        	        Distribution odd = new DiscretePermutationDistribution(N);
        	        NeighborFunction nf = new SwapNeighbor();
        	        MutationFunction mf = new SwapMutation();
        	        CrossoverFunction cf = new TravelingSalesmanCrossOver(tspef);
        	        HillClimbingProblem hcp = new GenericHillClimbingProblem(tspef, odd, nf);
        	        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(tspef, odd, mf, cf);
    		        
    		        
    		        //Genetic Algorithm
    		        start = System.nanoTime();
    		        StandardGeneticAlgorithm ga = new StandardGeneticAlgorithm(200, cross_over, mutate, gap);
    		        FixedIterationTrainer fit = new FixedIterationTrainer(ga, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10,9);
    		        runTimeForGA = runTimeForGA +  time;
    		        evalGA = evalGA + tspef.value(ga.getOptimal());
    		        		        

    	        
        		}
        		
        		//Averaging the total optimal fitness values and total time for test runs
        		runTimeForGA /= runForAvg;

        		
        		evalGA /= runForAvg;
        		listRHCeval[i][k] = df.format(evalGA);
        		}
        		

        }
        	
            
            String fileName = "output/travelling_salesman_result_ga_cross_over_mutate.csv";
            writeToFile(fileName, cross_overs_n_mutate, listRHCeval);	
    
    }
    
    
        public static void simulated_annealing_tsp() throws IOException{

    	double cooling_exps[] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.0};

    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};

        	

        	
        	String listRHCeval[][] = new String[cooling_exps.length][numOfRestarts.length];

        
        	for(int i=0;i < cooling_exps.length;i++)
        		
        	{
        		double cooling_exp = cooling_exps[i];
        		System.out.println("Cooling Exponent " + cooling_exp);
        		
        		for(int h=0;h < numOfRestarts.length;h++)	
            	{
        			int iterCount = numOfRestarts[h];
        		
        		
        		double evalSA = 0;

        		
        		double runTimeForsa = 0;
      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
        	    	double start, end;
        			System.out.println("Start run" + currun);
        			
        			Random random = new Random();
        	        // create the random points
        	        double[][] points = new double[N][2];
        	        for (int j = 0; j < points.length; j++) {
        	            points[j][0] = random.nextDouble();
        	            points[j][1] = random.nextDouble();   
        	        }
    		        
        	        TravelingSalesmanEvaluationFunction tspef = new TravelingSalesmanRouteEvaluationFunction(points);
        	        Distribution odd = new DiscretePermutationDistribution(N);
        	        NeighborFunction nf = new SwapNeighbor();
        	        MutationFunction mf = new SwapMutation();
        	        CrossoverFunction cf = new TravelingSalesmanCrossOver(tspef);
        	        HillClimbingProblem hcp = new GenericHillClimbingProblem(tspef, odd, nf);
        	        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(tspef, odd, mf, cf);
    		        
    		        
        	        //Simulated Annealing
    		        start = System.nanoTime();
    		        SimulatedAnnealing sa = new SimulatedAnnealing(1E11, cooling_exp, hcp);
    		        FixedIterationTrainer fit = new FixedIterationTrainer(sa, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10,9);
    		        runTimeForsa = runTimeForsa +  time;
    		        evalSA = evalSA + tspef.value(sa.getOptimal());
    		        		        

    	        
        		}
        		
        		//Averaging the total optimal fitness values and total time for test runs
        		runTimeForsa /= runForAvg;

        		
        		evalSA /= runForAvg;

        		listRHCeval[i][h] = df.format(evalSA);
            	}

        }
        	
            
            String fileName = "output/travelling_salesman_result_sa_coolin_fitness.csv";
            writeToFile(fileName, numOfRestarts, listRHCeval);

    }
    
        
        public static void writeToFile(String fileName, int[] header, String[][] rhc) throws IOException {
        	
            File file1 = new File(fileName);
            
            // creates the file
            file1.delete();
//            FileWriter fw = new FileWriter(file);
            BufferedWriter writer1 = Files.newBufferedWriter(
                    Paths.get(fileName), 
                    StandardOpenOption.APPEND, 
                    StandardOpenOption.CREATE);
            CSVPrinter printer1 = new CSVPrinter(writer1, CSVFormat.DEFAULT
            		  .withHeader(Arrays.stream(header)
                              .mapToObj(String::valueOf)
                              .toArray(String[]::new)));

            for(String[] rh:rhc) {
            printer1.printRecord(rh);
            }

            
            
            printer1.flush();
            printer1.close();
            System.out.println("Writing to the file finished "+fileName);
        	
        }    
        
    
    public static void rhc_tsp() throws IOException{
	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};


    	

    	
    	String listRHCeval[] = new String[numOfRestarts.length];

    
    	for(int i=0;i < numOfRestarts.length;i++)
    		
    	{
    		int iterCount = numOfRestarts[i];
    		System.out.println("Start iteration" + iterCount);
    		double evalRhc = 0;

    		
    		double runTimeForrhc = 0;
  
    		double time = 0;
        	int runForAvg = 5;
    		for (int currun= 0; currun < runForAvg; currun++)
    		{
    	    	double start, end;
    			System.out.println("Start run" + currun);
    			
    			Random random = new Random();
    	        // create the random points
    	        double[][] points = new double[N][2];
    	        for (int j = 0; j < points.length; j++) {
    	            points[j][0] = random.nextDouble();
    	            points[j][1] = random.nextDouble();   
    	        }
		        
    	        TravelingSalesmanEvaluationFunction tspef = new TravelingSalesmanRouteEvaluationFunction(points);
    	        Distribution odd = new DiscretePermutationDistribution(N);
    	        NeighborFunction nf = new SwapNeighbor();
    	        MutationFunction mf = new SwapMutation();
    	        CrossoverFunction cf = new TravelingSalesmanCrossOver(tspef);
    	        HillClimbingProblem hcp = new GenericHillClimbingProblem(tspef, odd, nf);
    	        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(tspef, odd, mf, cf);
		        
		        
		        //Randomized hill climbing
		        start = System.nanoTime();
		        RandomizedHillClimbing rhc = new RandomizedHillClimbing(hcp);      
		        FixedIterationTrainer fit = new FixedIterationTrainer(rhc, iterCount);
		        fit.train();
		        end = System.nanoTime();
		        time = (end - start)/Math.pow(10,9);
		        runTimeForrhc = runTimeForrhc +  time;
		        evalRhc = evalRhc + tspef.value(rhc.getOptimal());
		        		        

	        
    		}
    		
    		//Averaging the total optimal fitness values and total time for test runs
    		evalRhc /= runForAvg;

    		
    		runTimeForrhc /= runForAvg;

    		listRHCeval[i] = df.format(evalRhc);


    }
    	
        
        String fileName = "output/travelling_salesman_result_rhc_fitness.csv";
        writeToFile(fileName, numOfRestarts, listRHCeval);

        
        


               
    }
    
    
    
    public static void mimic_tsp() throws IOException{

    	int numOfRestarts[] = {10, 100, 250, 500, 1000, 2000};
    	
    	int samples[] = {10,100, 200, 300, 400, 500, 600};


        	

        	
        	String listRHCeval[][] = new String[samples.length][numOfRestarts.length];
        	for(int z = 0;z<samples.length;z++) {
        		int keep = samples[z];
        	for(int i=0;i < numOfRestarts.length;i++)	
        	{
        		int iterCount = numOfRestarts[i];
        		System.out.println("Start iteration" + iterCount);
        		double evalMIMIC = 0;

        		
        		double runTimeForMIMIC = 0;
      
        		double time = 0;
            	int runForAvg = 5;
        		for (int currun= 0; currun < runForAvg; currun++)
        		{
        	    	double start, end;
        			System.out.println("Start run" + currun);
        			
        			Random random = new Random();
        	        // create the random points
        	        double[][] points = new double[N][2];
        	        for (int j = 0; j < points.length; j++) {
        	            points[j][0] = random.nextDouble();
        	            points[j][1] = random.nextDouble();   
        	        }
    		        
        	        TravelingSalesmanEvaluationFunction tspef = new TravelingSalesmanRouteEvaluationFunction(points);
        	        Distribution odd = new DiscretePermutationDistribution(N);
        	        NeighborFunction nf = new SwapNeighbor();
        	        MutationFunction mf = new SwapMutation();
        	        CrossoverFunction cf = new TravelingSalesmanCrossOver(tspef);
        	        HillClimbingProblem hcp = new GenericHillClimbingProblem(tspef, odd, nf);
        	        GeneticAlgorithmProblem gap = new GenericGeneticAlgorithmProblem(tspef, odd, mf, cf);
    		        
    		        
    		        //Mimic
    		        start = System.nanoTime();
    		        tspef = new TravelingSalesmanSortEvaluationFunction(points);
    		        int[] ranges = new int[N];
    		        Arrays.fill(ranges, N);
    		        odd = new  DiscreteUniformDistribution(ranges);
    		        Distribution df = new DiscreteDependencyTree(.1, ranges); 
    		        ProbabilisticOptimizationProblem pop = new GenericProbabilisticOptimizationProblem(tspef, odd, df);
    		        
    		        MIMIC mimic = new MIMIC(1500, keep, pop);
    		        FixedIterationTrainer fit = new FixedIterationTrainer(mimic, iterCount);
    		        fit.train();
    		        end = System.nanoTime();
    		        time = (end - start)/Math.pow(10, 9);
    		        runTimeForMIMIC = runTimeForMIMIC +  time;
    		        evalMIMIC = evalMIMIC + tspef.value(mimic.getOptimal());
    		        		        

    	        
        		}
        		//Averaging the total optimal fitness values and total time for test runs
        		evalMIMIC /= runForAvg;
        		listRHCeval[z][i] = df.format(evalMIMIC);
        		runTimeForMIMIC /= runForAvg;
        	}
       
        }
        	
            
            String fileName = "output/travelling_salesman_result_mimic_keep_count.csv";
            writeToFile(fileName, numOfRestarts, listRHCeval);
    }
    
    
 public static void writeToFile(String fileName, int[] header, String[] rhc) throws IOException {
    	
        File file1 = new File(fileName);
        
        // creates the file
        file1.delete();
//        FileWriter fw = new FileWriter(file);
        BufferedWriter writer1 = Files.newBufferedWriter(
                Paths.get(fileName), 
                StandardOpenOption.APPEND, 
                StandardOpenOption.CREATE);
        CSVPrinter printer1 = new CSVPrinter(writer1, CSVFormat.DEFAULT
        		  .withHeader(Arrays.stream(header)
                          .mapToObj(String::valueOf)
                          .toArray(String[]::new)));


        printer1.printRecord(rhc);

        
        
        printer1.flush();
        printer1.close();
    	
    }
 
 
 
 public static void writeToFile(String fileName, double[] header, String[] rhc) throws IOException {
 	
     File file1 = new File(fileName);
     
     // creates the file
     file1.delete();
//     FileWriter fw = new FileWriter(file);
     BufferedWriter writer1 = Files.newBufferedWriter(
             Paths.get(fileName), 
             StandardOpenOption.APPEND, 
             StandardOpenOption.CREATE);
     CSVPrinter printer1 = new CSVPrinter(writer1, CSVFormat.DEFAULT
     		  .withHeader(Arrays.stream(header)
                       .mapToObj(String::valueOf)
                       .toArray(String[]::new)));


     printer1.printRecord(rhc);

     
     
     printer1.flush();
     printer1.close();
 	
 }

    
    
    public static void writeToFile(String fileName, int[] header, String[] rhc, String[] sa, String[] ga, String[] mimic) throws IOException {
    	
        File file1 = new File(fileName);
        
        // creates the file
        file1.delete();
//        FileWriter fw = new FileWriter(file);
        BufferedWriter writer1 = Files.newBufferedWriter(
                Paths.get(fileName), 
                StandardOpenOption.APPEND, 
                StandardOpenOption.CREATE);
        CSVPrinter printer1 = new CSVPrinter(writer1, CSVFormat.DEFAULT
        		  .withHeader(Arrays.stream(header)
                          .mapToObj(String::valueOf)
                          .toArray(String[]::new)));


        printer1.printRecord(rhc);
        printer1.printRecord(sa);
        printer1.printRecord(ga);
        printer1.printRecord(mimic);
        
        
        printer1.flush();
        printer1.close();
    	
    }
    
    
    
    
}
