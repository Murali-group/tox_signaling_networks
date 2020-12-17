import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import org.apache.commons.cli.*;

public class EdgeLinker {
	public static void main(String[] args)  throws IOException {
		
		Options options = new Options();
		
		Option network = new Option("n", "network", true, "Input interactome");
		network.setRequired(true);
		options.addOption(network);
		
		Option startends = new Option("s", "start-ends", true, "File specifying the receptors and tfs or nodes to use as sources and targets");
		startends.setRequired(true);
		options.addOption(startends);

		Option outprefix = new Option("o", "out-prefix", true, "path/to/prefix for the output files");
		outprefix.setRequired(true);
		options.addOption(outprefix);

		Option klimit = new Option("k", "max-k", true, "Maximum number of paths to output. Default is to write all paths and edges");
		options.addOption(klimit);

		Option multirun = new Option("m", "multi-run", false, "Option to read multiple start-ends file from the --start-ends option");
		options.addOption(multirun);

		Option edgepenalty = new Option("e", "edge-penalty", true, "Add the natural log of the specified penalty to the cost of each edge. This will effectively increase the cost of each path by the length * edge penalty");
		options.addOption(edgepenalty);
		
		Option startEndsPenaltyOpt = new Option("p", "start-ends-penalty", false, "Set the cost specified in the third column of the start-ends file to be the cost of the super-source->start and end->super-target edges. Otherwise will be 0");
		options.addOption(startEndsPenaltyOpt);
		
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("EdgeLinker", options);

            System.exit(1);
            return;
        }

		//Read in arguments
		String graphFileName = cmd.getOptionValue("network");
		String startEndsArg = cmd.getOptionValue("start-ends");
		String outputPrefixesArg = cmd.getOptionValue("out-prefix");
		ArrayList<String> startEnds = new ArrayList<String>();
		ArrayList<String> outputPrefixes = new ArrayList<String>();
        // default is to include a path for every edge in the network. (100 million should be big enough)
        int maxk = 100000000;
        if (cmd.hasOption("max-k")){
            maxk = Integer.parseInt(cmd.getOptionValue("max-k"));
        }
        long startTime;
        double inputTime; 
        double algorithmTime = 0;
        double fileWritingTime = 0;
        boolean verbose = true;
        double edgePenalty;
        if (cmd.hasOption("edge-penalty")){
			edgePenalty = Double.parseDouble(cmd.getOptionValue("edge-penalty"));
        }
        else {
        	// be default, an edge penalty of 1 will not add anything to the cost of an edge (log(1) is 0)
        	edgePenalty = 1;
        }
        boolean startEndsPenalty = cmd.hasOption("start-ends-penalty");
        boolean multiRun = cmd.hasOption("multi-run");
		
		//Read input; build graph.
        startTime = System.nanoTime();
		InputReader input = new InputReader(graphFileName, edgePenalty);
        inputTime = (System.nanoTime() - startTime) / 1000000;
        System.out.println("Time to parse input graph: " + inputTime / 1000 + " sec");

		// get the startEnds and outputPrefixes from the arguments
        if (multiRun){
        	System.out.println("Reading the startEnds and ouputPrefixes from the specified files");
        	// parse the startEnds list and outputPrefixes list from the specified input files
        	startEnds = input.ParseFileList(startEndsArg);
        	outputPrefixes = input.ParseFileList(outputPrefixesArg);
            if (startEnds.size() != outputPrefixes.size()){
                System.out.println("Error: # of startEnds does not equal # of outputPrefixes: " + startEnds.size() + " " + outputPrefixes.size());
                System.exit(1);
            }

            // don't print everything out to the log file
            verbose = false;
            for (int i = 0; i < startEnds.size(); i++){
                // first make sure all of the startEnds have sources and targets in the network
                // read the start end file, which also adds the edges from super-source to sources and from targets to super-target
                input.AddStartEnd(startEnds.get(i), startEndsPenalty, verbose);
                // remove the sources and targets from the graph for the next run 
                input.RemoveStartEnd();
            }
            // we already printed the stats abou the # of sources and targets in the network
            System.out.println("Running EdgeLinker on " + startEnds.size() + " sets of sources and targets");
        }
        else{
        	System.out.println("Reading the sources and targets from the file: " + startEndsArg);
        	// put the single startEnd and outputPrefixes into a list
        	startEnds.add(startEndsArg);
        	outputPrefixes.add(outputPrefixesArg);
        }

        for (int i = 0; i < startEnds.size(); i++){
            // read the start end file
            input.AddStartEnd(startEnds.get(i), startEndsPenalty, verbose);
		
            //Execute the algorithm.
            startTime = System.nanoTime();
            Algorithm execute = new Algorithm(input, maxk);
            execute.run();
            // keep adding up the time taken to run the algorithm
            algorithmTime += (System.nanoTime() - startTime) / 1000000;
            //System.out.println("Time (ms) taken to run: " + (System.nanoTime() - startTime) / 1000000000);

            //Output results to file.
            startTime = System.nanoTime();
            OutputWriter print = new OutputWriter(execute, outputPrefixes.get(i));
            print.printToFile();		
            fileWritingTime += (System.nanoTime() - startTime) / 1000000;
            //System.out.println("Time (ms) taken to write output: " + (System.nanoTime() - startTime) / 1000000000);

            // remove the sources and targets from the graph for the next run 
            input.RemoveStartEnd();
        }
        // divide by 1000 to get the time in seconds 
        algorithmTime = algorithmTime / 1000;
        fileWritingTime = fileWritingTime / 1000;
        System.out.println("Total time to run algorithm: " + algorithmTime + " sec. Avg per run: " + algorithmTime / startEnds.size() + " sec");
        System.out.println("Total time to write output: " + fileWritingTime + " sec. Avg per run: " + fileWritingTime / startEnds.size() + " sec");
	}
}
