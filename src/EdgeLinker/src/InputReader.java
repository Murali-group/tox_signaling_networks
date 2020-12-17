import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

public class InputReader {
	Scanner graphIn;
	Scanner startEnd;

	// These hashmaps map strings for proteins like "P04355" to integers, and
	// vice versa
	// Overflow warning: These integers should never be multiplied, without
	// converting to longs.
	HashMap<String, Integer> mapToInt;
	HashMap<Integer, String> reverseMap;
	
	// List of the starting (source) nodes and ending (target) nodes (in integer space)
	ArrayList<Integer> starts;
	ArrayList<Integer> ends;

	// Map of edge number to the edge start point and end point, as well as
	// cost. (order of file read)
	ArrayList<Integer> edgeEndSet;
	ArrayList<Integer> edgeStartSet;
	ArrayList<Double> weights;
	ArrayList<Edge>[] edges;
	ArrayList<Edge>[] reverseEdges;

	// Negative log transform map of edge costs.
	// Key for edge is a perfect hash of the start and end point. (use given
	// function. Works for up to 999,999,999 nodes)
	HashMap<Long, Double> edgeCost;
	
	// Add the log of the specified penalty to the cost of each edge. 
	// This will effectively increase the cost of each path by the length * log(edge penalty)
	double edgePenalty;

    // Option to split family nodes (proteins joined by a comma) that contain a source or target 
    // into the individual proteins. 
    // For example: suppose a source node X is present in the family node X,Y,Z.
    // All edges containing X,Y,Z will be split (such as X,Y,Z -> A split to X->A, Y->A and Z->A)
    //boolean splitFamilyNodes;

	public InputReader(String graphFile, double edgePenalty)
			throws IOException {
		graphIn = new Scanner(new File(graphFile));

		// These hashmaps map strings for proteins like "P04355" to integers,
		// and vice versa
		// Overflow warning: These integers should never be multiplied, without
		// converting to longs.
		mapToInt = new HashMap<String, Integer>();
		reverseMap = new HashMap<Integer, String>();
		
		// Map of edge number to the edge start point and end point, as well as
		// weight of the edge.
		edgeEndSet = new ArrayList<Integer>();
		edgeStartSet = new ArrayList<Integer>();
		weights = new ArrayList<Double>();

		// Negative log transform map of edge costs.
		// Key for edge is a perfect hash of the start and end point. (use given
		// function. Works for up to 999,999,999 nodes)
		edgeCost = new HashMap<Long, Double>();
		
		// Add the log of the specified penalty to the cost of each edge. 
		// This will effectively increase the cost of each path by the length * edge penalty
		this.edgePenalty = edgePenalty;

        // Option to split family nodes (proteins joined by a comma) that contain a source or target 
        //this.splitFamilyNodes = splitFamilyNodes;
		read();
	}

	@SuppressWarnings("unchecked")
	private void read() {
		int index = 2;
        // 'receptor' acts as the super source, 'tf' as the super target
		mapToInt.put("receptor", 0);
		reverseMap.put(0, "receptor");
		mapToInt.put("tf", 1);
		reverseMap.put(1, "tf");

        // read in the network
		while (graphIn.hasNext()) {
			String start = graphIn.next();
			// skip lines that are commented out like the header line
			if (start.startsWith("#")){
	            graphIn.nextLine();
				continue;
			}
			String end = graphIn.next();
			double weight = graphIn.nextDouble();
			graphIn.nextLine();
			if (!mapToInt.containsKey(start)) {
				mapToInt.put(start, index);
				reverseMap.put(index, start);
				index++;
			}
			if (!mapToInt.containsKey(end)) {
				mapToInt.put(end, index);
				reverseMap.put(index, end);
				index++;
			}
			edgeStartSet.add(mapToInt.get(start));
			edgeEndSet.add(mapToInt.get(end));
			weights.add(weight);
		}

        // initialize the edges and reverseEdges lists for each node
		int num_nodes = mapToInt.size();
		int num_edges = weights.size();
		edges = new ArrayList[num_nodes];
		reverseEdges = new ArrayList[num_nodes];
		for (int i = 0; i < num_nodes; i++) {
			edges[i] = new ArrayList<Edge>();
			reverseEdges[i] = new ArrayList<Edge>();
		}

		for (int i = 0; i < num_edges; i++) {
			int start = edgeStartSet.get(i);
			int end = edgeEndSet.get(i);
            // set the cost of the edge as the -log of the edge weight
            // log is the natural log by default
			double cost = (Math.log(weights.get(i)) * -1.0);
            if (cost < 0){
                System.out.println("Error: invalid weight for edge " + reverseMap.get(start) + "->" + reverseMap.get(end) + ": " + weights.get(i));
                System.out.println("Must be between 0 and 1. Quitting.");
                System.exit(1);
            }
            // add the edge penalty. By default, it is 1 (which will be 0 after the log)
            cost = cost + Math.log(edgePenalty);
			edges[start].add(new Edge(end, cost));
			edgeCost.put(hash(start, end), cost);
			reverseEdges[end].add(new Edge(start, cost));
		}

	}

    public void AddStartEnd(String startEndFile, boolean startEndsPenalty, boolean verbose)
			throws IOException {
		// Use these scanners to read files
		//System.out.println(startEndFile);
		startEnd = new Scanner(new File(startEndFile));
		starts = new ArrayList<Integer>();
		ends = new ArrayList<Integer>();
        int num_sources = 0;
        int num_targets = 0;
		while (startEnd.hasNext()) {
            // first column (start) contains the node name
            // second column (end) contains either 'receptor' or 'tf
			String nodeStr = startEnd.next();
			// skip lines that are commented out like the header line
			if (nodeStr.startsWith("#")){
	            startEnd.nextLine();
				continue;
			}
			String nodeTypeStr = startEnd.next();
			// by default, the cost of a super-source or super-target edge is 0
			double cost = .00000000000000001;
			if (startEndsPenalty) {
				cost = Double.parseDouble(startEnd.next());
			}
            // end of line. skip to next line
            startEnd.nextLine();
            if (mapToInt.containsKey(nodeTypeStr)) {
                if (mapToInt.get(nodeTypeStr) == 0) {
                    num_sources += 1;
                }
                else {
                    num_targets += 1;
                }
            }
			if (mapToInt.containsKey(nodeStr) && mapToInt.containsKey(nodeTypeStr)) {
                // start is the node name
				int node = mapToInt.get(nodeStr);
                // end is either 'receptor' or 'tf' 
				int rec_or_tf = mapToInt.get(nodeTypeStr);

                // 'receptor' was given the integer value 0, 'tf' was given a 1
                // 'receptor' acts as the super source, 'tf' as the super target
				if (rec_or_tf == 0) {
                    // add an edge from "receptor" to the source (node)
					edges[rec_or_tf].add(new Edge(node, cost));
                    // add a reverse edge from the source to 'receptor'
					reverseEdges[node].add(new Edge(rec_or_tf, cost));

					edgeCost.put(hash(node, rec_or_tf), cost);
                    starts.add(node);
				} 
				// or from the target (node) to "tf"
				else {
					edges[node].add(new Edge(rec_or_tf, cost));
					reverseEdges[rec_or_tf].add(new Edge(node, cost));

					edgeCost.put(hash(node, rec_or_tf), cost);
                    ends.add(node);
				}
			}
		}
        if (verbose){
            System.out.println("Reading sources and targets from file: " + startEndFile);
            System.out.println(starts.size() + "/" + num_sources + " receptors and " + 
                    ends.size() + "/" + num_targets + " tfs were in the network");
        }
        if (starts.size() == 0) {
            System.out.println("Error: No sources were found to connect to the super-source for file: " + startEndFile + ". Quitting");
            System.exit(1);
        }
        if (ends.size() == 0) {
            System.out.println("Error: No targets were found to connect to the super-target for file: " + startEndFile + ". Quitting");
            System.exit(1);
        }
    }
    
    // remove the "receptor" and "tf" super-source and super-target edges from the edges and reverseEdges lists.
    public void RemoveStartEnd(){
    	edges[0] = new ArrayList<Edge>();
    	reverseEdges[1] = new ArrayList<Edge>();
        // also remove the reverse edge to the super source from the source
        for (int i = 0; i < starts.size(); i++){
            int to_super_source = -1;
            for (int e = 0; e < reverseEdges[starts.get(i)].size(); e++){
                if (reverseEdges[starts.get(i)].get(e).end == 0){
                    to_super_source = e;
                }
            }
            reverseEdges[starts.get(i)].remove(to_super_source);
        }
        for (int i = 0; i < ends.size(); i++){
            int to_super_target = -1;
            for (int e = 0; e < edges[ends.get(i)].size(); e++){
                if (edges[ends.get(i)].get(e).end == 1){
                    to_super_target = e;
                }
            }
            edges[ends.get(i)].remove(to_super_target);
        }
    }
    
    // function to parse the input file and return it as a list of each line.
    public ArrayList<String> ParseFileList(String argFile)
			throws IOException {
    	Scanner s = new Scanner( new File(argFile));
    	ArrayList<String> list = new ArrayList<String>();
    	// read each line and add it to the list
    	while (s.hasNext()){
    		list.add(s.next());
    	}
    	s.close();
    	return list;
    }

	public static long hash(long startNode, long endNode) {
		return startNode * 1000000000l + endNode;
	}

	

}
