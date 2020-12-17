import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;

public class Algorithm {
	// These hashmaps map strings for proteins like "P04355" to integers, and
	// vice versa
	// Overflow warning: These integers should never be multiplied, without
	// converting to longs.
	HashMap<String, Integer> mapToInt;
	HashMap<Integer, String> reverseMap;

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

	// Used internally in dijkstra's, and in CriticalPath for path printing.
	// Stores the dijkstra's path in symbolic tree form.
	int[] pathForward;
	int[] pathBackward;
	final long INF = Long.MAX_VALUE;

	// Stores the output to pass to the outputWriter class.
	StringBuilder pathOutput;
	StringBuilder edgeOutput;

    //maximum number of paths to write
    int maxk;

	public Algorithm(InputReader graph, int maxk) {
		// set up our variables
		mapToInt = graph.mapToInt;
		reverseMap = graph.reverseMap;

		edgeEndSet = graph.edgeEndSet;
		edgeStartSet = graph.edgeStartSet;
		weights = graph.weights;
		edges = graph.edges;
		reverseEdges = graph.reverseEdges;

		edgeCost = graph.edgeCost;
		pathOutput = new StringBuilder();
		edgeOutput = new StringBuilder();

        //maxk = graph.maxk;
        this.maxk = maxk;
	}

	public void run() {
		int start = 0; // super source = 0
		int end = 1; // super sink = 1
		int n = mapToInt.size(); // number of nodes
		int e = weights.size(); // number of edges

		// Dijkstra's from the start, and the end.
		// Stores the cost at from the start to all points.
		final double[] startFromAllNodes = dijkstraFront(n, start, edges);
		// Stores the cost from the end to all points.
		final double[] endFromAllNodes = dijkstra(n, end, reverseEdges);

		// Each edge has a shortest path; We will denote this as a
		// 'CriticalPath'.
		// A critical path does not contain only redundant edges. However, a
		// CriticalPath might be the same as another CriticalPath.
		TreeSet<CriticalPath> potentialPaths = new TreeSet<CriticalPath>();
		for (int a = 0; a < e; a++) {
			// CriticalPath has 4 parts; Source -> edgeStart -> edgeEnd -> Sink
			CriticalPath tempPath = new CriticalPath(edgeStartSet.get(a),
					edgeEndSet.get(a), startFromAllNodes, endFromAllNodes,
					edges);
			potentialPaths.add(tempPath);
		}

		ArrayList<CriticalPath> outputPaths = new ArrayList<CriticalPath>();
		edgeOutput.append("#tail	head	KSP index\n"); // this was in pathlinker
														// output so I do it as
														// well
		int count2 = 0;
		int countPath = 0;

		// This hashmap stores a blacklist of edges that have appeared on
		// earlier paths
		HashSet<Long> ReWriteThisWithEdgeClassLater = new HashSet<Long>();
		// go through all of the potential paths.
		while (!potentialPaths.isEmpty()) {
			count2++;
			CriticalPath get = potentialPaths.first();
			potentialPaths.remove(get);
			long startNode = get.startNode;
			long endNode = get.endNode;
			outputPaths.add(get);

			// output the edge
			if (startNode != 0 && endNode != 1 && countPath < maxk) {
			//if (startNode != 0 && endNode != 1) {
				String outputEdge = reverseMap.get((int) startNode) + "\t"
						+ reverseMap.get((int) endNode) + "\t" + count2;
				edgeOutput.append(outputEdge + "\n");
			}

			// get the shortest path that uses that edge.
			ArrayList<Integer> pathTemp = get.getPath();

			// figure out if this 'criticaledge' is new
			boolean newEdge = false;
			for (int b = 0; b < pathTemp.size() - 1; b++) {
				long hash = hash(pathTemp.get(b), pathTemp.get(b + 1));
				if (!ReWriteThisWithEdgeClassLater.contains(hash)) {
					newEdge = true;
					ReWriteThisWithEdgeClassLater.add(hash);
				}
			}

			// if it is new, output edge
			if (newEdge && countPath < maxk) {
			//if (newEdge) {
				countPath++;
				pathOutput.append(countPath + "\t"
						+ Math.pow(Math.E, -1 * get.totalCost) + "\t"
						+ getString(pathTemp, reverseMap) + "\n");
		        //System.out.println(countPath);
            }

		}

	}

	// Converts a list of node ID's to a list of node names.
	// Ex: [123, 4123] -> "P03422|Q02312"
	private static String getString(ArrayList<Integer> path,
			HashMap<Integer, String> reverseMap) {
		StringBuilder node_names = new StringBuilder();
		for (Integer a : path) {
			if (a != 0 && a != 1) {
                node_names.append(reverseMap.get(a) + "|");
            }
		}
        // remove the last '|'
        String output = node_names.toString();
        output = output.substring(0, output.length()-1);
		return output;
	}

	// Returns a unique 64-bit integer representation of an edge. (hash with no
	// collisions)
	// No two edges can share the same representation.
	private static long hash(long startNode, long endNode) {
		// As long as endNode is below 1 billion and above or equal to 0, this
		// function will never have a collision
		return startNode * 1000000000l + endNode;
	}

	class CriticalPath implements Comparable<CriticalPath> {
		// Start node of the edge defining the 'CriticalPath'
		int startNode;
		// End node of the edge defining the 'CriticalPath'
		int endNode;
		// Cost of the shortest path that uses the edge defined by startNode and
		// endNode.
		double totalCost;

		public CriticalPath(int startNode, int endNode,
				double[] startFromAllNodes, double[] endFromAllNodes,
				ArrayList<Edge>[] edges) {
			this.startNode = startNode;
			this.endNode = endNode;

			this.totalCost = startFromAllNodes[startNode]
					+ endFromAllNodes[endNode]
					+ edgeCost.get(hash(startNode, endNode));
		}

		// If the cost is the same, break ties arbitrarily (we use our hash)
		// compares by cost.
		@Override
		public int compareTo(CriticalPath arg0) {
			if (totalCost == arg0.totalCost) {

				return Long.compare(hash(startNode, endNode),
						hash(arg0.startNode, arg0.endNode));
			}
			return Double.compare(totalCost, arg0.totalCost);
		}

		// Computes the entire path, by traversing the dijkstra's algorithm
		// trees.
		public ArrayList<Integer> getPath() {
			ArrayList<Integer> output = new ArrayList<Integer>();
			ArrayList<Integer> firstHalf = new ArrayList<Integer>();
			int start = startNode;
			while (start != -1) {
				firstHalf.add(start);
				start = pathForward[start];
				if (start == -3) {
					return new ArrayList<Integer>();
				}
			}

			for (int a = firstHalf.size() - 1; a >= 0; a--) {
				output.add(firstHalf.get(a));
			}
			int end = endNode;
			while (end != -1) {
				output.add(end);
				end = pathBackward[end];
				if (end == -3) {
					return new ArrayList<Integer>();
				}
			}
			return output;
		}

		public String toString() {
			return startNode + " " + endNode + " " + totalCost + " "
					+ " Path: " + getString(this.getPath(), reverseMap);
		}
	}

	// /////////////////////////////////////////////////////////////////////////
	// Anything below here is a slightly modified dijkstras SSAD //
	// ////////////////////////////////////////////////////////////////////////
	static class VertexDist implements Comparable<VertexDist> {
		private static final double EPS = 1E-14;
		int vertex;
		double distance;

		public VertexDist(int vertex, double distance) {
			this.vertex = vertex;
			this.distance = distance;
		}

		public String toString() {
			return "(" + vertex + " " + distance + ")";
		}

		@Override
		public int compareTo(VertexDist arg0) {
			VertexDist temp = arg0;
			if (Math.abs(distance - temp.distance) < EPS) {
				return vertex - temp.vertex;
			}
			if (distance < temp.distance) {
				return -1;
			} else {
				return 1;
			}
		}

	}

	/*
	 * Extensively modified.
	 * 
	 * @author Godmar Back
	 */
	double[] dijkstra(int N, int start, ArrayList<Edge>[] edges) {
		final double[] dist = new double[N];
		Arrays.fill(dist, INF);
		dist[start] = 0;
		pathBackward = new int[N];
		Arrays.fill(pathBackward, -3);
		pathBackward[start] = -1;
		TreeSet<VertexDist> frontier = new TreeSet<VertexDist>();

		frontier.add(new VertexDist(start, 0));

		boolean[] optimal = new boolean[N];
		while (!frontier.isEmpty()) {
			VertexDist u2 = frontier.pollFirst();
			int u = u2.vertex;
			if (optimal[u])
				continue;

			optimal[u] = true;
			for (Edge e : edges[u]) {
				double uv = e.dist;
				int v = e.end;
				if (uv != INF) {
					if (dist[u] + uv < dist[v]) {
						dist[v] = dist[u] + uv;
						pathBackward[v] = u;
						frontier.add(new VertexDist(v, dist[v]));
					}
				}
			}
		}
		return dist;
		// 'dist' contains the shortest distance from start to all nodes
	}

	/*
	 * Extensively modified.
	 * 
	 * @author Godmar Back
	 */
	double[] dijkstraFront(int N, int start, ArrayList<Edge>[] edges) {
		final double[] dist = new double[N];
		pathForward = new int[N];
		Arrays.fill(pathForward, -3);

		Arrays.fill(dist, INF);
		dist[start] = 0;
		pathForward[start] = -1;
		TreeSet<VertexDist> frontier = new TreeSet<VertexDist>();

		frontier.add(new VertexDist(start, 0));
		boolean[] optimal = new boolean[N];
		while (!frontier.isEmpty()) {
			VertexDist u2 = frontier.pollFirst();
			int u = u2.vertex;
			if (optimal[u])
				continue;

			optimal[u] = true;
			for (Edge e : edges[u]) {
				double uv = e.dist;
				int v = e.end;
				if (uv != INF) {
					if (dist[u] + uv < dist[v]) {
						dist[v] = dist[u] + uv;
						pathForward[v] = u;
						frontier.add(new VertexDist(v, dist[v]));
					}
				}
			}
		}
		return dist;
		// 'dist' contains the shortest distance from start to all nodes
	}
}
