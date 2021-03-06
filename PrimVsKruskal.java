/* PrimVsKruskal.java
   CSC 226 - Summer 2018
   Assignment 2 - Prim MST versus Kruskal MST Template

   The file includes the "import edu.princeton.cs.algs4.*;" so that yo can use
   any of the code in the algs4.jar file. You should be able to compile your program
   with the command

	javac -cp .;algs4.jar PrimVsKruskal.java

   To conveniently test the algorithm with a large input, create a text file
   containing a test graphs (in the format described below) and run
   the program with

	java -cp .;algs4.jar PrimVsKruskal file.txt

   where file.txt is replaced by the name of the text file.

   The input consists of a graph (as an adjacency matrix) in the following format:

    <number of vertices>
	<adjacency matrix row 1>
	...
	<adjacency matrix row n>

   Entry G[i][j] >= 0.0 of the adjacency matrix gives the weight (as type double) of the edge from
   vertex i to vertex j (if G[i][j] is 0.0, then the edge does not exist).
   Note that since the graph is undirected, it is assumed that G[i][j]
   is always equal to G[j][i].


   R. Little - 06/22/2018
*/

 import edu.princeton.cs.algs4.*;
 import java.util.Scanner;
 import java.io.File;

//Do not change the name of the PrimVsKruskal class
public class PrimVsKruskal{
	/* PrimVsKruskal(G)
		Given an adjacency matrix for connected graph G, with no self-loops or parallel edges,
		determine if the minimum spanning tree of G found by Prim's algorithm is equal to
		the minimum spanning tree of G found by Kruskal's algorithm.

		If G[i][j] == 0.0, there is no edge between vertex i and vertex j
		If G[i][j] > 0.0, there is an edge between vertices i and j, and the
		value of G[i][j] gives the weight of the edge.
		No entries of G will be negative.
	*/

    /**
     * Prim's Algorithm
     */
    // public class PrimMST {
    //     private Edge[] edgeTo; // shortest edge from tree to vertex
    //     private double[] distTo; // distTo[w] = edgeTo[w].weight()
    //     private boolean[] marked; // true if v in mst
    //     private IndexMinPQ<Edge> pq; // eligible crossing edges
    //
        // public PrimMST(WeightedGraph G) {
        //     edgeTo = new Edge[G.V()];
        //     distTo = new double[G.V()];
        //     marked = new boolean[G.V()];
        //
        //     for(int v = 0; v < G.V(); v++) {
        //         distTo[v] = Double.POSITIVE_INFINITY;
        //     }
        //
        //     pq = new IndexMinPQ<Double>(G.V());
        //     distTo[0[ = 0.0;
        //     pq.insert(0, 0.0);
            // while(!pq.isEmpty()) {
            //     visit(G, pq.delMin());
            // }
        //
        // }
    //
    //     private void visit(WeightedGraph G, int v) {
    //         marked[v] = true;
    //         for (Edge e : G.adj(v)) {
    //             int w = e.other(v);
    //             if (marked[w]) continue;
    //             if (e.weight() < distTo[w]) {
    //                 edgeTo[w] = e;
    //                 distTo[w] = e.weight();
    //                 if (pq.contains(w)) pq.changeKey(w, distTo[w]);
    //                 else pq.insert(w, distTo[w]);
    //             }
    //         }
    //     }
    // }

    static boolean parallelPrimVsKruskal(double[][] G) {
        int N = G.length;
        boolean[][] key = new boolean[N][N];
        EdgeWeightedGraph EG = buildGraph(G);
        UF KPuf = new UF(N);

        /**
         * Prim Vars
         */
        Edge[] edgeTo = new Edge[N]; // shortest edge from tree to vertex
        double[] distTo = new double[N]; // distTo[w] = edgeTo[w].weight()
        boolean[] marked = new boolean[N]; // true if v in mst
        IndexMinPQ<Double> Ppq; // eligible crossing edges
        Edge ep;
        double weight;

        // Initialize Prim's disTo array
        for(int v = 0; v < N; v++) {
            distTo[v] = Double.POSITIVE_INFINITY;
        }

        // Initialize Prim's pq
        Ppq = new IndexMinPQ<Double>(N);
        distTo[0] = 0.0;
        Ppq.insert(0, 0.0);

        /**
         * Kruskal Vars
         */
        Queue<Edge> mst = new Queue<Edge>();
        MinPQ<Edge> Kpq = new MinPQ<Edge>();
        UF uf;

        // Initialize kruskal pq
        for(Edge e : EG.edges()) {
            Kpq.insert(e);
        }

        // Initialize Disjoint Set datastructure
        uf = new UF(N);

        /**
         * Concurrent Prim's and Kruskal's
         */
         String KEdges = "";
         String PEdges = "";
         while( !Ppq.isEmpty() && ( !Kpq.isEmpty() && mst.size() < N-1 ) ) {
             // prim operation
             if(!Ppq.isEmpty()) {
                 // visit the top of the heap
                int vp = Ppq.delMin();
                marked[vp] = true;

                for (Edge e : EG.adj(vp)) {
                    int wp = e.other(vp);
                    if (marked[wp]) continue;
                    if (e.weight() < distTo[wp]) {
                        edgeTo[wp] = e;
                        distTo[wp] = e.weight();
                        if (Ppq.contains(wp)) Ppq.changeKey(wp, distTo[wp]);
                        else Ppq.insert(wp, distTo[wp]);
                    }

                }

                // kruskal
                if(!Kpq.isEmpty() && mst.size() < N-1) {
                    Edge ek = Kpq.delMin();
                    int vk = ek.either();
                    int wk = ek.other(vk);
                    if(uf.connected(vk, wk)) continue;
                    uf.union(vk, wk);
                    mst.enqueue(ek);
                    StdOut.printf("Kruskal <- (%s)\n", ek);
                    KEdges += "("+ek.toString()+"), ";
                    if(!key[vk][wk] || !key[wk][vk]) {
                        if(KPuf.connected(vk, wk)) {
                            StdOut.printf("Divergent Edge: %d-%d\n", vk, wk);
                            return false;
                        }
                    }
                    KPuf.union(vk, wk);
                    key[vk][wk] = true;
                    key[wk][vk] = true;
                }

                if(!Ppq.isEmpty()) {
                    Edge primEdge = edgeTo[Ppq.minIndex()];
                    StdOut.printf("Prim <- (%s)\n", primEdge);
                    vp = primEdge.either();
                    int wp = primEdge.other(vp);
                    if(!key[vp][wp] || !key[wp][vp]) {
                        if(KPuf.connected(vp, wp)) {
                            StdOut.printf("Divergent Edge: %d-%d\n", vp, wp);
                            return false;
                        }
                    }
                    KPuf.union(vp, wp);
                    key[vp][wp] = true;
                    key[wp][vp] = true;
                }

             }
         }

         StdOut.printf("Kruskal's: \n %s\n", KEdges);
         StdOut.printf("Prims's:\n");
         StdOut.println();
         for(Edge e : edgeTo) {
             if(e != null && e.weight() != 0.0)
                StdOut.printf("(%s), ", e.toString());
         }

         StdOut.println();
         return true;
    }

    public static boolean cpPrimOpToKruskal(UF uf, boolean[][] key, Edge min) {
        int vp = min.either();
        int wp = min.other(vp);
        // prim added an edge that kruskals hasn't
        if(!key[vp][wp]) {
            // are the ednpoints of prim's edge already in kruskal's?
            if(uf.connected(vp, wp)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Build a MinPQ of edges from an adjacency matrix
     * @param  G adjacency matrix
     * @return pq  MinPQ
     */
    public static MinPQ<Edge> buildEdgeMinPQ(double[][] G) {
        int V = G.length;
        Edge e;
        double weight;
        MinPQ<Edge> pq = new MinPQ<Edge>();
        for(int i = 0; i < V; i++) {
            for(int j = i; j < V; j++) {
                weight = G[i][j];
                if(weight != 0) {
                    e = new Edge(i, j, weight);
                    pq.insert(e);
                }
            }
        }
        return pq;
    }

	static boolean PrimVsKruskal(double[][] G){
		int n = G.length;

		/* Build the MST by Prim's and the MST by Kruskal's */
		/* (You may add extra methods if necessary) */

		/* ... Your code here ... */
        EdgeWeightedGraph EG = buildGraph(G);
        outputEdgeData(EG);

        StdOut.println("Concurrent result: ");
        if(parallelPrimVsKruskal(G)) {
            StdOut.println("true");
        } else {
            StdOut.println("false");
        }


		/* Determine if the MST by Prim equals the MST by Kruskal */
		boolean pvk = true;
		/* ... Your code here ... */
        pvk = naiveSolution(EG, n);

		return pvk;
	}

    /**
     * @author Avery K.
     *
     * Prints out edges from an Iterable edgelist
     * @param Edges Iterable edgelist
     */
    public static void outputEdgeData(EdgeWeightedGraph EG) {
        StdOut.printf("Edge List: ");
        printEdges(EG.edges());

		PrimMST PMst = new PrimMST(EG);
        KruskalMST KMST = new KruskalMST(EG);

        StdOut.printf("Prim Edge List: ");
        printEdges(PMst.edges());

        StdOut.printf("Kruksal's Edge List: ");
        printEdges(KMST.edges());
    }

    /**
     * @author Avery K.
     *
     * Prints out edges from an Iterable edgelist
     * @param Edges Iterable edgelist
     */
    public static void printEdges(Iterable<Edge> Edges) {
        StdOut.println();
        for(Edge e : Edges) {
            StdOut.printf("(%s), ", e.toString());
        }

        StdOut.println();
    }

    /**
     * [naiveSolution description]
     * @param  EG [description]
     * @param  N  number of nodes in the graph
     * @return    [description]
     */
    public static boolean naiveSolution(EdgeWeightedGraph EG, int N) {
        boolean[][] key = new boolean[N][N];
        int v, w;

        PrimMST PMst = new PrimMST(EG);
        KruskalMST KMst = new KruskalMST(EG);
        Iterable<Edge> PMstEdges = PMst.edges();
        Iterable<Edge> KMstEdges = KMst.edges();

        for(Edge e : PMst.edges()) {
            v = e.either();
            w = e.other(v);
            key[v][w] = true;
        }

        for(Edge e : KMst.edges()) {
            v = e.either();
            w = e.other(v);
            if(!key[v][w]) {
                return false;
            }
        }
        return true;
    }

    /**
     * [generateRandomGraph description]
     * @param  N [description]
     * @return   [description]
     */
    public static double[][] generateRandomGraph(int N) {
        double[][] G = new double[N][N];
        double r;
        for(int i = 0; i < N; i++) {
            for(int j = i; j < N; j++) {
                if(j != i) {
                    if(G[i][j] == 0.0) {
                        r = StdRandom.uniform(N-1);
                        G[i][j] = r;
                    }
                }
            }
        }

        return G;
    }

    /**
     * @author Avery K.
     *
     * Builds an EdgeWeightedGraph from an adjacency matrix as input
     * @param  G adjacency matrix
     * @return EG  EdgeWeightedGraph
     */
    public static EdgeWeightedGraph buildGraph(double[][] G) {
        int V = G.length;
        Edge e;
        double weight;

        EdgeWeightedGraph EG = new EdgeWeightedGraph(V);
        for(int i = 0; i < V; i++) {
            for(int j = i; j < V; j++) {
                weight = G[i][j];
                if(weight != 0) {
                    e = new Edge(i, j, weight);
                    EG.addEdge(e);
                }
            }
        }
        return EG;
    }

	/* main()
	   Contains code to test the PrimVsKruskal function. You may modify the
	   testing code if needed, but nothing in this function will be considered
	   during marking, and the testing process used for marking will not
	   execute any of the code below.
	*/
   public static void main(String[] args) {
		Scanner s;
		if (args.length > 0){
			try{
				s = new Scanner(new File(args[0]));
			} catch(java.io.FileNotFoundException e){
				System.out.printf("Unable to open %s\n",args[0]);
				return;
			}
			System.out.printf("Reading input values from %s.\n",args[0]);
		}else{
			s = new Scanner(System.in);
			System.out.printf("Reading input values from stdin.\n");
		}

		int n = s.nextInt();
		double[][] G = new double[n][n];
		int valuesRead = 0;
		for (int i = 0; i < n && s.hasNextDouble(); i++){
			for (int j = 0; j < n && s.hasNextDouble(); j++){
				G[i][j] = s.nextDouble();
				if (i == j && G[i][j] != 0.0) {
					System.out.printf("Adjacency matrix contains self-loops.\n");
					return;
				}
				if (G[i][j] < 0.0) {
					System.out.printf("Adjacency matrix contains negative values.\n");
					return;
				}
				if (j < i && G[i][j] != G[j][i]) {
					System.out.printf("Adjacency matrix is not symmetric.\n");
					return;
				}
				valuesRead++;
			}
		}

		if (valuesRead < n*n){
			System.out.printf("Adjacency matrix for the graph contains too few values.\n");
			return;
		}
        // G = generateRandomGraph(20);
        boolean pvk = PrimVsKruskal(G);
        System.out.printf("Does Prim MST = Kruskal MST? %b\n", pvk);
    }
}
