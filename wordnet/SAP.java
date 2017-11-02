import edu.princeton.cs.algs4.BreadthFirstDirectedPaths;
import edu.princeton.cs.algs4.Digraph;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

import java.util.Collections;

public class SAP {

    private Digraph digraph;

    // constructor takes a digraph (not necessarily a DAG)
    public SAP(Digraph G) {
        this.digraph = new Digraph(G);
    }

    // do unit testing of this class
    public static void main(String[] args) {
        In in = new In(args[0]);
        Digraph G = new Digraph(in);
        SAP sap = new SAP(G);
        while (!StdIn.isEmpty()) {
            int v = StdIn.readInt();
            int w = StdIn.readInt();
            int length = sap.length(v, w);
            int ancestor = sap.ancestor(v, w);
            StdOut.printf("length = %d, ancestor = %d\n", length, ancestor);
        }
    }

    // length of shortest ancestral path between v and w; -1 if no such path
    public int length(int v, int w) {
        return -1;
    }

    private int[] sap(Iterable<Integer> v, Iterable<Integer> w) {
        BreadthFirstDirectedPaths pathsV = new BreadthFirstDirectedPaths(this.digraph, v);
        BreadthFirstDirectedPaths pathsW = new BreadthFirstDirectedPaths(this.digraph, w);
        int numNode = this.digraph.V();
        int[] info = new int[2];
        info[0] = -1;
        info[1] = Integer.MAX_VALUE;
        for (int n = 0; n < numNode; n++) {
            if (pathsV.hasPathTo(n) && pathsW.hasPathTo(n)) {
                int distN = pathsV.distTo(n) + pathsW.distTo(n);
                if (distN < info[1]) {
                    info[0] = n;
                    info[1] = distN;
                }
            }
        }
        return info;
    }

    // a common ancestor of v and w that participates in a shortest ancestral path; -1 if no such path
    public int ancestor(int v, int w) {
        int[] sap = sap(Collections.singletonList(v), Collections.singletonList(w));
        return sap[0];
    }

    // length of shortest ancestral path between any vertex in v and any vertex in w; -1 if no such path
    public int length(Iterable<Integer> v, Iterable<Integer> w) {
        int[] sap = sap(v, w);
        if (sap[0] == -1) {
            return -1;
        } else {
            return sap[1];
        }
    }

    // a common ancestor that participates in shortest ancestral path; -1 if no such path
    public int ancestor(Iterable<Integer> v, Iterable<Integer> w) {
        int[] sap = sap(v, w);
        return sap[0];
    }
}
