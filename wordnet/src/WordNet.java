import edu.princeton.cs.algs4.BreadthFirstDirectedPaths;
import edu.princeton.cs.algs4.Digraph;
import edu.princeton.cs.algs4.DirectedCycleX;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Queue;

public class WordNet {

    private HashMap<String, List<Integer>> wordMap = new HashMap<>();
    private List<String> synsets = new ArrayList<>();
    private Digraph digraph = null;

    private static boolean isRootedDAG(Digraph digraph) {
        if (null == digraph || digraph.V() < 1) {
            return false;
        }
        DirectedCycleX cycleX = new DirectedCycleX(digraph);
        if (cycleX.hasCycle()) {
            return false;
        }
        int numNode = digraph.V();
        int numRoot = 0;
        for (int v = 0; v < numNode; v++) {
            if (digraph.outdegree(v) == 0) {
                numRoot++;
            }
        }
        if (numRoot > 1) {
            return false;
        }
        return true;
    }

    // constructor takes the name of the two input files
    public WordNet(String synsets, String hypernyms) {
        if (null == synsets || null == hypernyms) {
            throw new IllegalArgumentException();
        }
        int nodeNum = 0;
        In synsetFile = new In(synsets);
        while (synsetFile.hasNextLine()) {
            nodeNum += 1;
            String line = synsetFile.readLine();
            String[] fields = line.split(",");
            int sid = Integer.parseInt(fields[0]);
            String wordsField = fields[1];
            String[] words = wordsField.split(" ");
            for (String word : words) {
                List<Integer> sids = wordMap.getOrDefault(word, new ArrayList<>());
                sids.add(sid);
                wordMap.put(word, sids);
            }
            this.synsets.add(wordsField);
        }

        digraph = new Digraph(nodeNum);
        In hypernymFile = new In(hypernyms);
        while (hypernymFile.hasNextLine()) {
            String line = hypernymFile.readLine();
            String[] sidStr = line.split(",");
            int length = sidStr.length;
            int[] sids = new int[length];
            for (int i = 0; i < length; i++) {
                sids[i] = Integer.parseInt(sidStr[i]);
            }
            int sid = sids[0];
            for (int i = 1; i < length; i++) {
                digraph.addEdge(sid, sids[i]);
            }
        }
        if (!isRootedDAG(digraph)) {
            throw new IllegalArgumentException();
        }
    }

    // do unit testing of this class
    public static void main(String[] args) {
        WordNet wordnet = new WordNet(args[0], args[1]);
        while (!StdIn.isEmpty()) {
            String v = StdIn.readString();
            String w = StdIn.readString();
            String ancestor = wordnet.sap(v, w);
            StdOut.printf("ancestor = %s\n", ancestor);
        }
    }

    // returns all WordNet nouns
    public Iterable<String> nouns() {
        return wordMap.keySet();
    }

    // is the word a WordNet noun?
    public boolean isNoun(String word) {
        if (null == word) {
            throw new IllegalArgumentException();
        }
        return wordMap.containsKey(word);
    }

    private void explore(Iterable<Integer> nodes, List<Integer> hyps, boolean[] marked, Queue<Integer> q) {
        for (Integer node : nodes) {
            if (!marked[node]) {
                marked[node] = true;
                q.add(node);
                hyps.add(node);
            }
        }
    }

    private List<Integer> getShortestAncestorInfo(String nounA, String nounB) {
        List<Integer> nodesA = wordMap.get(nounA);
        List<Integer> nodesB = wordMap.get(nounB);
        if (null == nodesA || null == nodesB) {
            throw new IllegalArgumentException();
        }
        BreadthFirstDirectedPaths pathsA = new BreadthFirstDirectedPaths(this.digraph, nodesA);
        BreadthFirstDirectedPaths pathsB = new BreadthFirstDirectedPaths(this.digraph, nodesB);
        List<Integer> hyps = new ArrayList<>();
        Queue<Integer> q = new ArrayDeque<>();
        boolean[] marked = new boolean[digraph.V()];
        explore(nodesA, hyps, marked, q);
        while (!q.isEmpty()) {
            Integer node = q.poll();
            Iterable<Integer> adj = digraph.adj(node);
            explore(adj, hyps, marked, q);
        }
        List<Integer> info = new ArrayList<>();
        info.add(-1);
        info.add(Integer.MAX_VALUE);
        for (Integer hyp : hyps) {
            if (pathsB.hasPathTo(hyp)) {
                int numHop = pathsA.distTo(hyp) + pathsB.distTo(hyp);
                if (numHop < info.get(1)) {
                    info.set(1, numHop);
                    info.set(0, hyp);
                }
            }
        }
        return info;
    }

    // distance between nounA and nounB (defined below)
    public int distance(String nounA, String nounB) {
        if (!isNoun(nounA) || !isNoun(nounB)) {
            throw new IllegalArgumentException();
        }
        return getShortestAncestorInfo(nounA, nounB).get(1);
    }


    // a synset (second field of synsets.txt) that is the common ancestor of nounA and nounB
    // in a shortest ancestral path (defined below)
    public String sap(String nounA, String nounB) {
        if (!isNoun(nounA) || !isNoun(nounB)) {
            throw new IllegalArgumentException();
        }
        List<Integer> info = getShortestAncestorInfo(nounA, nounB);
        return synsets.get(info.get(0));
    }
}
