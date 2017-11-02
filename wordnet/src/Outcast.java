import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

public class Outcast {

    private WordNet wordNet;

    public Outcast(WordNet wordnet)         // constructor takes a WordNet object
    {
        this.wordNet = wordnet;
    }

    public static void main(String[] args) {
        WordNet wordnet = new WordNet(args[0], args[1]);
        Outcast outcast = new Outcast(wordnet);
        for (int t = 2; t < args.length; t++) {
            In in = new In(args[t]);
            String[] nouns = in.readAllStrings();
            StdOut.println(args[t] + ": " + outcast.outcast(nouns));
        }
    }

    private int sumDist(String noun, String[] nouns) {
        int dist = 0;
        for (String n : nouns) {
            if (!n.equals(noun)) {
                dist += wordNet.distance(n, noun);
            }
        }
        return dist;
    }

    public String outcast(String[] nouns)   // given an array of WordNet nouns, return an outcast
    {
        int maxDist = Integer.MIN_VALUE;
        String out = null;
        for (String n :
                nouns) {
            int dist = sumDist(n, nouns);
            if (dist > maxDist) {
                out = n;
                maxDist = dist;
            }
        }
        return out;
    }
}
