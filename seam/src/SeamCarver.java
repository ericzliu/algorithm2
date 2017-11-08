import edu.princeton.cs.algs4.AcyclicSP;
import edu.princeton.cs.algs4.DirectedEdge;
import edu.princeton.cs.algs4.EdgeWeightedDigraph;
import edu.princeton.cs.algs4.Picture;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class SeamCarver {
    private Picture origin;
    private Picture current;

    private void assertNonNull(Object object) {
        if (null == object) {
            throw new IllegalArgumentException();
        }
    }
    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        assertNonNull(picture);
        this.origin = picture;
        this.current = new Picture(picture);
    }

    // current picture
    public Picture picture() {
        return this.current;
    }

    // width of current picture
    public int width() {
        return this.current.width();
    }

    // height of current picture
    public int height() {
        return this.current.height();
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        guardColumn(x);
        guardRow(y);
        if (x == 0 || y == 0 || x == this.width() - 1 || y == this.height() - 1) {
            return 1000;
        }
        int l = x - 1;
        int r = x + 1;
        Picture pic = picture();
        Color cl = pic.get(l, y);
        Color cr = pic.get(r, y);
        int u = y - 1;
        int d = y + 1;
        Color cu = pic.get(x, u);
        Color cd = pic.get(x, d);
        int dxr = cl.getRed() - cr.getRed();
        int dxg = cl.getGreen() - cr.getGreen();
        int dxb = cl.getBlue() - cr.getBlue();
        int dyr = cu.getRed() - cd.getRed();
        int dyg = cl.getGreen() - cd.getGreen();
        int dyb = cl.getBlue() - cd.getBlue();
        int sum = dxr * dxr + dxg * dxg + dxb * dxb + dyr * dyr + dyg * dyg + dyb * dyb;
        return Math.sqrt(sum);
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        return null;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        return null;
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (height() <= 1) {
            throw new IllegalArgumentException();
        }
        assertNonNull(seam);
        assertValidSeam(seam);
        int length = seam.length;
        if (length != width()) {
            throw new IllegalArgumentException();
        }
        for (int i = 0; i < length; i++) {
            guardRow(seam[i]);
        }
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (width() <= 1) {
            throw new IllegalArgumentException();
        }
        assertNonNull(seam);
        assertValidSeam(seam);
        int length = seam.length;
        if (length != height()) {
            throw new IllegalArgumentException();
        }
        for (int i = 0; i < length; i++) {
            guardColumn(seam[i]);
        }
        for (int x = 0; x < width(); x++) {

        }
    }

    private Seam verticalSeam(int x) {
        int height = height();
        int width = width();
        int id = 0;
        List<Pt> pts = new ArrayList<>();
        List<DirectedEdge> edges = new ArrayList<>();
        HashMap<Integer, Pt> lastLine = new HashMap<>();
        HashMap<Integer, Pt> line = new HashMap<>();
        for (int y = 0; y < height; y++) {
            int lo = Math.max(x - y, 0);
            int hi = Math.min(x + y, width - 1);
            line = new HashMap<>();
            for (int c = lo; c <= hi; c++) {
                Pt pt = new Pt(id++, c, y);
                pts.add(pt);
                line.put(c, pt);
                double energy = energy(c, y);
                if (lastLine.containsKey(c - 1)) {
                    edges.add(new DirectedEdge(lastLine.get(c - 1).id, pt.id, energy));
                }
                if (lastLine.containsKey(c)) {
                    edges.add(new DirectedEdge(lastLine.get(c).id, pt.id, energy));
                }
                if (lastLine.containsKey(c + 1)) {
                    edges.add(new DirectedEdge(lastLine.get(c + 1).id, pt.id, energy));
                }
            }
            lastLine = line;
        }
        EdgeWeightedDigraph digraph= new EdgeWeightedDigraph(pts.size());
        for (DirectedEdge edge: edges) {
            digraph.addEdge(edge);
        }
        double dist = Double.MAX_VALUE;
        int dest = -1;
        AcyclicSP sp = new AcyclicSP(digraph, 0);
        for (Pt pt : line.values()) {
            int ptId = pt.getId();
            if (sp.hasPathTo(ptId)) {
                double distTo = sp.distTo(ptId);
                if (distTo < dist) {
                    dist = distTo;
                    dest = ptId;
                }
            }
        }
        if (dest != -1) {
            int[] seam = new int[height];
            seam[0] = x;
            Iterable<DirectedEdge> path = sp.pathTo(dest);
            int loc = 1;
            for (DirectedEdge e : path) {
                int to = e.to();
                Pt pt = pts.get(to);
                seam[loc++] = pt.getId();
            }
            return new Seam(dist, seam);
        }
        throw new RuntimeException("No shortest path found.");
    }

    private void assertValidSeam(int[] seam) {
        int length = seam.length;
        for (int i = 1; i < length; i++) {
            if (Math.abs(seam[i] - seam[i - 1]) > 1) {
                throw new IllegalArgumentException();
            }
        }
    }

    private void guardColumn(int x) {
        if (x < 0 || x >= this.width()) {
            throw new IllegalArgumentException("Column out of range.");
        }
    }

    private void guardRow(int y) {
        if (y < 0 || y >= this.height()) {
            throw new IllegalArgumentException("Row out of range.");
        }
    }

    private static class Pt {
        private int id;
        private int x; // column
        private int y; // row

        public Pt(int id, int x, int y) {
            this.id = id;
            this.x = x;
            this.y = y;
        }

        public int getId() {
            return id;
        }

        public int getX() {
            return x;
        }

        public int getY() {
            return y;
        }
    }

    private static class Seam {
        private double dist;
        private int[] seam;

        public Seam(double dist, int[] seam) {
            this.dist = dist;
            this.seam = seam;
        }

        public double getDist() {
            return dist;
        }

        public int[] getSeam() {
            return seam;
        }
    }
}
