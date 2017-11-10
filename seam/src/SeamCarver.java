/**
 *  Copyright Murex S.A.S., 2003-2017. All Rights Reserved.
 *
 *  This software program is proprietary and confidential to Murex S.A.S and its affiliates ("Murex") and, without limiting the generality of the foregoing reservation of rights, shall not be accessed, used, reproduced or distributed without the
 *  express prior written consent of Murex and subject to the applicable Murex licensing terms. Any modification or removal of this copyright notice is expressly prohibited.
 */
import java.awt.Color;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

import edu.princeton.cs.algs4.AcyclicSP;
import edu.princeton.cs.algs4.DirectedEdge;
import edu.princeton.cs.algs4.EdgeWeightedDigraph;
import edu.princeton.cs.algs4.Picture;


public class SeamCarver {

    //~ ----------------------------------------------------------------------------------------------------------------
    //~ Enums
    //~ ----------------------------------------------------------------------------------------------------------------

    enum Direction {
        Right, Down
    }

    //~ ----------------------------------------------------------------------------------------------------------------
    //~ Instance fields
    //~ ----------------------------------------------------------------------------------------------------------------

    private Picture current;
    private Pixel[][] pixels;

    //~ ----------------------------------------------------------------------------------------------------------------
    //~ Constructors
    //~ ----------------------------------------------------------------------------------------------------------------

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        assertNonNull(picture);
        this.current = new Picture(picture);
        this.pixels = createPixels(picture);
    }

    //~ ----------------------------------------------------------------------------------------------------------------
    //~ Methods
    //~ ----------------------------------------------------------------------------------------------------------------

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
        return energy(picture(), x, y);
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        final int height = height();
        double dist = Double.MAX_VALUE;
        Seam seam = null;
        for (int y = 0; y < height; y++) {
            Pixel from = pixels[0][y];
            final Seam tmp = getSeam(from, Direction.Right);
            if (tmp.dist < dist) {
                dist = tmp.dist;
                seam = tmp;
            }
        }
        if (null == seam) {
            throw new IllegalArgumentException();
        }
        return seam.getRowArray();
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        final int width = width();
        double dist = Double.MAX_VALUE;
        Seam seam = null;
        for (int x = 0; x < width; x++) {
            Pixel from = pixels[x][0];
            final Seam tmp = getSeam(from, Direction.Down);
            if (tmp.dist < dist) {
                dist = tmp.dist;
                seam = tmp;
            }
        }
        if (null == seam) {
            throw new IllegalArgumentException();
        }
        return seam.getColumnArray();
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        final int width = width();
        final int height = height();
        if (height <= 1) {
            throw new IllegalArgumentException();
        }
        assertNonNull(seam);
        assertValidSeam(seam);
        int length = seam.length;
        if (length != width) {
            throw new IllegalArgumentException();
        }
        for (int i = 0; i < length; i++) {
            guardRow(seam[i]);
        }
        final Picture picture = new Picture(width, height - 1);
        for (int col = 0; col < width; col++) {
            int row = 0;
            for (int y = 0; y < height; y++) {
                if (y != seam[col]) {
                    picture.setRGB(col, row++, current.get(col, y).getRGB());
                }
            }
        }
        this.current = picture;
        this.pixels = createPixels(picture);
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        final int width = width();
        final int height = height();
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
        final Picture picture = new Picture(width - 1, height);
        for (int row = 0; row < height; row++) {
            int col = 0;
            for (int x = 0; x < width; x++) {
                if (x != seam[row]) {
                    picture.setRGB(col++, row, current.getRGB(col, x));
                }
            }
        }
        this.current = picture;
        this.pixels = createPixels(picture);
    }

    private static double energy(Picture pic, int x, int y) {
        if ((x == 0) || (y == 0) || (x == (pic.width() - 1)) || (y == (pic.height() - 1))) {
            return 1000;
        }
        int l = x - 1;
        int r = x + 1;
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
        int sum = (dxr * dxr) + (dxg * dxg) + (dxb * dxb) + (dyr * dyr) + (dyg * dyg) + (dyb * dyb);
        return Math.sqrt(sum);
    }

    private static Pixel[][] createPixels(Picture picture) {
        final int width = picture.width();
        final int height = picture.height();
        Pixel[][] pixels = new Pixel[width][height];
        for (int col = 0; col < width; col++) {
            for (int row = 0; row < height; row++) {
                pixels[col][row] = new Pixel((row * width) + col, col, row, energy(picture, col, row));
            }
        }
        return pixels;
    }

    private void assertNonNull(Object object) {
        if (null == object) {
            throw new IllegalArgumentException();
        }
    }

    private Seam getSeam(Pixel first, Direction direction) {
        final int height = height();
        final int width = width();
        EdgeWeightedDigraph digraph = new EdgeWeightedDigraph(width * height);
        final ConnectPixels connectPixels = new ConnectPixels(digraph, direction, pixels);
        List<Pixel> buffer = Arrays.asList(first);
        List<Pixel> boundary = buffer;
        while (!buffer.isEmpty()) {
            boundary = buffer;
            buffer = connectPixels.apply(buffer);
        }
        int dest = -1;
        AcyclicSP sp = new AcyclicSP(digraph, getNodeId(first));
        double dist = Double.MAX_VALUE;
        for (Pixel pixel : boundary) {
            final int node = getNodeId(pixel);
            if (sp.hasPathTo(node) && (sp.distTo(node) < dist)) {
                dist = sp.distTo(node);
                dest = node;
            }
        }
        final Iterable<DirectedEdge> edges = sp.pathTo(dest);
        final List<Pixel> path = new ArrayList<>();
        path.add(first);
        for (DirectedEdge edge : edges) {
            int row = edge.to() / width;
            int col = edge.to() % width;
            path.add(pixels[col][row]);
        }
        return new Seam(dist, path);
    }

    private int getNodeId(Pixel pixel) {
        return (pixel.getY() * width()) + pixel.getX();
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
        if ((x < 0) || (x >= this.width())) {
            throw new IllegalArgumentException("Column out of range.");
        }
    }

    private void guardRow(int y) {
        if ((y < 0) || (y >= this.height())) {
            throw new IllegalArgumentException("Row out of range.");
        }
    }

    //~ ----------------------------------------------------------------------------------------------------------------
    //~ Nested Classes
    //~ ----------------------------------------------------------------------------------------------------------------

    private class ConnectPixels implements Function<List<Pixel>, List<Pixel>> {

        private Direction direction;
        private EdgeWeightedDigraph digraph;
        private Pixel[][] pixels;

        ConnectPixels(EdgeWeightedDigraph digraph, Direction direction, Pixel[][] pixels) {
            this.direction = direction;
            this.digraph = digraph;
            this.pixels = pixels;
        }

        @Override
        public List<Pixel> apply(List<Pixel> fromBuffer) {
            final ArrayList<Pixel> toBuffer = new ArrayList<>();
            switch (this.direction) {

            case Down:
                int row = fromBuffer.get(0).getY() + 1;
                if (row < height()) {
                    int lo = Math.max(0, fromBuffer.get(0).getX() - 1);
                    int hi = Math.min(width() - 1, fromBuffer.get(fromBuffer.size() - 1).getX() + 1);
                    for (int col = lo; col <= hi; col++) {
                        toBuffer.add(pixels[col][row]);
                    }
                    for (Pixel from : fromBuffer) {
                        safeAddEdge(toBuffer, from, (from.getX() - 1) - lo);
                        safeAddEdge(toBuffer, from, from.getX() - lo);
                        safeAddEdge(toBuffer, from, (from.getX() + 1) - lo);
                    }
                }
                break;

            case Right:
                int col = fromBuffer.get(0).getX() + 1;
                if (col < width()) {
                    int lo = Math.max(0, fromBuffer.get(0).getY() - 1);
                    int hi = Math.min(height() - 1, fromBuffer.get(fromBuffer.size() - 1).getY() + 1);
                    for (int r = lo; r <= hi; r++) {
                        toBuffer.add(pixels[col][r]);
                    }
                    for (Pixel from : fromBuffer) {
                        safeAddEdge(toBuffer, from, (from.getY() - 1) - lo);
                        safeAddEdge(toBuffer, from, from.getY() - lo);
                        safeAddEdge(toBuffer, from, (from.getY() + 1) - lo);
                    }
                }
                break;
            }
            return toBuffer;
        }

        private void safeAddEdge(ArrayList<Pixel> buffer, Pixel from, int loc) {
            if ((loc >= 0) && (loc < buffer.size())) {
                final Pixel to = buffer.get(loc);
                digraph.addEdge(new DirectedEdge(from.getId(), to.getId(), to.getEnergy()));
            }
        }
    }

    private static class Pixel {
        private int id;
        private int x; // column
        private int y; // row
        private double energy;

        Pixel(int id, int x, int y, double energy) {
            this.x = x;
            this.y = y;
            this.energy = energy;
        }

        int getId() {
            return id;
        }

        double getEnergy() {
            return energy;
        }

        int getY() {
            return y;
        }

        int getX() {
            return x;
        }
    }

    private static class Seam {
        private double dist;
        private List<Pixel> nodes;

        Seam(double dist, List<Pixel> nodes) {
            this.dist = dist;
            this.nodes = nodes;
        }

        int[] getColumnArray() {
            final int length = nodes.size();
            int[] result = new int[length];
            for (int i = 0; i < length; i++) {
                result[i] = nodes.get(i).getX();
            }
            return result;
        }

        int[] getRowArray() {
            final int length = nodes.size();
            int[] result = new int[length];
            for (int i = 0; i < length; i++) {
                result[i] = nodes.get(i).getY();
            }
            return result;
        }

    }

}
