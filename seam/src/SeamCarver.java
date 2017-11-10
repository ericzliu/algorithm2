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
                    picture.setRGB(col++, row, current.getRGB(x, row));
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
        int dyg = cu.getGreen() - cd.getGreen();
        int dyb = cu.getBlue() - cd.getBlue();
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
        for (int col = 0; col < width; col++) {
            for (int row = 0; row < height; row++) {
                Pixel pixel = pixels[col][row];
                if (row + 1 < height) {
                    if (col - 1 >= 0) {
                        pixel.setLowerLeft(pixels[col - 1][row + 1]);
                    }
                    pixel.setLower(pixels[col][row + 1]);
                    if (col + 1 < width) {
                        pixel.setLowerRight(pixels[col + 1][row + 1]);
                    }
                }
                if (col + 1 < width) {
                    pixel.setRight(pixels[col + 1][row]);
                    if (row - 1 >= 0) {
                        pixel.setUpperRight(pixels[col + 1][row - 1]);
                    }
                }
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
        List<Pixel> layer = Arrays.asList(first);
        List<Pixel> lastLayer = layer;
        while (!layer.isEmpty()) {
            lastLayer = layer;
            layer = connectPixels.apply(lastLayer);
        }

        int dest = -1;
        AcyclicSP sp = new AcyclicSP(digraph, first.getId());
        double dist = Double.MAX_VALUE;
        for (Pixel pixel : lastLayer) {
            final int node = pixel.getId();
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
            if (!fromBuffer.isEmpty()) {
                Pixel first = getFirst(fromBuffer.get(0));
                if (first != null) {
                    toBuffer.add(first);
                }
                for (Pixel pixel : fromBuffer) {
                    Pixel middle = getMiddle(pixel);
                    if (middle == null) {
                        break;
                    }
                    toBuffer.add(middle);
                }
                Pixel last = getLast(fromBuffer.get(fromBuffer.size() - 1));
                if (last != null) {
                    toBuffer.add(last);
                }
            }
            for (Pixel pixel : fromBuffer) {
                addEdge(pixel, getFirst(pixel));
                addEdge(pixel, getMiddle(pixel));
                addEdge(pixel, getLast(pixel));
            }
            return toBuffer;
        }

        private void addEdge(Pixel pixel, Pixel p) {
            if (p != null) {
                digraph.addEdge(new DirectedEdge(pixel.getId(), p.getId(), p.getEnergy()));
            }
        }

        private Pixel getFirst(Pixel pixel) {
            switch (direction) {
                case Right:
                    return pixel.getUpperRight();
                case Down:
                    return pixel.getLowerLeft();
            }
            return null;
        }

        private Pixel getMiddle(Pixel pixel) {
            switch (direction) {
                case Right:
                    return pixel.getRight();
                case Down:
                    return pixel.getLower();
            }
            return null;
        }

        private Pixel getLast(Pixel pixel) {
            return pixel.getLowerRight();
        }

    }

    private static class Pixel {
        private int id;
        private int x; // column
        private int y; // row
        private double energy;
        private Pixel lowerLeft;
        private Pixel lower;
        private Pixel lowerRight;
        private Pixel right;
        private Pixel upperRight;

        Pixel(int id, int x, int y, double energy) {
            this.id = id;
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

        Pixel getLowerLeft() {
            return lowerLeft;
        }

        void setLowerLeft(Pixel lowerLeft) {
            this.lowerLeft = lowerLeft;
        }

        Pixel getLower() {
            return lower;
        }

        void setLower(Pixel lower) {
            this.lower = lower;
        }

        Pixel getLowerRight() {
            return lowerRight;
        }

        void setLowerRight(Pixel lowerRight) {
            this.lowerRight = lowerRight;
        }

        Pixel getRight() {
            return right;
        }

        void setRight(Pixel right) {
            this.right = right;
        }

        Pixel getUpperRight() {
            return upperRight;
        }

        void setUpperRight(Pixel upperRight) {
            this.upperRight = upperRight;
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
