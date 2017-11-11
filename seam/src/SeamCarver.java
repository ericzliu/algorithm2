import java.awt.Color;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.princeton.cs.algs4.DirectedEdge;
import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Stack;


public class SeamCarver {

    private enum Direction {
        Right, Down
    }

    private Color[][] colorMatrix;
    private double[][] energyMatrix;
    private int width;
    private int height;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        assertNonNull(picture);
        updateInternals(getColorMatrix(picture));
    }

    private void updateInternals(Color[][] colors) {
        colorMatrix = colors;
        energyMatrix = getEnergyMatrix(colors);
        width = colorMatrix.length;
        height = colorMatrix[0].length;
    }

    //~ ----------------------------------------------------------------------------------------------------------------
    //~ Methods
    //~ ----------------------------------------------------------------------------------------------------------------

    // current picture
    public Picture picture() {
        Picture picture = new Picture(width, height);
        for (int c = 0; c < width; c++) {
            for (int r = 0; r < height; r++) {
                picture.set(c, r, colorMatrix[c][r]);
            }
        }
        return picture;
    }

    // width of current picture
    public int width() {
        return width;
    }

    // height of current picture
    public int height() {
        return height;
    }

    // getEnergy of pixel at column x and row y
    public double energy(int x, int y) {
        guardColumn(x);
        guardRow(y);
        return energyMatrix[x][y];
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        double dist = Double.MAX_VALUE;
        Seam seam = null;
        for (int y = 0; y < height; y++) {
            final Seam tmp = getSeam(getId(0, y), Direction.Right);
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
        double dist = Double.MAX_VALUE;
        Seam seam = null;
        for (int x = 0; x < width; x++) {
            final Seam tmp = getSeam(getId(x, 0), Direction.Down);
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
        if (height <= 1) {
            throw new IllegalArgumentException();
        }
        assertNonNull(seam);
        assertValidSeam(seam);
        int length = seam.length;
        if (length != width) {
            throw new IllegalArgumentException();
        }
        for (int row : seam) {
            guardRow(row);
        }
        Color[][] pixels = new Color[width][height - 1];
        for (int col = 0; col < width; col++) {
            int row = 0;
            for (int y = 0; y < height; y++) {
                if (y != seam[col]) {
                    pixels[col][row++] = colorMatrix[col][y];
                }
            }
        }
        updateInternals(pixels);
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
        for (int col : seam) {
            guardColumn(col);
        }
        Color[][] pixels = new Color[width - 1][height];
        for (int row = 0; row < height; row++) {
            int col = 0;
            for (int x = 0; x < width; x++) {
                if (x != seam[row]) {
                    pixels[col++][row] = colorMatrix[x][row];
                }
            }
        }
        updateInternals(pixels);

    }

    private static double getEnergy(Color[][] pic, int x, int y) {
        int width = pic.length;
        int height = pic[0].length;
        if ((x == 0) || (y == 0) || (x == (width - 1)) || (y == (height - 1))) {
            return 1000;
        }
        int g = x - 1;
        int r = x + 1;
        Color cl = pic[g][y];
        Color cr = pic[r][y];
        int u = y - 1;
        int d = y + 1;
        Color cu = pic[x][u];
        Color cd = pic[x][d];
        int dxr = cl.getRed() - cr.getRed();
        int dxg = cl.getGreen() - cr.getGreen();
        int dxb = cl.getBlue() - cr.getBlue();
        int dyr = cu.getRed() - cd.getRed();
        int dyg = cu.getGreen() - cd.getGreen();
        int dyb = cu.getBlue() - cd.getBlue();
        int sum = (dxr * dxr) + (dxg * dxg) + (dxb * dxb) + (dyr * dyr) + (dyg * dyg) + (dyb * dyb);
        return Math.sqrt(sum);
    }

    private static Color[][] getColorMatrix(Picture picture) {
        final int width = picture.width();
        final int height = picture.height();
        Color[][] colors = new Color[width][height];
        for (int col = 0; col < width; col++) {
            for (int row = 0; row < height; row++) {
                colors[col][row] = picture.get(col, row);
            }
        }
        return colors;
    }

    private static double[][] getEnergyMatrix(Color[][] colors) {
        int width = colors.length;
        int height = colors[0].length;
        double[][] matrix = new double[width][height];
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                matrix[x][y] = getEnergy(colors, x, y);
            }
        }
        return matrix;
    }

    private void assertNonNull(Object object) {
        if (null == object) {
            throw new IllegalArgumentException();
        }
    }

    private Seam getSeam(int first, Direction direction) {
        final ConnectPixels connect = new ConnectPixels(width, height, direction, first);
        List<java.lang.Integer> layer = Collections.singletonList(first);
        List<java.lang.Integer> lastLayer = layer;
        while (!layer.isEmpty()) {
            lastLayer = layer;
            layer = connect.apply(lastLayer);
        }
        int dest = -1;
        double dist = Double.MAX_VALUE;
        for (java.lang.Integer pixel : lastLayer) {
            if (connect.hasPathTo(pixel) && (connect.distTo(pixel) < dist)) {
                dist = connect.distTo(pixel);
                dest = pixel;
            }
        }
        final Iterable<DirectedEdge> edges = connect.pathTo(dest);
        final List<java.lang.Integer> path = new ArrayList<>();
        path.add(first);
        for (DirectedEdge edge : edges) {
            path.add(edge.to());
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

    private class ConnectPixels {

        private Direction direction;
        private double[] distTo;
        private DirectedEdge[] edgeTo;
        private int width;
        private int height;

        ConnectPixels(int width, int height, Direction direction, int from) {
            this.width = width;
            this.height = height;
            this.direction = direction;
            final int count = width * height;
            this.distTo = new double[count];
            this.edgeTo = new DirectedEdge[count];
            for (int i = 0; i < count; i++) {
                distTo[i] = Double.MAX_VALUE;
                edgeTo[i] = null;
            }
            distTo[from] = energyMatrix[getX(from)][getY(from)];
        }

        boolean hasPathTo(int node) {
            return this.distTo[node] != Double.MAX_VALUE;
        }

        double distTo(int node) {
            return this.distTo[node];
        }

        Iterable<DirectedEdge> pathTo(int v) {
            Stack<DirectedEdge> path = new Stack<>();
            for (DirectedEdge e = edgeTo[v]; e != null; e = edgeTo[e.from()]) {
                path.push(e);
            }
            return path;
        }

        List<java.lang.Integer> apply(List<java.lang.Integer> fromBuffer) {
            final ArrayList<java.lang.Integer> toBuffer = new ArrayList<>();
            if (!fromBuffer.isEmpty()) {
                int first = getFirst(fromBuffer.get(0));
                if (first >= 0) {
                    toBuffer.add(first);
                }
                for (java.lang.Integer pixel : fromBuffer) {
                    int middle = getMiddle(pixel);
                    if (middle < 0) {
                        break;
                    }
                    toBuffer.add(middle);
                }
                java.lang.Integer end = fromBuffer.get(fromBuffer.size() - 1);
                int last = getLast(end);
                if (last >= 0) {
                    toBuffer.add(last);
                }
            }
            for (java.lang.Integer pixel : fromBuffer) {
                relax(pixel, getFirst(pixel));
                relax(pixel, getMiddle(pixel));
                relax(pixel, getLast(pixel));
            }
            return toBuffer;
        }

        private void relax(int from, int to) {
            if (to >= 0) {
                double energy = energyMatrix[getX(to)][getY(to)];
                double newDist = distTo[from] + energy;
                if (newDist < distTo[to]) {
                    distTo[to] = newDist;
                    edgeTo[to] = new DirectedEdge(from, to, energy);
                }
            }
        }

        private boolean isOutside(int x, int y) {
            return (x < 0 || x >= width || y < 0 || y >= height);
        }

        private int getFirst(int id) {
            int x = getX(id);
            int y = getY(id);
            switch (direction) {
                case Right:
                    if (!isOutside(x + 1, y - 1)) {
                        return getId(x + 1, y - 1);
                    }
                    break;
                case Down:
                    if (!isOutside(x - 1, y + 1)) {
                        return getId(x - 1, y + 1);
                    }
                    break;
            }
            return -1;
        }

        private int getMiddle(int id) {
            int x = getX(id);
            int y = getY(id);
            switch (direction) {
                case Right:
                    if (!isOutside(x + 1, y)) {
                        return getId(x + 1, y);
                    }
                    break;
                case Down:
                    if (!isOutside(x, y + 1)) {
                        return getId(x, y + 1);
                    }
                    break;
            }
            return -1;
        }

        private int getLast(int id) {
            int x = getX(id);
            int y = getY(id);
            if (!isOutside(x + 1, y + 1)) {
                return getId(x + 1, y + 1);
            }
            return -1;
        }
    }


    private int getId(int x, int y) {
        return y * width + x;
    }

    private int getX(int id) {
        return id % width;
    }

    private int getY(int id) {
        return id / width;
    }

    private class Seam {
        private double dist;
        private List<Integer> nodes;

        Seam(double dist, List<Integer> nodes) {
            this.dist = dist;
            this.nodes = nodes;
        }

        int[] getColumnArray() {
            final int length = nodes.size();
            int[] result = new int[length];
            for (int i = 0; i < length; i++) {
                result[i] = getX(nodes.get(i));
            }
            return result;
        }

        int[] getRowArray() {
            final int length = nodes.size();
            int[] result = new int[length];
            for (int i = 0; i < length; i++) {
                result[i] = getY(nodes.get(i));
            }
            return result;
        }

    }

}
