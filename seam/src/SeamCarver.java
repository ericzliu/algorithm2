import java.util.ArrayList;
import java.util.List;

import edu.princeton.cs.algs4.DirectedEdge;
import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Stack;


public class SeamCarver {

    private enum Direction {
        Right, Down
    }

    private int[][] colorMatrix;
    private double[][] energyMatrix;
    private int width;
    private int height;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        assertNonNull(picture);
        updateInternals(getColorMatrix(picture));
    }

    private void updateInternals(int[][] colors) {
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
                picture.setRGB(c, r, colorMatrix[c][r]);
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
        Seam seam = getSeam(Direction.Right);
        return seam.getRowArray();
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        Seam seam = getSeam(Direction.Down);
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
        int[][] pixels = new int[width][height - 1];
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
        int[][] pixels = new int[width - 1][height];
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

    private static int getRed(int color) {
        return (color >> 16) & 0xFF;
    }

    private static int getGreen(int color) {
        return (color >> 8) & 0xFF;
    }

    private static int getBlue(int color) {
        return (color >> 0) & 0xFF;
    }

    private static double getEnergy(int[][] pic, int x, int y) {
        int width = pic.length;
        int height = pic[0].length;
        if ((x == 0) || (y == 0) || (x == (width - 1)) || (y == (height - 1))) {
            return 1000;
        }
        int g = x - 1;
        int r = x + 1;
        int cl = pic[g][y];
        int cr = pic[r][y];
        int u = y - 1;
        int d = y + 1;
        int cu = pic[x][u];
        int cd = pic[x][d];
        int dxr = getRed(cl) - getRed(cr);
        int dxg = getGreen(cl) - getGreen(cr);
        int dxb = getBlue(cl) - getBlue(cr);
        int dyr = getRed(cu) - getRed(cd);
        int dyg = getGreen(cu) - getGreen(cd);
        int dyb = getBlue(cu) - getBlue(cd);
        int sum = (dxr * dxr) + (dxg * dxg) + (dxb * dxb) + (dyr * dyr) + (dyg * dyg) + (dyb * dyb);
        return Math.sqrt(sum);
    }

    private static int[][] getColorMatrix(Picture picture) {
        final int width = picture.width();
        final int height = picture.height();
        int[][] colors = new int[width][height];
        for (int col = 0; col < width; col++) {
            for (int row = 0; row < height; row++) {
                colors[col][row] = picture.getRGB(col, row);
            }
        }
        return colors;
    }

    private static double[][] getEnergyMatrix(int[][] colors) {
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

    private Seam getSeam(Direction direction) {
        int[] dests;
        if (direction == Direction.Right) {
            dests = new int[height];
            for (int i = 0; i < height; i++) {
                dests[i] = getId(width - 1, i);
            }
        }
        else {
            dests = new int[width];
            for (int i = 0; i < width; i++) {
                dests[i] = getId(i, height - 1);
            }
        }

        final ConnectPixels connect = new ConnectPixels(width, height, direction);
        int dest = -1;
        double dist = Double.MAX_VALUE;
        for (int pixel : dests) {
            if (connect.hasPathTo(pixel) && (connect.distTo(pixel) < dist)) {
                dist = connect.distTo(pixel);
                dest = pixel;
            }
        }
        final Iterable<DirectedEdge> edges = connect.pathTo(dest);
        final List<java.lang.Integer> path = new ArrayList<>();
        for (DirectedEdge edge : edges) {
            path.add(edge.to() - 1);
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

        ConnectPixels(int width, int height, Direction direction) {
            this.width = width;
            this.height = height;
            this.direction = direction;
            final int count = width * height + 1;
            this.distTo = new double[count];
            this.edgeTo = new DirectedEdge[count];
            for (int i = 0; i < count; i++) {
                distTo[i] = Double.MAX_VALUE;
                edgeTo[i] = null;
            }
            distTo[0] = 0;
            if (direction == Direction.Down) {
                for (int col = 0; col < width; col++) {
                    int id = getId(col, 0);
                    relax(-1, id);
                }
                for (int row = 0; row < height; row++) {
                    for (int col = 0; col < width; col++) {
                        int id = getId(col, row);
                        relax(id, getFirst(id));
                        relax(id, getMiddle(id));
                        relax(id, getLast(id));
                    }
                }
            } else {
                for (int row = 0; row < height; row++) {
                    int id = getId(0, row);
                    relax(-1, id);
                }

                for (int col = 0; col < width; col++) {
                    for (int row = 0; row < height; row++) {
                        int id = getId(col, row);
                        relax(id, getFirst(id));
                        relax(id, getMiddle(id));
                        relax(id, getLast(id));
                    }
                }
            }
        }

        boolean hasPathTo(int node) {
            return this.distTo[node + 1] != Double.MAX_VALUE;
        }

        double distTo(int node) {
            return this.distTo[node + 1];
        }

        Iterable<DirectedEdge> pathTo(int v) {
            Stack<DirectedEdge> path = new Stack<>();
            for (DirectedEdge e = edgeTo[v + 1]; e != null; e = edgeTo[e.from()]) {
                path.push(e);
            }
            return path;
        }

        private void relax(int from, int to) {
            if (to >= 0) {
                double energy = energyMatrix[getX(to)][getY(to)];
                double newDist = distTo[from + 1] + energy;
                if (newDist < distTo[to + 1]) {
                    distTo[to + 1] = newDist;
                    edgeTo[to + 1] = new DirectedEdge(from + 1, to + 1, energy);
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
