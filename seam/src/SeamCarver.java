import edu.princeton.cs.algs4.Picture;

import java.awt.*;

public class SeamCarver {
    private Picture origin;
    private Picture current;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
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

    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {

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
}
