import java.util.Arrays;

public class Toeplitz {

    double[] vector;
    double[][] matrix;

    public Toeplitz(double[] vec) {
        this.vector = vec;
    }
    
    public Toeplitz(double[][] mat) {
        this.matrix = mat;
    }

    public double[] convertToVector(){
        int rows = this.matrix.length;
        int cols = this.matrix[0].length;

        double[] vec = new double[rows+cols-1];
        int vidx = 0;

        for (int i = rows - 1; i >= 0; i--) {
            vec[vidx++] = this.matrix[i][0];
        }
        for (int i = 1; i < cols; i++) {
            vec[vidx++] = this.matrix[0][i];
        }

        if(this.vector == null) this.vector = vec;
        return vec;
    }

    public double[][] convertToMatrix(){

        int n = this.vector.length;
        int m = (n+1)/2;

        double[][] mat = new double[m][m];

        for (int i = 0; i < n; i++) {
            int d = i - m + 1;

            int row = d < 0 ? -d : 0;
            int col = d > 0 ? d : 0;

            for (int j = 0; (row + j < m)&&(col + j < m); j++) {
                mat[row+j][col+j] = this.vector[i];
            }
        }
        if(this.matrix == null) this.matrix = mat;
        return mat;
    }

    public Toeplitz add(Toeplitz toel){

        int n = this.vector.length;
        double[] result = new double[n];

        for (int i = 0; i < n; i++) 
            result[i] = this.vector[i] + toel.vector[i];

        return new Toeplitz(result);
    }

    public double[] multiply(double[] vec){

        int n = this.matrix.length;
        double[] newVec = reverseVector(vec);
        double[] result = new double[n];

        // for (int i = 0; i < n; i++)
        //     for (int j = 0; j < n; j++)
        //         // result[i] += this.matrix[i][j] * vec[j];
        //         // result[i] += this.vector[j - i + n - 1] * vec[j];
        //         result[i] += this.vector[j - i + n - 1] * newVec[n - 1 - j];

        // double[] c = convolve(this.vector, newVec);
        double[] c = convolveFFT(this.vector, newVec);
        for (int i = 0; i < n; i++)
            result[i] = c[2*n - 2 - i];
        return result;
    }
    

    static FFT.ComplexNumber[] turnToComplex(double[] v){
        FFT.ComplexNumber[] complexVec = new FFT.ComplexNumber[v.length];
        for (int i = 0; i < v.length; i++) {
            complexVec[i] = new FFT.ComplexNumber(v[i], 0);
        }
        return complexVec;
    }

    static double[] padVector(double[] v, int minLength) {
        int N = 1;
        while (N < Math.max(v.length, minLength)) N <<= 1;

        double[] newVec = new double[N];
        System.arraycopy(v, 0, newVec, 0, v.length);
        return newVec;
    }

    static double[] convolveFFT(double[] a, double[] b){
        b = padVector(b, a.length);
        a = padVector(a, b.length);

        FFT.ComplexNumber[] c = FFT.multiplyPolys(turnToComplex(a), turnToComplex(b));
        double[] result = new double[c.length];
        for (int i = 0; i < c.length; i++)
            result[i] = c[i].getReal();
        return result;
    }

    static double[] convolve(double[] a, double[] b){
        int n = a.length + b.length - 1;
        double[] result = new double[n];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < b.length; j++)
                result[i+j] += a[i] * b[j];

        return result;
    }

    public static double[] reverseVector(double[] v) {
        double[] newVec = new double[v.length];
        for (int i = 0; i < v.length; i++) 
            newVec[i] = v[v.length - 1 - i];
        return newVec; 
    }

    public static void printVector(double[] v) {
        System.out.println(Arrays.toString(v));
    }

    public static void printMatrix(double[][] mat) {
        for (double[] row : mat) {
            for (double x : row) System.out.printf("%6.1f", x);
            System.out.println();
        }
    }

    public static void main(String[] args) {
        double[][] mat1 = {
            {3, 4, 5},
            {2, 3, 4},
            {1, 2, 3}
        };

        double[] vec = {6, 7, 8};

        double[] expected = {86, 65, 44};

        Toeplitz t1 = new Toeplitz(mat1);

        t1.convertToVector();

        double[] t3 = t1.multiply(vec);

        System.out.println("Matrix 1:");
        printMatrix(t1.matrix);
        System.out.println("Vector 2:");
        printVector(vec);
        System.out.println("Product of Matrix 1 and Vector 2:");
        printVector(t3);
        System.out.println("Expected Product:");
        printVector(expected);

    }
}
