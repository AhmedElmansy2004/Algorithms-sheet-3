public class FFT {

	public static class ComplexNumber {

		private final double real;
		private final double imag;

		public ComplexNumber(double real, double imag) {
			this.real = real;
			this.imag = imag;
		}

		public ComplexNumber mul(ComplexNumber b) {
			return new ComplexNumber(this.real * b.real - this.imag * b.imag,
					this.real * b.imag + this.imag * b.real);
		}

		public ComplexNumber add(ComplexNumber b) {
			return new ComplexNumber(this.real + b.real, this.imag + b.imag);
		}

		public ComplexNumber subtract(ComplexNumber b) {
			return new ComplexNumber(this.real - b.real, this.imag - b.imag);
		}

		public ComplexNumber scale(double s) {
			return new ComplexNumber(this.real * s, this.imag * s);
		}

		public ComplexNumber conjugate() {
			return new ComplexNumber(this.real, -this.imag);
		}

		public String toString() {
			if(this.imag != 0)
			return String.format("%.2f", this.real)
					+ (this.imag >= 0 ? " + " : " - ")
					+ String.format("%.2f", Math.abs(this.imag)) + "i";
			else
				return String.format("%.2f", this.real); 	

		}

        public double getReal() {
            return this.real;
        }

        public double getImag() {
            return this.imag;
        }
	}

	// This computes the DFT of a given coefficient vector 'a'
	// This applies the recursive-FFT procedure discussed in lecture
	public static ComplexNumber[] fft(ComplexNumber[] a) {
		
		int n = a.length;
		// base case
		if (n == 1) {
			return new ComplexNumber[] { a[0] };
		}
			
		// Checking that n is even across all recursive calls implies a
		// check on whether n is a power of two or not
		if (n % 2 != 0) {
			throw new IllegalArgumentException(
					"This implementation assumes n to be a power of 2.");
		}

		// Note: this code is not space-efficient.

		ComplexNumber[] partition = new ComplexNumber[n / 2];
		// Even partition
		for (int k = 0; k < n / 2; k++) {
			partition[k] = a[2 * k];
		}
		ComplexNumber[] l = fft(partition); 

		// Odd partition
		for (int k = 0; k < n / 2; k++) {
			partition[k] = a[2 * k + 1];
		}
		ComplexNumber[] r = fft(partition);
		
		double angle = 2 * Math.PI / n;
		ComplexNumber w_n = new ComplexNumber(Math.cos(angle), Math.sin(angle));
		ComplexNumber w = new ComplexNumber(1.0, 0);

		// combine
		ComplexNumber[] y = new ComplexNumber[n];
		for (int k = 0; k < n / 2; k++) {
			y[k] = l[k].add(w.mul(r[k]));
			y[k + n / 2] = l[k].subtract(w.mul(r[k]));
			w = w.mul(w_n);
		}
		return y;
	}

	/*
	 * This method is one possible way for implementing the inverse DFT. It uses
	 * the FFT in the forward direction as a black box. 
	 */
	public static ComplexNumber[] ifft1(ComplexNumber[] y) {

		int n = y.length;
		ComplexNumber[] yConjugate = new ComplexNumber[n];
		for (int i = 0; i < n; i++) {
			yConjugate[i] = y[i].conjugate();
		}

		// Calling DFT in the *forward* direction
		ComplexNumber[] a = fft(yConjugate);

		for (int i = 0; i < n; i++) {
			a[i] = a[i].conjugate();
		}

		for (int i = 0; i < n; i++) {
			a[i] = a[i].scale(1.0 / n);
		}
		
		return a;
	}

	
	// TODO: This method should compute the inverse DFT according to the method discussed in
	// the book/lecture. 
	public static ComplexNumber[] myIfft(ComplexNumber[] y) {
		int n = y.length;
		// base case
		if (n == 1) {
			return new ComplexNumber[] { y[0] };
		}
			
		// Checking that n is even across all recursive calls implies a
		// check on whether n is a power of two or not
		if (n % 2 != 0) {
			throw new IllegalArgumentException(
					"This implementation assumes n to be a power of 2.");
		}

		// Note: this code is not space-efficient.

		ComplexNumber[] partition = new ComplexNumber[n / 2];
		// Even partition
		for (int k = 0; k < n / 2; k++) {
			partition[k] = y[2 * k];
		}
		ComplexNumber[] l = myIfft(partition); 

		// Odd partition
		for (int k = 0; k < n / 2; k++) {
			partition[k] = y[2 * k + 1];
		}
		ComplexNumber[] r = myIfft(partition);
		
        // modified the angle
		double angle =  -2 * Math.PI / n;
		ComplexNumber w_n = new ComplexNumber(Math.cos(angle), Math.sin(angle));
		ComplexNumber w = new ComplexNumber(1.0, 0);

		// combine
		ComplexNumber[] a = new ComplexNumber[n];
		for (int k = 0; k < n / 2; k++) {
			a[k] = l[k].add(w.mul(r[k]));        
			a[k + n / 2] = l[k].subtract(w.mul(r[k]));
			w = w.mul(w_n);
		}
		return a;
	}

    public static ComplexNumber[] ifft2(ComplexNumber[] y) {

		int n = y.length;

		// Calling DFT but with (w_n)^(-k*j)
		ComplexNumber[] a = myIfft(y);

		for (int i = 0; i < n; i++) {
			a[i] = a[i].scale(1.0 / n);
		}
		
		return a;
	}

	// multiply the polynomials represented by the coefficient vectors a, b
	public static ComplexNumber[] multiplyPolys(ComplexNumber[] a,
			ComplexNumber[] b) {
		ComplexNumber zero = new ComplexNumber(0, 0);

		if (a.length != b.length) {
			throw new IllegalArgumentException("Length mismatch");
		}

		int n = a.length;
		if ((n & (n - 1)) != 0 || n <= 0) {
			throw new IllegalArgumentException("n is not a valid power of two");
		}

		// padding a
		ComplexNumber[] a2 = new ComplexNumber[2 * n];
		for (int i = 0; i < n; i++)
			a2[i] = a[i];
		for (int i = n; i < 2 * n; i++)
			a2[i] = zero;

		// padding b
		ComplexNumber[] b2 = new ComplexNumber[2 * b.length];
		for (int i = 0; i < n; i++)
			b2[i] = b[i];
		for (int i = n; i < 2 * n; i++)
			b2[i] = zero;

		// note: a2 and b2 are the padded versions of a and b
		ComplexNumber[] ya = fft(a2);
		ComplexNumber[] yb = fft(b2);

		// point-wise multiplication
		ComplexNumber[] yc = new ComplexNumber[2 * n];
		for (int i = 0; i < 2 * n; i++) {
			yc[i] = ya[i].mul(yb[i]);
		}

		// compute the inverse FFT
		return ifft2(yc); // replace with a call to ifft2 after you implement it
	}

	public static void main(String[] args) {
		int n = 4;

		ComplexNumber[] a = new ComplexNumber[n];
		ComplexNumber[] b = new ComplexNumber[n];

		// In this example, a and b are polynomials with integer coefficients

		a[0] = new ComplexNumber(3, 0);
		a[1] = new ComplexNumber(5, 0);
		a[2] = new ComplexNumber(7, 0);
		a[3] = new ComplexNumber(8, 0);
			
		b[0] = new ComplexNumber(6, 0);
		b[1] = new ComplexNumber(7, 0);
		b[2] = new ComplexNumber(1, 0);
		b[3] = new ComplexNumber(2, 0);
		
		ComplexNumber[] c = multiplyPolys(a, b);
		System.out.println("The coefficients of the output polynomial = ");
		for (int i = 0; i < c.length; i++) {
			System.out.println("c[" + i + "] = " + c[i]);
		}

	}

}