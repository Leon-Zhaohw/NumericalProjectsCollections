

/**
 * Class for efficient manipulation of sparse matrices. A matrix is stored as
 * non-zero only values, the rows are stored in jagged two-dimensional arrays as
 * suggested by Geir Gundersen and Trond Steihaug in "Data structures in Java
 * for Matrix computations" Each row is stored in a single-dimensional array
 * (that grows as necessary), the column indexes are stored accordingly. <p/>
 * For example, the matrix
 * <pre>
 * 0 1 0 0
 * 2 0 3 0
 * 0 0 0 0
 * 4 1 0 0
 * </pre>
 * is stored as
 * <pre>
 * [1] [1]
 * [2, 3] [0, 2]
 * [] []
 * [4, 1] [0, 1]
 * </pre>
 *
 * @author Kirill Grouchnikov
 */


public class SparseMatrix {
    protected double[][] nzValues;

    protected int[][] columnIndices;
    protected int[] nzCounters;
    public int colCount=0;
    public int rowCount=0;
    public static void main (String[] Argv) {
	SparseMatrix a = new SparseMatrix(4,4);
	
	double []v = new double[4];

	a.set(3,0,1.0);
	a.set(3,1,-2.0);
	a.set(3,2,3.0);
	a.set(3,3,-4.0);

	v[0] = 10.0;
	v[1] = 10.0;
	v[2] = 10.0;
	v[3] = 10.0;

	
	a.dump();

	double []x = a.mult(v);
	for (int i=0;i<4;i++) 
	    System.out.printf("%f\n",x[i]);

    }


    /**
     * @param colCount number of columns
     * @param rowCount number of rows
     */

    public SparseMatrix(int colCount, int rowCount) {
	this.rowCount=rowCount;
	this.colCount=colCount;
        this.nzValues = new double[rowCount][];

        this.columnIndices = new int[rowCount][];

        this.nzCounters = new int[rowCount];

    }



    /**
     * Gets values at specified location
     *
     * @param column column index
     * @param row    row index
     * @return value
     */
    public double get(int column, int row) {
        if (this.columnIndices[row] == null) {
            return 0.0;
        }
        int columnIndex = binarySearch(this.columnIndices[row], 0,
				       this.nzCounters[row] - 1, column);

        if (columnIndex < 0) {
            return 0.0;
        }

        return this.nzValues[row][columnIndex];

    }



    /**

     * Performs a binary search for a given value in sorted integer array. The

     * only difference from <b>Arrays.binarySearch()</b> is that this function

     * gets <i>start</i> and <i>end</i> indexes. The array <strong>must</strong>

     * be sorted prior to making this call.  If it is not sorted, the results

     * are undefined.

     *

     * @param array      array to scan

     * @param startIndex start index of sub-array

     * @param endIndex   end index of sub-array

     * @param value      key to find

     * @return index of the search key, if it is contained in the list;

     *         otherwise, <tt>(-(<i>insertion point</i>) - 1)</tt>.  The

     *         <i>insertion point</i> is defined as the point at which the key

     *         would be inserted into the list: the index of the first element

     *         greater than the key, or <tt>list.size()</tt>, if all elements in

     *         the list are less than the specified key.  Note that this

     *         guarantees that the return value will be &gt;= 0 if and only if

     *         the key is found.

     */

    private static int binarySearch(int[] array, int startIndex, int endIndex,

                                    int value) {

        if (value < array[startIndex]) {

            return (-startIndex - 1);

        }

        if (value > array[endIndex]) {

            return (-(endIndex + 1) - 1);

        }



        if (startIndex == endIndex) {

            if (array[startIndex] == value) {

                return startIndex;

            }

            else {

                return (-(startIndex + 1) - 1);

            }

        }

        int midIndex = (startIndex + endIndex) / 2;

        if (value == array[midIndex]) {

            return midIndex;

        }



        if (value < array[midIndex]) {

            return binarySearch(array, startIndex, midIndex - 1, value);

        }

        else {

            return binarySearch(array, midIndex + 1, endIndex, value);

        }

    }




    /**

     * Sets value at specified location

     *

     * @param column column index

     * @param row    row index

     * @param value  value

     */

    public void set(int column, int row, double value) {

        if (this.columnIndices[row] == null) {

            // first value in this row

            this.columnIndices[row] = new int[2];

            this.nzValues[row] = new double[2];

            this.columnIndices[row][0] = column;

            this.nzValues[row][0] = value;

            this.nzCounters[row] = 1;

            return;

        }



        // search for it

        int columnIndex = binarySearch(this.columnIndices[row], 0,

				       this.nzCounters[row] - 1, column);

        if (columnIndex >= 0) {

            // already setLocation, just change

            this.nzValues[row][columnIndex] = value;

            return;

        }

        else {

            // columnIndex = (-(insertion point) - 1)

            int insertionPoint = -(columnIndex + 1);

            // allocate new arrays

            int oldLength = this.nzCounters[row];

            int newLength = oldLength + 1;

            // check if need to allocate

            if (newLength <= this.columnIndices[row].length) {

                // just copy

                if (insertionPoint != oldLength) {

                    for (int i = oldLength; i > insertionPoint; i--) {

                        this.nzValues[row][i] = this.nzValues[row][i - 1];

                        this.columnIndices[row][i] =

			    this.columnIndices[row][i - 1];

                    }

                }

                this.columnIndices[row][insertionPoint] = column;

                this.nzValues[row][insertionPoint] = value;

                this.nzCounters[row]++;

                return;

            }



            int[] newColumnIndices = new int[2 * oldLength];

            double[] newNzValues = new double[2 * oldLength];



            if (insertionPoint == oldLength) {

                // special case - new column is the last

                System.arraycopy(this.columnIndices[row], 0, newColumnIndices,

				 0, oldLength);

                System.arraycopy(this.nzValues[row], 0, newNzValues, 0,

				 oldLength);

            }

            else {

                System.arraycopy(this.columnIndices[row], 0, newColumnIndices,

				 0, insertionPoint);

                System.arraycopy(this.nzValues[row], 0, newNzValues, 0,

				 insertionPoint);

                System.arraycopy(this.columnIndices[row], insertionPoint,

				 newColumnIndices, insertionPoint + 1, oldLength

				 - insertionPoint);

                System.arraycopy(this.nzValues[row], insertionPoint,

				 newNzValues, insertionPoint + 1, oldLength

				 - insertionPoint);

            }

            newColumnIndices[insertionPoint] = column;

            newNzValues[insertionPoint] = value;

            this.columnIndices[row] = null;

            this.columnIndices[row] = newColumnIndices;

            this.nzValues[row] = null;

            this.nzValues[row] = newNzValues;

            this.nzCounters[row]++;



        }

    }





    /**

     * Dump to standard output

     */

    public void dump() {

        System.out.println("MATRIX " + this.rowCount + "*" + this.colCount);

        for (int row = 0; row < this.rowCount; row++) {

            int[] columnIndices = this.columnIndices[row];

            if (columnIndices == null) {

                for (int col = 0; col < this.colCount; col++) {

                    System.out.print("0.0 ");

                }

            }

            else {

                int prevColumnIndex = 0;

                for (int colIndex = 0;

                     colIndex < this.nzCounters[row]; colIndex++) {

                    int currColumnIndex = columnIndices[colIndex];

                    // put zeroes

                    for (int col = prevColumnIndex;

                         col < currColumnIndex; col++) {

                        System.out.print("0.0 ");

                    }

                    System.out.print(this.nzValues[row][colIndex] + " ");

                    prevColumnIndex = currColumnIndex + 1;

                }

                // put trailing zeroes

                for (int col = prevColumnIndex; col < this.colCount; col++) {

                    System.out.print("0.0 ");

                }

            }



            System.out.println();

        }

    }



    /**

     * Dump to standard output as integer values

     */

    public void dumpInt() {

        System.out.println("MATRIX " + this.rowCount + "*" + this.colCount);

        for (int row = 0; row < this.rowCount; row++) {

            int[] columnIndices = this.columnIndices[row];

            if (columnIndices == null) {

                for (int col = 0; col < this.colCount; col++) {

                    System.out.print("0 ");

                }

            }

            else {

                int prevColumnIndex = 0;

                for (int colIndex = 0;

                     colIndex < this.nzCounters[row]; colIndex++) {

                    int currColumnIndex = columnIndices[colIndex];

                    // put zeroes

                    for (int col = prevColumnIndex;

                         col < currColumnIndex; col++) {

                        System.out.print("0 ");

                    }

                    System.out.print((int) this.nzValues[row][colIndex] + " ");

                    prevColumnIndex = currColumnIndex + 1;

                }

                // put trailing zeroes

                for (int col = prevColumnIndex; col < this.colCount; col++) {

                    System.out.print("0 ");

                }

            }



            System.out.println();

        }

    }



    /**

     * Add empty (zero) columns to this matrix

     *

     * @param columns number of columns to add

     */

    public void addEmptyColumns(int columns) {

        // just as easy as that

        this.colCount += columns;

    }


    /**

     * Multiply this matrix by the specified vector

     *

     * @param vector vector

     * @return vector result

     */

    public double[] mult(double[] vector) {

        if (this.colCount != vector.length) {
            return null;
        }


        int n = this.rowCount;

        double[] result = new double[n];

        for (int row = 0; row < n; row++) {

            double sum = 0.0;

            // go over all non-zero column of this row

            int[] nzIndexes = this.columnIndices[row];

            int nzLength = nzCounters[row];

            if (nzLength == 0) {

                continue;

            }

            for (int colIndex = 0; colIndex < nzLength; colIndex++) {

                double c = vector[nzIndexes[colIndex]];

                sum += (this.nzValues[row][colIndex] * c);

            }

            result[row] = sum;

        }

        return result;

    }



    /**

     * Compute the sum of all elements in given row

     *

     * @param row row index

     * @return sum of all elements in this row

     */

    public double getSum(int row) {

        double sum = 0.0;

        // go over all non-zero column of this row

        int nzLength = nzCounters[row];

        if (nzLength == 0) {

            return 0.0;

        }

        for (int colIndex = 0; colIndex < nzLength; colIndex++) {

            sum += this.nzValues[row][colIndex];

        }

        return sum;

    }



    /**

     * Compute the sum of all elements in given row. Only entries with 'true' in

     * the corresponding position of <b>toConsider</b> array are summed up.

     *

     * @param row        row index

     * @param toConsider boolean array

     * @return sum of all elements in this row

     */

    public double getSum(int row, boolean[] toConsider) {

        double sum = 0.0;

        // go over all non-zero column of this row

        int nzLength = nzCounters[row];

        if (nzLength == 0) {

            return 0.0;

        }

        for (int colIndex = 0; colIndex < nzLength; colIndex++) {

            if (toConsider[this.columnIndices[row][colIndex]]) {

                sum += this.nzValues[row][colIndex];

            }

        }

        return sum;

    }





    public int getNzCount() {

        int allNz = 0;

        for (int i : this.nzCounters) {

            allNz += i;

        }

        return allNz;

    }



    /**

     * Normalize all values to be in 0.0-1.0 range

     */

    public void normalize() {

        double minValue = 0.0;

        double maxValue = 0.0;

        boolean isFirst = true;

        for (int i = 0; i < this.rowCount; i++) {

            if (this.nzCounters[i] == 0) {

                continue;

            }

            for (int j = 0; j < this.nzCounters[i]; j++) {

                double val = this.nzValues[i][j];

                if (isFirst) {

                    minValue = val;

                    maxValue = val;

                    isFirst = false;

                }

                else {

                    if (val < minValue) {

                        minValue = val;

                    }

                    if (val > maxValue) {

                        maxValue = val;

                    }

                }

            }

        }



        for (int i = 0; i < this.rowCount; i++) {

            if (this.nzCounters[i] == 0) {

                continue;

            }

            for (int j = 0; j < this.nzCounters[i]; j++) {

                this.nzValues[i][j] /= maxValue;

            }

        }

    }

}
