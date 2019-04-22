(ns numthy.linear-algebra)

(defn transpose
  [matrix]
  (apply mapv vector matrix))

(defn- row
  [matrix idx]
  (get matrix idx))

(defn- column
  [matrix idx]
  (mapv #(get % idx) matrix))

(defn- swap
  [matrix a b]
  (assoc matrix a (get matrix b) b (get matrix a)))

(defn- scale-row
  [row n]
  (mapv #(*' % n) row))

(defn- add-rows
  [a b]
  (mapv + a b))

(defn- find-pivot-row
  "The index of the first row whose pivot is in the given column."
  [matrix col]
  (letfn [(pivot-index [row]
            (first (keep-indexed (fn [i c]
                                   (when-not (zero? c) i))
                                 row)))]
    (first (keep-indexed (fn [i row]
                           (when (= (pivot-index row) col) i))
                         matrix))))

(defn- make-row-have-pivot-1
  "Scale a row by the inverse of its pivot (its first non-zero value),
  so that its pivot is 1."
  [row]
  (if-let [pivot (some #(when-not (zero? %) %) row)]
    (scale-row row (/ pivot))
    row))

(defn- zero-out-remaining-columns
  "Applies add & scale row operations so that all entries in column #col-idx,
  apart from the one in row #row-idx are zero."
  [row-idx col-idx matrix]
  (loop [m matrix
         r 0]
    (condp = r
      (count m)       m                   ;; Finished all rows
      row-idx         (recur m (inc r))   ;; Ignore the pivot row
      (let [cur-row   (row m r)
            pivot-row (row m row-idx)
            top       (- (get cur-row col-idx))
            bot       (get pivot-row col-idx)
            scalar    (/ top bot)]
        (recur (assoc m r (add-rows cur-row (scale-row pivot-row scalar)))
               (inc r))))))

(defn reduced-row-echelon
  "Applies Gauss-Jordan elimination to returns the reduced-row echelon form
  of the given matrix."
  [matrix]
  (let [num-cols (count (first matrix))
        num-rows (count matrix)]
    (loop [m matrix target-row 0 target-col 0]
      (let [pivot-row (find-pivot-row m target-col)
            target-el (get-in m [target-row target-col])]
        (cond
          ;; Finished eliminating? Return the matrix
          (= target-col num-cols) m
          ;; No pivot row -> Column is all zeroes? Go to next column, but same row
          (nil? pivot-row) (recur m target-row (inc target-col))
          ;; Just this element is zero? Swap with next row that has a pivot in column at idx
          (zero? target-el) (recur (swap m target-row pivot-row) target-row target-col)
          ;; Has a pivot? Make that pivot 1 and zero out the remaining columns;
          ;; Proceed to next column and row
          :else
          (recur (->> (assoc m pivot-row (make-row-have-pivot-1 (row m pivot-row)))
                      (zero-out-remaining-columns pivot-row target-col))
                 (inc target-row)
                 (inc target-col)))))))

(defn mmul
  "Dot product of a and b, where b can be a matrix or vector."
  [a b]
  ;; Rows on right have to equal cols on left
  (when (= (count (first a)) (count b))
    (if (get-in b [0 0])
      ;; b is matrix
      (reduce (fn [res left-row]
                (conj res (mapv (fn [right-col]
                                  (reduce + (mapv * left-row right-col)))
                                (transpose b)))) ;; right-cols = (transpose b)
              [] a)
      ;; b is vector
      (mapv (fn [left-row]
              (reduce + (mapv * left-row b)))
            a))))

(defn- augment
  "Augmented matrix of a and b, where b is a vector."
  [a b]
  (mapv conj a b))

(defn least-squares
  "Least-squares solution of Ax = b, where A is a matrix and b is a vector."
  [a b]
  ;; With reference to https://textbooks.math.gatech.edu/ila/least-squares.html
  (let [ata  (mmul (transpose a) a)
        atb  (mmul (transpose a) b)
        aug  (augment ata atb)
        rref (reduced-row-echelon aug)]
    (mapv last rref)))
