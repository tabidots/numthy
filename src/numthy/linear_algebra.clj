(ns numthy.linear-algebra
  "Matrix multiplication, Gauss-Jordan elimination, and linear least-squares in ℚ
  with exact arithmetic; Gaussian elimination and linear equation solver over GF(2)."
  (:require [clojure.core.reducers :as r]
            [clojure.core.matrix :as m]))

(m/set-current-implementation :vectorz)

(defn- transpose
  [matrix]
  (apply mapv vector matrix))

(defn- row
  [matrix idx]
  (get matrix idx))

(defn- column
  [matrix idx]
  (into [] (r/map #(get % idx) matrix)))

(defn- swap ;; m/swap-rows seems a tiny bit faster
  [matrix a b]
  (assoc matrix a (get matrix b) b (get matrix a)))

(defn- scale-row
  [row n]
  (into [] (r/map #(* % n) row)))

(defn- add-rows
  [a b]
  (into [] (pmap #(+ %1 %2) a b)))

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
    (scale-row row (/ pivot))  ;; m/scale is a bit slower on native Clojure vectors
    row))

(defn- zero-out-rest-of-column
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
        (recur (->> (scale-row pivot-row scalar) ;; m/scale is a bit slower on native Clojure vectors
                    (m/add cur-row)
                    (assoc m r)) ;; assoc seems to work much faster than m/set-row for native Clojure vectors
               (inc r))))))

(defn reduced-row-echelon
  "Applies Gauss-Jordan elimination to return the reduced-row echelon form of
  the given matrix. For large matrices, this will be very slow because of the repeated
  division operations."
  [matrix]
  (let [num-cols (count (first matrix))
        num-rows (count matrix)]
    (loop [m matrix target-row 0 target-col 0]
      ;; Finished eliminating? Return the matrix
      (if (or (= target-row num-rows)
              (= target-col num-cols)) m
        (let [pivot-row (find-pivot-row m target-col)
              target-el (m/mget m target-row target-col)]
          (cond
            ;; No pivot row -> Column is all zeroes? Go to next column, but same row
            (nil? pivot-row) (recur m target-row (inc target-col))
            ;; Just this element is zero? Swap with next row that has a pivot in column at idx
            (zero? target-el) (recur (m/swap-rows m target-row pivot-row) target-row target-col)
            ;; Has a pivot? Make that pivot 1 and zero out the remaining columns;
            ;; Proceed to next column and row
            :else
            (recur (->> (row m pivot-row)
                        (make-row-have-pivot-1)
                        (assoc m pivot-row)     ;; assoc seems to work much faster than m/set-row for native Clojure vectors
                        (zero-out-rest-of-column pivot-row target-col))
                   (inc target-row)
                   (inc target-col))))))))

(defn mmul
  "Dot product of A and b, where b can be a matrix or vector."
  [A b]
  ;; Rows on right have to equal cols on left
  (when (= (count (first A)) (count b))
    (if (get-in b [0 0])
      ;; b is matrix
      (reduce (fn [res left-row]
                (conj res (mapv (fn [right-col]
                                  (reduce + (mapv * left-row right-col)))
                                (transpose b)))) ;; right-cols = (transpose b)
              [] A)
      ;; b is vector
      (mapv (fn [left-row]
              (reduce + (mapv * left-row b)))
            A))))

(defn- augment
  "Augmented matrix of A and b, where A is a matrix and b is a vector."
  [A b]
  (mapv conj A b))

(defn least-squares
  "Least-squares solution of Ax = b, where A is a matrix and b is a vector."
  [A b]
  ;; With reference to https://textbooks.math.gatech.edu/ila/least-squares.html
  (let [at   (transpose A)
        ata  (mmul at A)      ;; m/mmul is slow on Clojure vectors
        atb  (mmul at b)      ;; m/mmul is slow on Clojure vectors
        aug  (augment ata atb)   ;; (augment ata atb) roughly same performance
        _    (println aug)
        rref (reduced-row-echelon aug)]
    (peek (m/columns rref))))    ;; (mapv peek rref) roughly same performance

;; GF(2)-specific functions

(defn- xor-rows
  "XOR an origin row of a matrix in GF(2) to a destination row. If a column index
  is supplied, then XOR will be applied only from that column rightwards, and the
  old values to the left of that column will be kept. If no column index is supplied,
  the whole row will be XOR'ed."
  ([origin destination]
   (xor-rows 0 origin destination))
  ([from-column origin destination]
   (vec (concat (take from-column destination)
                (drop from-column (mapv #(bit-xor %1 %2) origin destination))))))

(defn- zero-out-rest-of-column-gf2
  "Given a pivot row and column, XORs the pivot row with all other rows having a
  1 in that column, for all values in that column and rightwards."
  [row-idx col-idx matrix]
  (loop [m matrix
         r 0]
    (cond
      (= r (count m)) m                ;; Finished all rows
      (or (= r row-idx)
          (zero? (m/mget m r col-idx))) (recur m (inc r))  ;; Ignore pivot & zero entries
      :else (recur (->> (row m r)
                        (xor-rows col-idx (row m row-idx))
                        (assoc m r)) ;; assoc seems to work much faster than m/set-row for native Clojure vectors
                   (inc r)))))

(defn rref-gf2
  "Uses Gaussian elimination to return the row-reduced echelon form of a matrix in GF(2)."
  [matrix]
  (let [num-cols (count (first matrix))
        num-rows (count matrix)]
    (loop [m matrix target-row 0 target-col 0]
      ;; Finished eliminating? Return the matrix
      (if (or (= target-row num-rows)
              (= target-col num-cols)) m
        (let [pivot-row (find-pivot-row m target-col)
              target-el (m/mget m target-row target-col)]
          (cond
            ;; No pivot row -> Column is all zeroes? Go to next column, but same row
            (nil? pivot-row) (recur m target-row (inc target-col))
            ;; Just this element is zero? Swap with next row that has a pivot in column at idx
            (zero? target-el) (recur (m/swap-rows m target-row pivot-row) target-row target-col)
            ;; Has a pivot? Make that pivot 1 and zero out the remaining columns;
            ;; Proceed to next column and row
            :else
            (recur (zero-out-rest-of-column-gf2 pivot-row target-col m)
                   (inc target-row)
                   (inc target-col))))))))

(comment
  "Keep these as test cases"
  (def A [[0 0 1 0 1 1 0 0 1 0]
          [0 1 0 1 0 0 0 0 1 0]
          [0 0 0 0 1 0 0 0 1 0]
          [0 0 1 0 1 1 0 0 0 0]
          [0 1 1 1 1 0 0 0 0 0]])
  (def b [0,1,0,1,0,0,0,0,0,0]))

(defn solve-gf2
  "Solves Ax = b, where A is a matrix and b is a vector in GF(2)."
  [A b]
  (when (= (count A) (count b))
    (let [num-entries (count (first A))]
      (->> (augment A b)
           (rref-gf2)
           (take num-entries) ;; because the result is zero-padded
           (mapv peek)))))
