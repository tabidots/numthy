(ns numthy.factorization.squares-utils
  "Workhorse functions for the family of congruence-of-squares factorization algorithms."
  (:require [clojure.math.numeric-tower :refer [gcd expt abs]]
            [clojure.core.reducers :as r]
            [clojure.core.matrix :as m]
            [clojure.string :as s]
            [numthy.modular-arithmetic.utils :refer [mod-mul]]))

(defn smoothness-bound
  "Ceiling of e ^ (1/2 sqrt(ln n ln ln n))."
  ;; https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html
  ;; and some papers as well
  [n]
  (-> (*' (Math/log n) (Math/log (Math/log n)))
      Math/sqrt (* 1/2) Math/exp))

(defn smooth?
  "Returns [exponent-vector] if n is smooth over the given empty hash-map of primes."
  [n primes prime-map]  ;; extra arg because iterating over keySeq is slow
  (loop [x      (abs n)
         m      prime-map
         primes primes]
    (if (= 1 x) (vec (vals m))
      (when-let [k (first primes)]
        (cond
          (neg? k)           (if (neg? n)
                               (recur x (update m -1 inc) (rest primes))
                               (recur x m (rest primes)))
          (zero? (mod x k)) (recur (/ x k) (update m k inc) primes)
          :else             (recur x m (rest primes)))))))

(defn row->mod2
  "Given a vector, reduces all of its entries mod 2."
  [row]
  (mapv #(mod % 2) row))

(defn gf2-matrix->column-wise-binary
  "Encodes each column of a matrix in GF(2) as a binary number. Uses BigIntegers
  to accommodate matrices of arbitrary numbers of columns and rows."
  [m]
  (apply map (comp #(BigInteger. % 2)
                   #(s/reverse %)
                   #(reduce str %)
                   vector)
             m))

(comment
 "Keep as test data"
  (def m [[1 1 0 0]
          [1 1 0 1]
          [0 1 1 1]
          [0 0 1 0]
          [0 0 0 1]]))

(defn- lowest-set-bit
  "Finds the lowest set bit of a binary number. Only works for native Clojure numbers.
  This function is left here for reference."
  [a]
  (when (pos? a)
    (loop [b (.and a (- a)) low-bit -1]
      (if (zero? b) low-bit
        (recur (bit-shift-right b 1) (inc low-bit))))))

(defn get-pivot
  "Finds the pivot (first '1') in a binary-encoded column of a GF(2) matrix. Since
  the columns are encoded so that the top entry of the column is the rightmost bit,
  the first 1 is the lowest set bit of the number. When the pivot is -1, there is no
  pivot, so returns nil. Uses Java interop to handle BigIntegers."
  [n]
  (let [lb (.getLowestSetBit n)]
    (when-not (neg? lb)
      (long lb))))

(defn find-winning-subsets
  "Uses Koç & Arachchige's algorithm to perform fast Gaussian elimination in GF(2) and
  find all subsets of rows that add to the zero vector."
  ;; With reference to https://github.com/skollmann/PyFactorise/blob/master/factorise.py
  ;; Original paper: https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
  [m]
  (let [num-cols          (count (first m))
        num-rows          (count m)
        columns           (to-array (gf2-matrix->column-wise-binary m))
        ^booleans marked? (boolean-array num-rows false)
        ^longs    pivots  (long-array num-cols -1)]
    (dotimes [j num-cols]
             (when-some [i (get-pivot (aget columns j))]
               (aset-long pivots j i)
               (aset-boolean marked? i true)
               (dotimes [k num-cols]
                        (when (and (not= k j) (.testBit (aget columns k) i))
                          (aset columns k (.xor (aget columns k) (aget columns j)))))))
    (for [i (range num-rows) :when (not (aget ^booleans marked? i))]
      (cons i
            (for [j (range num-cols) :when (.testBit (aget columns j) i)]
              (aget ^longs pivots j))))))

(defn make-x-from-idxs
  "Given a list of indices, extracts the corresponding z-values from the candidates
  and multiplies them together to give x."
  [idxs zs n]
  (->> (m/order zs idxs)
       (r/fold (fn ([] 1) ([a b] (mod-mul a b n))))))

(defn make-y-from-idxs
  "Given a list of indices, extracts the corresponding rows from the original
  exponent matrix and adds their values (which is equivalent to multiplying the
  integers they represent). This gives y^2. Then finds y by halving the entries
  of the sum vector, and converting that back to an integer."
  [idxs expm factor-base n]
  (->> (m/order expm idxs)
       (apply m/add)
       (zipmap factor-base)
       (r/fold (r/monoid (fn [a b] (mod-mul a b n)) (constantly 1))
               (fn ([res b e]
                    (-> (expt b (quot e 2))
                        (mod-mul res n)))))))

(defn find-factor
  "Separates a matrix of relations into linearly independent and dependent rows.
  For the submatrix of linearly independent rows A, and each dependent row b,
  solves the matrix equation Ax = b in GF(2) to find the vectors of the null space of A
  until a congruence of squares is found."
  [n factor-base relations]
  (let [xs               (vec (keys relations))
        exponent-vectors (vec (vals relations))
        matrix           (mapv row->mod2 exponent-vectors)]
    (->> (for [subset (find-winning-subsets matrix)]
           (let [x (make-x-from-idxs subset xs n)
                 y (make-y-from-idxs subset exponent-vectors factor-base n)]
             (gcd (- x y) n)))
         (filter #(< 1 % n))
         (first))))
